package org.reactome.orthoinference;

import org.gk.model.ReactomeJavaConstants;
import org.reactome.curation.controller.CurationController;
import org.reactome.curation.model.CurationAttribute;
import org.reactome.curation.model.InstanceList;
import org.reactome.curation.model.NamedReferrerList;
import org.reactome.curation.model.SimpleInstance;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import javax.servlet.http.HttpServletRequest;
import java.lang.reflect.Proxy;
import java.util.*;
import java.util.concurrent.atomic.AtomicLong;
import java.util.stream.Collectors;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 * Created 7/5/2026
 */
public class CuratorToolAPI {

    private static final Logger logger = LoggerFactory.getLogger(CuratorToolAPI.class);

    private static final int PAGE_SIZE = 500;

    private final CurationController controller;
    private final long personId;

    // The GO instances are loaded once, up-front, and kept in memory for the whole run. The server rejects a
    // commit whose "modified" InstanceEdit does not match the stored one (optimistic locking) and a commit
    // replaces *all* of the stored instance's attributes, so an in-memory copy must be re-read once this run has
    // written to it. These sets record what this run has written to and what it has removed.
    private final Set<Long> committedDbIds = new HashSet<>();
    private final Set<Long> deletedDbIds = new HashSet<>();

    // Source of the negative placeholder dbIds given to instances that are not in the database yet.
    private final AtomicLong placeholderDbIdCounter = new AtomicLong();

    // A species instance is wanted once per SwissProt entry, from the same dozen species names for the whole run, so
    // the instance found (or created) for a name is kept rather than queried again for every entry of that species.
    private final Map<String, SimpleInstance> speciesNameToInstance = new HashMap<>();

    // The names of each schema class's multi-valued attributes, looked up once per class. Every instance sent to
    // curator-tool-ws is checked against them, so they are not fetched again for each instance.
    private final Map<String, Set<String>> classNameToMultiValuedAttributeNames = new HashMap<>();

    // Number of calls to, and time spent in, each kind of curator-tool-ws call, keyed by a label for the call. Kept to
    // find where the run spends its time; see logCallStats.
    private final Map<String, CallStats> callStats = new LinkedHashMap<>();
    // The dbIds inflated so far, so that inflating one again can be counted separately ("inflate (repeat)").
    private final Set<Long> inflatedDbIds = new HashSet<>();

    // Inflated instances, by dbId, as they were read from the database; inflate hands out copies of them. The run
    // reads the same human instances again and again (an EWAS is inflated each time a reaction it takes part in is
    // inferred) while changing few of them, so most inflates can be answered from here. An entry is dropped when a
    // commit may have changed what reading the instance would return -- see forgetInstancesChangedBy -- and the
    // whole cache is dropped on a delete. The least recently used entries are dropped beyond the limit.
    private static final int INFLATED_INSTANCE_CACHE_LIMIT = 200_000;
    private final Map<Long, SimpleInstance> inflatedInstanceCache =
        new LinkedHashMap<Long, SimpleInstance>(16, 0.75f, true) {
            @Override
            protected boolean removeEldestEntry(Map.Entry<Long, SimpleInstance> eldest) {
                return size() > INFLATED_INSTANCE_CACHE_LIMIT;
            }
        };
    private final long startNanos = System.nanoTime();

    // Built by CuratorToolConfig, which takes the controller from the curator-tool-ws context it starts and closes.
    public CuratorToolAPI(CurationController controller, ConfigProperties configProperties) {
        this.controller = controller;
        this.personId = configProperties.getPersonId();
    }

    public SimpleInstance commit(SimpleInstance simpleInstance) {
        if (simpleInstance.getDefaultPersonId() == null) {
            simpleInstance.setDefaultPersonId(getPersonId());
        }

        boolean isNewInstance = simpleInstance.getDbId() == null || simpleInstance.getDbId() < 0;
        if (simpleInstance.getDbId() == null) {
            simpleInstance.setDbId(nextPlaceholderDbId());
        }

        prepareAttributeValues(simpleInstance);
        // Collected before the commit: it replaces the placeholder dbIds of new instances, and the cached copy of the
        // instance names what it referred to before.
        Set<Long> possiblyChangedDbIds = getDbIdsPossiblyChangedByCommit(simpleInstance);
        SimpleInstance committedInstance;
        try {
            // Orthoinference names each inferred instance's StableIdentifier after its human source
            // (R-MMU-<human dbId>), so curator-tool-ws is told to keep it rather than rename it after the inferred
            // instance's own dbId.
            committedInstance = timed("commit", () -> controller.commit(simpleInstance, true));
        } finally {
            // Also when the commit fails, since it may have written part of what it was given.
            possiblyChangedDbIds.forEach(inflatedInstanceCache::remove);
        }
        // Nothing was stored, so callers cannot go on to use the instance -- for a new one, the dbId it would have
        // been given is the very thing that is missing. Reported here rather than as a NullPointerException at
        // whichever caller dereferences the result first.
        if (committedInstance == null) {
            throw new IllegalStateException("Commit of " + simpleInstance + " returned no instance");
        }

        if (isNewInstance) {
            // The response carries the dbId the database assigned in place of the placeholder, and it is the only
            // place it is reported, so it is copied back onto the instance the caller holds.
            simpleInstance.setDbId(committedInstance.getDbId());
        }

        Long dbId = simpleInstance.getDbId();
        if (dbId != null && dbId > 0) {
            committedDbIds.add(dbId);
            deletedDbIds.remove(dbId);
        }
        return committedInstance;
    }

    /**
     * Returns the dbIds of the stored instances whose inflated form a commit of the instance may change: the
     * instance itself and every stored instance it refers to, both as given and as last read. Committing an instance
     * rewrites its relationships, and an instance at the other end of one can see it from its own side -- e.g. an
     * inferred instance's inferredFrom is its source's inferredTo, so committing the former changes the latter.
     * The same goes for the new instances committed along with it (e.g. a new ModifiedResidue in a new EWAS).
     */
    private Set<Long> getDbIdsPossiblyChangedByCommit(SimpleInstance instance) {
        Set<Long> dbIds = new HashSet<>();
        collectDbIdsPossiblyChangedByCommit(instance, dbIds, Collections.newSetFromMap(new IdentityHashMap<>()));
        SimpleInstance cachedInstance = instance.getDbId() != null ? inflatedInstanceCache.get(instance.getDbId()) : null;
        if (cachedInstance != null) {
            collectDbIdsPossiblyChangedByCommit(cachedInstance, dbIds, Collections.newSetFromMap(new IdentityHashMap<>()));
        }
        return dbIds;
    }

    private void collectDbIdsPossiblyChangedByCommit(SimpleInstance instance, Set<Long> dbIds,
                                                     Set<SimpleInstance> visitedInstances) {
        if (!visitedInstances.add(instance)) {
            return;
        }
        Long dbId = instance.getDbId();
        if (dbId != null && dbId > 0) {
            dbIds.add(dbId);
        }
        if (instance.getAttributes() == null) {
            return;
        }
        for (Object value : instance.getAttributes().values()) {
            Collection<?> values = value instanceof Collection ? (Collection<?>) value : Collections.singletonList(value);
            for (Object element : values) {
                if (!(element instanceof SimpleInstance)) {
                    continue;
                }
                SimpleInstance referredInstance = (SimpleInstance) element;
                Long referredDbId = referredInstance.getDbId();
                if (referredDbId == null || referredDbId < 0) {
                    // A new instance, committed along with this one.
                    collectDbIdsPossiblyChangedByCommit(referredInstance, dbIds, visitedInstances);
                } else {
                    dbIds.add(referredDbId);
                }
            }
        }
    }

    /**
     * Returns a placeholder dbId for an instance that is not in the database yet. curator-tool-ws identifies a new
     * instance by a NEGATIVE dbId: it is what makes a commit store the instance (recording the author in its
     * "created" slot rather than in "modified") and replace the placeholder with a real dbId. A null dbId is not a
     * substitute -- curator-tool-ws unboxes the dbId without a null check while working out what to store
     * (DatabaseObjectInstanceConverter.convert and CurationService.grepNewInstances), so committing an instance
     * with a null dbId fails with a NullPointerException.
     *
     * @return a negative dbId, unused by any other instance created during this run.
     */
    private long nextPlaceholderDbId() {
        return -placeholderDbIdCounter.incrementAndGet();
    }

    /**
     * Returns a copy of the instance that is up to date with the database if this run has already committed it,
     * and the instance itself otherwise. Returns null if this run has deleted the instance, in which case the
     * caller must stop using it -- committing it would re-create the deleted instance.
     *
     * An instance that this run has committed is stale in two ways: its "modified" InstanceEdit no longer matches
     * the stored one, so a further commit is rejected with an InstanceChangedException, and any attribute written
     * by that commit (or by a commit of a fresh copy of the same instance elsewhere in the run) still holds its
     * pre-commit value, which the next commit would write back over the stored value.
     *
     * @param instance - the instance to refresh.
     * @return an up-to-date instance, or null if the instance has been deleted by this run.
     */
    public SimpleInstance refresh(SimpleInstance instance) {
        Long dbId = instance.getDbId();
        if (dbId == null) {
            return instance; // Never committed, so there is nothing stored to be out of date with.
        }
        if (deletedDbIds.contains(dbId)) {
            return null;
        }
        return committedDbIds.contains(dbId) ? inflate(instance) : instance;
    }

    public SimpleInstance findByDbId(long dbId) {
        // findByDdIdInInstance reads the instance's attributes straight off its node. Loading it as a DatabaseObject
        // instead maps it -- and every instance it refers to -- through graph-core's domain model, which fails on
        // any node whose stored values do not match the model's types (e.g. a GO term with a single string name,
        // where the model expects a list).
        if (!timed("existsByDbId", () -> controller.existsByDbId(dbId))) {
            return null;
        }
        return timed("findByDbId", () -> controller.findByDdIdInInstance(dbId));
    }

    public SimpleInstance findByDisplayName(String className, String displayName) {
        // NB: the controller takes the display name first and a comma-separated list of class names second.
        return timed("findByDisplayName", () -> controller.findByDisplayName(displayName, className));
    }

    public void deleteInstance(SimpleInstance instance) {
        if (instance.getDefaultPersonId() == null) {
            instance.setDefaultPersonId(getPersonId());
        }

        try {
            timed("delete", () -> controller.delete(instance));
        } finally {
            // A delete also changes every instance that referred to the deleted one, which nothing here keeps track
            // of. Orthoinference deletes rarely, so the whole cache is dropped.
            inflatedInstanceCache.clear();
        }
        if (instance.getDbId() != null) {
            deletedDbIds.add(instance.getDbId());
            committedDbIds.remove(instance.getDbId());
        }
    }

    public void deleteByDbId(long noReferrerDbId) {
        SimpleInstance instance = timed("findByDbId", () -> controller.findByDdIdInInstance(noReferrerDbId));
        // Nothing stored under that db id -- there is nothing to delete, and passing the null on would fail inside
        // curator-tool-ws instead of saying which db id was asked for.
        if (instance == null) {
            logger.warn("No instance found for db id " + noReferrerDbId + " -- nothing to delete");
            return;
        }
        deleteInstance(instance);
    }

    public SimpleInstance inflate(SimpleInstance shellInstance) {
        Long dbId = shellInstance.getDbId();
        SimpleInstance cachedInstance = inflatedInstanceCache.get(dbId);
        if (cachedInstance != null) {
            return timed("inflate (cached)", () -> copyOf(cachedInstance));
        }
        String label = inflatedDbIds.add(dbId) ? "inflate (first)" : "inflate (repeat)";
        SimpleInstance inflatedInstance = timed(label, () -> controller.findByDdIdInInstance(dbId));
        if (inflatedInstance != null && dbId != null && dbId > 0) {
            inflatedInstanceCache.put(dbId, copyOf(inflatedInstance));
        }
        return inflatedInstance;
    }

    // A copy that can be changed without changing the original: callers add values to the instances they inflate
    // before committing them. Its attribute lists are copied too; the instances in them are shells (dbId,
    // displayName and schemaClass only), which callers do not change.
    private static SimpleInstance copyOf(SimpleInstance instance) {
        SimpleInstance copy = instance.cloneInstance();
        if (copy.getAttributes() != null) {
            copy.getAttributes().replaceAll((attributeName, value) ->
                value instanceof List ? new ArrayList<>((List<?>) value) : value);
        }
        return copy;
    }

    /**
     * Returns the stored instances whose defining attributes match those of the instance, inflated. curator-tool-ws
     * only matches a reference attribute on the dbIds of the instances it refers to, so a referred instance that
     * has not been committed yet is left out of the match.
     *
     * @param instance - the instance whose defining attributes are searched on.
     * @return the matching instances, or an empty list if there are none or no defining attribute has a value.
     */
    public List<SimpleInstance> fetchIdenticalInstances(SimpleInstance instance) {
        // The matches come back as shells -- dbId, displayName and schemaClass only. They are inflated because
        // callers go on to add attribute values to them and commit them, and committing a shell would clear every
        // attribute of the stored instance.
        // The match converts the instance the same way a commit does, so a single value given for a multi-valued
        // attribute would be dropped from the match rather than matched on.
        prepareAttributeValues(instance);
        return timed("matchInstances", () -> controller.matchInstances(instance))
            .stream()
            .map(this::inflate)
            .filter(Objects::nonNull)
            .collect(Collectors.toList());
    }

    public List<SimpleInstance> getReferrers(SimpleInstance instance, String referrerAttributeName) throws Exception {
        return timed("getReferrers", () -> controller.getReferrers(instance.getDbId()))
            .stream()
            .filter(g -> referrerAttributeName.equals(g.getAttributeName()))
            .findFirst()
            .map(NamedReferrerList::getReferrers)
            .orElse(Collections.emptyList());
    }

    public Collection<NamedReferrerList> getReferrers(SimpleInstance instance) throws Exception {
        return timed("getReferrers", () -> controller.getReferrers(instance.getDbId()));
    }

    /**
     * Returns the names of every attribute the schema class has, including those it inherits.
     *
     * @param className - the schema class, e.g. "Compartment".
     * @return the attribute names, or an empty list if the class is not in the schema.
     */
    public List<String> getAttributeNames(String className) {
        return timed("getAttributes", () -> controller.getAttributes(className))
            .stream()
            .map(CurationAttribute::getName)
            .collect(Collectors.toList());
    }

    /**
     * Returns the names of the schema class's defining attributes (ALL_DEFINING or ANY_DEFINING), including those it
     * inherits -- the attributes curator-tool-ws matches on in fetchIdenticalInstances.
     *
     * @param className - the schema class, e.g. "Complex".
     * @return the defining attribute names, or an empty list if the class is not in the schema.
     */
    public List<String> getDefiningAttributeNames(String className) {
        return timed("getAttributes", () -> controller.getAttributes(className))
            .stream()
            .filter(attribute -> attribute.getDefiningType() == CurationAttribute.DefiningType.ALL_DEFINING ||
                attribute.getDefiningType() == CurationAttribute.DefiningType.ANY_DEFINING)
            .map(CurationAttribute::getName)
            .collect(Collectors.toList());
    }

    public long getPersonId() {
        return this.personId;
    }

    /**
     * Puts the instance's attribute values, and those of every instance it refers to, in the form curator-tool-ws's
     * DatabaseObjectInstanceConverter needs. The instances are changed in place:
     * <ul>
     *   <li>An attribute with a null value is removed. The converter looks up a setter by the value's class, so
     *   a null value fails with a NullPointerException. A commit replaces every attribute of the stored instance,
     *   so a removed attribute is stored the same way a null one would be.</li>
     *   <li>Null elements are removed from a multi-valued attribute's list, which the converter also cannot handle.</li>
     *   <li>A single value given for a multi-valued attribute is wrapped in a list. The graph model's setter for
     *   such an attribute takes a List, and the converter does not wrap the value itself: it logs "Cannot find a
     *   set method" and skips the attribute, so the value is not stored or matched on.</li>
     * </ul>
     *
     * @param instance - the instance about to be sent to curator-tool-ws.
     */
    private void prepareAttributeValues(SimpleInstance instance) {
        prepareAttributeValues(instance, Collections.newSetFromMap(new IdentityHashMap<>()));
    }

    private void prepareAttributeValues(SimpleInstance instance, Set<SimpleInstance> visitedInstances) {
        // Instances can refer to each other in cycles, so each is prepared only once.
        if (!visitedInstances.add(instance) || instance.getAttributes() == null) {
            return;
        }

        Set<String> multiValuedAttributeNames = getMultiValuedAttributeNames(instance.getSchemaClassName());
        Iterator<Map.Entry<String, Object>> attributeIterator = instance.getAttributes().entrySet().iterator();
        while (attributeIterator.hasNext()) {
            Map.Entry<String, Object> attribute = attributeIterator.next();
            Object value = attribute.getValue();
            if (value == null) {
                attributeIterator.remove();
                continue;
            }

            if (value instanceof Collection) {
                Collection<?> values = (Collection<?>) value;
                if (values.contains(null)) {
                    // A new list rather than removing from this one, which may be shared or unmodifiable.
                    List<Object> nonNullValues = values.stream().filter(Objects::nonNull).collect(Collectors.toList());
                    attribute.setValue(nonNullValues);
                    values = nonNullValues;
                }
                for (Object element : values) {
                    if (element instanceof SimpleInstance) {
                        prepareAttributeValues((SimpleInstance) element, visitedInstances);
                    }
                }
            } else {
                if (multiValuedAttributeNames.contains(attribute.getKey())) {
                    attribute.setValue(new ArrayList<>(Collections.singletonList(value)));
                }
                if (value instanceof SimpleInstance) {
                    prepareAttributeValues((SimpleInstance) value, visitedInstances);
                }
            }
        }
    }

    /**
     * Returns the names of the schema class's multi-valued attributes, including those it inherits. Whether an
     * attribute is multi-valued comes from the graph model, so it agrees with whether the setter the converter
     * uses takes a List.
     *
     * @param className - the schema class, e.g. "Reaction".
     * @return the multi-valued attribute names, or an empty set if the class is null or not in the schema.
     */
    private Set<String> getMultiValuedAttributeNames(String className) {
        if (className == null) {
            return Collections.emptySet();
        }
        return classNameToMultiValuedAttributeNames.computeIfAbsent(className, name ->
            timed("getAttributes", () -> controller.getAttributes(name))
                .stream()
                .filter(attribute -> attribute.getProperties() != null &&
                    "+".equals(attribute.getProperties().getCardinality()))
                .map(CurationAttribute::getName)
                .collect(Collectors.toSet())
        );
    }

    private void updateDisplayName(SimpleInstance referenceSequence) {
        // A search returns shell instances -- dbId, displayName and schemaClass only -- so the instance has to be
        // inflated before there are any attributes to build a display name out of. Committing a shell would also
        // clear every attribute of the stored instance, since a commit replaces all of them.
        SimpleInstance inflatedReferenceSequence = inflate(referenceSequence);
        if (inflatedReferenceSequence == null) {
            logger.warn("No instance found for db id " + referenceSequence.getDbId() +
                " -- skipping its display name update");
            return;
        }

        String currentDisplayName = inflatedReferenceSequence.getDisplayName();
        String newDisplayName = getReferenceSequenceDisplayName(inflatedReferenceSequence);

        if (!newDisplayName.equals(currentDisplayName)) {
            inflatedReferenceSequence.setDisplayName(newDisplayName);
            commit(inflatedReferenceSequence);
        }
    }

    public String getReferenceSequenceDisplayName(SimpleInstance referenceSequence) {
        String dbName = null;
        SimpleInstance refDB =
            (SimpleInstance) referenceSequence.getAttribute(ReactomeJavaConstants.referenceDatabase);
        if (refDB != null) {
            dbName = refDB.getDisplayName();
        }
        if (dbName == null) {
            dbName = "Unknown";
        }

        String identifier = getFirstAttributeValueAsString(referenceSequence, ReactomeJavaConstants.variantIdentifier);

        if (identifier == null) {
            identifier = getFirstAttributeValueAsString(referenceSequence, ReactomeJavaConstants.identifier);
        }
        if (identifier == null) {
            identifier = "Unknown";
        }

        // geneName and name are multi-valued attributes, so their values are Lists rather than the Strings this used
        // to cast them to; the first value is the one that belongs in the display name.
        String name = getFirstAttributeValueAsString(referenceSequence, ReactomeJavaConstants.geneName);
        if (name == null) {
            name = getFirstAttributeValueAsString(referenceSequence, ReactomeJavaConstants.name);
        }
        if (name == null) {
            name = "Unknown";
        }
        return dbName + ":" + identifier + " " + name;
    }

    public SimpleInstance fetchSpecies(String speciesName) {
        InstanceList speciesInstances = timed("searchInstances", () -> controller.searchInstances(
            ReactomeJavaConstants.Species,
            0,
            1,
            Optional.of(ReactomeJavaConstants.name),
            Optional.of("equal"),
            searchKeysRequest(speciesName)
        ));

        return !speciesInstances.isEmpty() ? speciesInstances.getInstances().get(0) : null;
    }

    /**
     * Returns the ReferenceDatabase instance with the name, as a shell instance.
     *
     * @param referenceDatabaseName - the name of the reference database, e.g. "UniProt".
     * @return the ReferenceDatabase shell instance, or null if there is none.
     */
    public SimpleInstance fetchReferenceDatabase(String referenceDatabaseName) {
        InstanceList referenceDatabaseInstances = timed("searchInstances", () -> controller.searchInstances(
            ReactomeJavaConstants.ReferenceDatabase,
            0,
            1,
            Optional.of(ReactomeJavaConstants.name),
            Optional.of("equal"),
            searchKeysRequest(referenceDatabaseName)
        ));

        return !referenceDatabaseInstances.isEmpty() ? referenceDatabaseInstances.getInstances().get(0) : null;
    }

    /**
     * Returns the StableIdentifier instance with the identifier, as a shell instance.
     *
     * @param identifier - the stable identifier, e.g. "R-MMU-123456".
     * @return the StableIdentifier shell instance, or null if there is none.
     */
    public SimpleInstance fetchStableIdentifier(String identifier) {
        InstanceList stableIdentifierInstances = timed("searchInstances", () -> controller.searchInstances(
            ReactomeJavaConstants.StableIdentifier,
            0,
            1,
            Optional.of(ReactomeJavaConstants.identifier),
            Optional.of("equal"),
            searchKeysRequest(identifier)
        ));

        return !stableIdentifierInstances.isEmpty() ? stableIdentifierInstances.getInstances().get(0) : null;
    }

    /**
     * Returns every ReactionlikeEvent (of any subclass) whose species is Homo sapiens, as shell instances -- dbId,
     * displayName and schemaClass only -- which have to be inflated before any other attribute can be read.
     *
     * @return the human ReactionlikeEvent shell instances, or an empty list if there are none.
     */
    public List<SimpleInstance> fetchHumanReactionlikeEvents() {
        return fetchInstancesBySpecies("ReactionLikeEvent", "Homo sapiens");
    }

    /**
     * Returns every instance of the class (or one of its subclasses) whose species has the given name, as shell
     * instances. curator-tool-ws matches an instance-valued attribute on the referred instance's display name, so
     * the species is given by name rather than by dbId.
     *
     * @param className - the schema class to search.
     * @param speciesName - the display name of the species, e.g. "Homo sapiens".
     * @return the matching shell instances, or an empty list if there are none.
     */
    public List<SimpleInstance> fetchInstancesBySpecies(String className, String speciesName) {
        // Results are ordered by display name, which is not unique, so paging could skip or repeat instances sharing
        // a name across a page boundary. The total count is read first and everything is then fetched in one request.
        Integer totalCount = searchBySpecies(className, speciesName, 0, 1).getTotalCount();
        if (totalCount == null || totalCount == 0) {
            return Collections.emptyList();
        }

        InstanceList instanceList = searchBySpecies(className, speciesName, 0, totalCount);
        return instanceList.isEmpty() ? Collections.emptyList() : instanceList.getInstances();
    }

    private InstanceList searchBySpecies(String className, String speciesName, int skip, int limit) {
        return timed("searchInstances", () -> controller.searchInstances(
            className,
            skip,
            limit,
            Optional.of(ReactomeJavaConstants.species),
            Optional.of("equal"),
            searchKeysRequest(speciesName)
        ));
    }

    /**
     * Returns a stand-in for the HTTP request that CurationController.searchInstances reads its search keys from.
     * The controller takes them off the request as repeated "searchKeys" parameters, rather than as an argument,
     * so that a key containing commas (such as a regex quantifier) is not split -- and that parameter is all it
     * reads. Any other use of the request fails loudly, rather than handing back a null the controller would
     * then trip over.
     *
     * @param searchKeys - the search keys, one per attribute searched on.
     * @return a request carrying the search keys as its "searchKeys" parameter values.
     */
    private HttpServletRequest searchKeysRequest(String... searchKeys) {
        return (HttpServletRequest) Proxy.newProxyInstance(
            HttpServletRequest.class.getClassLoader(),
            new Class<?>[] {HttpServletRequest.class},
            (proxy, method, args) -> {
                if (method.getName().equals("getParameterValues") && "searchKeys".equals(args[0])) {
                    return searchKeys.clone();
                }
                if (method.getName().equals("getParameterValues")) {
                    return null;
                }
                if (method.getName().equals("toString")) {
                    return "searchKeys request " + Arrays.toString(searchKeys);
                }
                throw new UnsupportedOperationException(
                    "CuratorToolAPI's search request does not support " + method.getName());
            }
        );
    }

    /**
     * Returns the instance's first value for the attribute as a string, whether the attribute is multi-valued (its
     * value a list) or single-valued.
     *
     * @param instance - the instance to read.
     * @param attributeName - the attribute, e.g. "name".
     * @return the first value as a string, or null if the attribute has no value.
     */
    public static String getFirstAttributeValueAsString(SimpleInstance instance, String attributeName) {
        List<Object> values = getAttributeValues(instance, attributeName);
        return values.isEmpty() ? null : values.get(0).toString();
    }

    private static List<Object> getAttributeValues(SimpleInstance instance, String attributeName) {
        Object value = instance.getAttribute(attributeName);
        if (value == null) {
            return Collections.emptyList();
        }
        return value instanceof List ? new ArrayList<>((List<?>) value) : Collections.singletonList(value);
    }

    /**
     * Logs, for each kind of curator-tool-ws call made so far, how many times it was called and the time spent in it,
     * most time first, with the time elapsed since this API was created for comparison.
     *
     * @param header - a line logged above the figures, e.g. how far through the run they were taken.
     */
    public void logCallStats(String header) {
        long elapsedMillis = (System.nanoTime() - startNanos) / 1_000_000;
        StringBuilder report = new StringBuilder(header)
            .append(" -- curator-tool-ws calls after ").append(elapsedMillis).append(" ms:");
        callStats.entrySet().stream()
            .sorted((a, b) -> Long.compare(b.getValue().totalNanos, a.getValue().totalNanos))
            .forEach(entry -> {
                CallStats stats = entry.getValue();
                long totalMillis = stats.totalNanos / 1_000_000;
                report.append(String.format("%n  %-20s calls=%-9d total=%-9d ms (%5.1f%%) avg=%.2f ms",
                    entry.getKey(), stats.count, totalMillis,
                    elapsedMillis > 0 ? 100.0 * totalMillis / elapsedMillis : 0.0,
                    stats.count > 0 ? stats.totalNanos / 1_000_000.0 / stats.count : 0.0));
            });
        logger.info(report.toString());
    }

    // Runs the call, adding its duration to the figures kept under the label. A call made inside another timed call
    // (e.g. a nested commit) is counted under both labels.
    private <T, E extends Exception> T timed(String label, ThrowingSupplier<T, E> call) throws E {
        long start = System.nanoTime();
        try {
            return call.get();
        } finally {
            CallStats stats = callStats.computeIfAbsent(label, key -> new CallStats());
            stats.count++;
            stats.totalNanos += System.nanoTime() - start;
        }
    }

    @FunctionalInterface
    private interface ThrowingSupplier<T, E extends Exception> {
        T get() throws E;
    }

    private static final class CallStats {
        private long count;
        private long totalNanos;
    }
}
