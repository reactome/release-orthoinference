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

        SimpleInstance committedInstance = controller.commit(simpleInstance);
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
        if (!controller.existsByDbId(dbId)) {
            return null;
        }
        return controller.findByDdIdInInstance(dbId);
    }

    public SimpleInstance findByDisplayName(String className, String displayName) {
        // NB: the controller takes the display name first and a comma-separated list of class names second.
        return controller.findByDisplayName(displayName, className);
    }

    public void deleteInstance(SimpleInstance instance) {
        if (instance.getDefaultPersonId() == null) {
            instance.setDefaultPersonId(getPersonId());
        }

        controller.delete(instance);
        if (instance.getDbId() != null) {
            deletedDbIds.add(instance.getDbId());
            committedDbIds.remove(instance.getDbId());
        }
    }

    public void deleteByDbId(long noReferrerDbId) {
        SimpleInstance instance = controller.findByDdIdInInstance(noReferrerDbId);
        // Nothing stored under that db id -- there is nothing to delete, and passing the null on would fail inside
        // curator-tool-ws instead of saying which db id was asked for.
        if (instance == null) {
            logger.warn("No instance found for db id " + noReferrerDbId + " -- nothing to delete");
            return;
        }
        deleteInstance(instance);
    }

    public SimpleInstance inflate(SimpleInstance shellInstance) {
        return controller.findByDdIdInInstance(shellInstance.getDbId());
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
        return controller.matchInstances(instance)
            .stream()
            .map(this::inflate)
            .filter(Objects::nonNull)
            .collect(Collectors.toList());
    }

    public List<SimpleInstance> getReferrers(SimpleInstance instance, String referrerAttributeName) throws Exception {
        return controller.getReferrers(instance.getDbId())
            .stream()
            .filter(g -> referrerAttributeName.equals(g.getAttributeName()))
            .findFirst()
            .map(NamedReferrerList::getReferrers)
            .orElse(Collections.emptyList());
    }

    public Collection<NamedReferrerList> getReferrers(SimpleInstance instance) throws Exception {
        return controller.getReferrers(instance.getDbId());
    }

    /**
     * Returns the names of every attribute the schema class has, including those it inherits.
     *
     * @param className - the schema class, e.g. "Compartment".
     * @return the attribute names, or an empty list if the class is not in the schema.
     */
    public List<String> getAttributeNames(String className) {
        return controller.getAttributes(className)
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
        return controller.getAttributes(className)
            .stream()
            .filter(attribute -> attribute.getDefiningType() == CurationAttribute.DefiningType.ALL_DEFINING ||
                attribute.getDefiningType() == CurationAttribute.DefiningType.ANY_DEFINING)
            .map(CurationAttribute::getName)
            .collect(Collectors.toList());
    }

    public long getPersonId() {
        return this.personId;
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
        InstanceList speciesInstances = controller.searchInstances(
            ReactomeJavaConstants.Species,
            0,
            1,
            Optional.of(ReactomeJavaConstants.name),
            Optional.of("equal"),
            searchKeysRequest(speciesName)
        );

        return !speciesInstances.isEmpty() ? speciesInstances.getInstances().get(0) : null;
    }

    /**
     * Returns the ReferenceDatabase instance with the name, as a shell instance.
     *
     * @param referenceDatabaseName - the name of the reference database, e.g. "UniProt".
     * @return the ReferenceDatabase shell instance, or null if there is none.
     */
    public SimpleInstance fetchReferenceDatabase(String referenceDatabaseName) {
        InstanceList referenceDatabaseInstances = controller.searchInstances(
            ReactomeJavaConstants.ReferenceDatabase,
            0,
            1,
            Optional.of(ReactomeJavaConstants.name),
            Optional.of("equal"),
            searchKeysRequest(referenceDatabaseName)
        );

        return !referenceDatabaseInstances.isEmpty() ? referenceDatabaseInstances.getInstances().get(0) : null;
    }

    /**
     * Returns the StableIdentifier instance with the identifier, as a shell instance.
     *
     * @param identifier - the stable identifier, e.g. "R-MMU-123456".
     * @return the StableIdentifier shell instance, or null if there is none.
     */
    public SimpleInstance fetchStableIdentifier(String identifier) {
        InstanceList stableIdentifierInstances = controller.searchInstances(
            ReactomeJavaConstants.StableIdentifier,
            0,
            1,
            Optional.of(ReactomeJavaConstants.identifier),
            Optional.of("equal"),
            searchKeysRequest(identifier)
        );

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
        return controller.searchInstances(
            className,
            skip,
            limit,
            Optional.of(ReactomeJavaConstants.species),
            Optional.of("equal"),
            searchKeysRequest(speciesName)
        );
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

    private static String getFirstAttributeValueAsString(SimpleInstance instance, String attributeName) {
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
}
