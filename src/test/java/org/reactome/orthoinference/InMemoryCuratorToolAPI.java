package org.reactome.orthoinference;

import com.fasterxml.jackson.core.type.TypeReference;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.gk.model.ReactomeJavaConstants;
import org.reactome.curation.model.SimpleInstance;
import org.reactome.server.graph.domain.model.DatabaseObject;

import java.io.IOException;
import java.io.InputStream;
import java.io.UncheckedIOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * A CuratorToolAPI that keeps its "database" in memory, so that the inference classes can be run end to end without
 * curator-tool-ws or a graph database.
 *
 * <p>It behaves as curator-tool-ws does in the ways the inference relies on:
 * <ul>
 *   <li>Stored instances refer to other instances as shells (dbId, displayName and schemaClass only), and
 *   {@link #inflate} hands out copies, so changing an inflated instance does not change what is stored.</li>
 *   <li>A commit gives a new instance a positive dbId, copies it back onto the caller's instance, and also stores the
 *   new instances it refers to (those with no dbId or a negative placeholder dbId, e.g. a new StableIdentifier).</li>
 *   <li>{@link #fetchIdenticalInstances} matches on the schema's defining attributes, which are read from the
 *   curation_schema_attributes.json that curator-tool-ws itself uses.</li>
 * </ul>
 *
 * <p>Calls the inference does not make are left to the superclass, whose controller is null, so a call the fake does
 * not support fails loudly rather than returning something misleading.
 */
public class InMemoryCuratorToolAPI extends CuratorToolAPI {

	private static final long FIRST_ASSIGNED_DB_ID = 900_000_000L;
	private static final Map<String, List<SchemaAttribute>> SCHEMA = loadSchema();

	private final Map<Long, SimpleInstance> storedInstances = new LinkedHashMap<>();
	private final List<SimpleInstance> commits = new ArrayList<>();
	private long nextDbId = FIRST_ASSIGNED_DB_ID;
	private long nextPlaceholderDbId = -1;

	public InMemoryCuratorToolAPI(ConfigProperties configProperties) {
		super(null, configProperties);
	}

	// ---- Test helpers -------------------------------------------------------------------------------------------

	/**
	 * Stores a fixture instance with the dbId it already has, as if it were in the database before the run.
	 */
	public SimpleInstance add(SimpleInstance instance) {
		if (instance.getDbId() == null || instance.getDbId() <= 0) {
			throw new IllegalArgumentException("Fixture instances need a positive dbId: " + instance);
		}
		storedInstances.put(instance.getDbId(), toStoredCopy(instance));
		return instance;
	}

	/**
	 * @return a copy of the stored instance, with the instances it refers to as shells, or null if none is stored.
	 */
	public SimpleInstance stored(long dbId) {
		SimpleInstance storedInstance = storedInstances.get(dbId);
		return storedInstance != null ? deepCopy(storedInstance) : null;
	}

	/**
	 * @return copies of the instances committed (not added as fixtures) of exactly the schema class, in commit order.
	 */
	public List<SimpleInstance> committedInstancesOfClass(String schemaClassName) {
		return storedInstances.values().stream()
			.filter(instance -> instance.getDbId() >= FIRST_ASSIGNED_DB_ID)
			.filter(instance -> schemaClassName.equals(instance.getSchemaClassName()))
			.map(InMemoryCuratorToolAPI::deepCopy)
			.collect(Collectors.toList());
	}

	/**
	 * @return the instances passed to {@link #commit}, in order (the caller's objects, not copies).
	 */
	public List<SimpleInstance> getCommits() {
		return Collections.unmodifiableList(commits);
	}

	// ---- CuratorToolAPI ----------------------------------------------------------------------------------------

	@Override
	public SimpleInstance commit(SimpleInstance instance) {
		if (instance.getDefaultPersonId() == null) {
			instance.setDefaultPersonId(getPersonId());
		}
		store(instance, Collections.newSetFromMap(new IdentityHashMap<>()));
		commits.add(instance);
		return instance;
	}

	@Override
	public void assignPlaceholderDbId(SimpleInstance instance) {
		if (instance.getDbId() == null) {
			instance.setDbId(nextPlaceholderDbId--);
		}
	}

	@Override
	public SimpleInstance refresh(SimpleInstance instance) {
		if (instance.getDbId() == null || !storedInstances.containsKey(instance.getDbId())) {
			return instance;
		}
		return inflate(instance);
	}

	@Override
	public SimpleInstance findByDbId(long dbId) {
		return stored(dbId);
	}

	@Override
	public SimpleInstance inflate(SimpleInstance shellInstance) {
		return shellInstance.getDbId() != null ? stored(shellInstance.getDbId()) : null;
	}

	@Override
	public List<SimpleInstance> fetchIdenticalInstances(SimpleInstance instance) {
		String className = instance.getSchemaClassName();
		List<String> allDefining = getAllDefiningAttributeNames(className);
		List<String> anyDefining = getDefiningAttributeNames(className).stream()
			.filter(attributeName -> !allDefining.contains(attributeName))
			.collect(Collectors.toList());

		boolean hasDefiningValue = getDefiningAttributeNames(className).stream()
			.anyMatch(attributeName -> !matchKeys(instance, attributeName).isEmpty());
		if (!hasDefiningValue) {
			return new ArrayList<>();
		}

		Class<? extends DatabaseObject> modelClass = instance.getGraphModelClass();
		return storedInstances.values().stream()
			.filter(candidate -> modelClass.isAssignableFrom(candidate.getGraphModelClass()))
			.filter(candidate -> allDefining.stream().allMatch(attributeName ->
				matchKeys(instance, attributeName).equals(matchKeys(candidate, attributeName))))
			.filter(candidate -> anyDefining.stream().allMatch(attributeName -> {
				List<String> keys = matchKeys(instance, attributeName);
				return keys.isEmpty() || !Collections.disjoint(keys, matchKeys(candidate, attributeName));
			}))
			.map(InMemoryCuratorToolAPI::deepCopy)
			.collect(Collectors.toList());
	}

	@Override
	public List<SimpleInstance> getReferrers(SimpleInstance instance, String referrerAttributeName) {
		return storedInstances.values().stream()
			.filter(referrer -> valuesOf(referrer, referrerAttributeName).stream()
				.anyMatch(value -> value instanceof SimpleInstance &&
					instance.getDbId().equals(((SimpleInstance) value).getDbId())))
			.map(InMemoryCuratorToolAPI::shellOf)
			.collect(Collectors.toList());
	}

	@Override
	public List<String> getAttributeNames(String className) {
		return attributesOf(className).stream().map(SchemaAttribute::getName).collect(Collectors.toList());
	}

	@Override
	public List<String> getDefiningAttributeNames(String className) {
		return attributesOf(className).stream()
			.filter(attribute -> attribute.isDefining("ALL_DEFINING") || attribute.isDefining("ANY_DEFINING"))
			.map(SchemaAttribute::getName)
			.collect(Collectors.toList());
	}

	@Override
	public List<String> getAllDefiningAttributeNames(String className) {
		return attributesOf(className).stream()
			.filter(attribute -> attribute.isDefining("ALL_DEFINING"))
			.map(SchemaAttribute::getName)
			.collect(Collectors.toList());
	}

	@Override
	public SimpleInstance fetchSpecies(String speciesName) {
		return findShell(ReactomeJavaConstants.Species, ReactomeJavaConstants.name, speciesName);
	}

	@Override
	public SimpleInstance fetchReferenceDatabase(String referenceDatabaseName) {
		return findShell(ReactomeJavaConstants.ReferenceDatabase, ReactomeJavaConstants.name, referenceDatabaseName);
	}

	@Override
	public SimpleInstance fetchStableIdentifier(String identifier) {
		return findShell(ReactomeJavaConstants.StableIdentifier, ReactomeJavaConstants.identifier, identifier);
	}

	// As curator-tool-ws does, the class matches its subclasses, the species is matched on its display name, and the
	// matches are returned as shells.
	@Override
	public List<SimpleInstance> fetchInstancesBySpecies(String className, String speciesName) {
		SimpleInstance classProbe = new SimpleInstance();
		classProbe.setSchemaClassName(className);
		Class<?> modelClass = classProbe.getGraphModelClass();
		return storedInstances.values().stream()
			.filter(instance -> modelClass.isAssignableFrom(instance.getGraphModelClass()))
			.filter(instance -> valuesOf(instance, ReactomeJavaConstants.species).stream()
				.anyMatch(species -> species instanceof SimpleInstance &&
					speciesName.equals(((SimpleInstance) species).getDisplayName())))
			.map(InMemoryCuratorToolAPI::shellOf)
			.collect(Collectors.toList());
	}

	@Override
	public void logCallStats(String header) {
		// Nothing is timed.
	}

	// ---- Storage ------------------------------------------------------------------------------------------------

	// New referred instances are stored first, so the stored copy of this instance refers to their real dbIds.
	private void store(SimpleInstance instance, Set<SimpleInstance> visited) {
		if (!visited.add(instance)) {
			return;
		}
		if (instance.getAttributes() != null) {
			for (Object value : instance.getAttributes().values()) {
				for (Object element : asList(value)) {
					if (element instanceof SimpleInstance && isNew((SimpleInstance) element)) {
						store((SimpleInstance) element, visited);
					}
				}
			}
		}
		if (isNew(instance)) {
			instance.setDbId(nextDbId++);
		}
		storedInstances.put(instance.getDbId(), toStoredCopy(instance));
	}

	private static boolean isNew(SimpleInstance instance) {
		return instance.getDbId() == null || instance.getDbId() < 0;
	}

	private SimpleInstance findShell(String className, String attributeName, String value) {
		return storedInstances.values().stream()
			.filter(instance -> className.equals(instance.getSchemaClassName()))
			.filter(instance -> valuesOf(instance, attributeName).stream().anyMatch(v -> value.equals(String.valueOf(v))))
			.findFirst()
			.map(InMemoryCuratorToolAPI::shellOf)
			.orElse(null);
	}

	// An attribute's values as compared by fetchIdenticalInstances: instances by dbId, other values as strings,
	// sorted so that the order of a multi-valued attribute does not matter. Referred instances that are not stored
	// yet are left out, as curator-tool-ws leaves them out of the match.
	private static List<String> matchKeys(SimpleInstance instance, String attributeName) {
		return valuesOf(instance, attributeName).stream()
			.filter(Objects::nonNull)
			.filter(value -> !(value instanceof SimpleInstance) || !isNew((SimpleInstance) value))
			.map(value -> value instanceof SimpleInstance ?
				"#" + ((SimpleInstance) value).getDbId() : String.valueOf(value))
			.sorted()
			.collect(Collectors.toList());
	}

	private static List<Object> valuesOf(SimpleInstance instance, String attributeName) {
		return asList(instance.getAttribute(attributeName));
	}

	private static List<Object> asList(Object value) {
		if (value == null) {
			return Collections.emptyList();
		}
		return value instanceof List ? new ArrayList<>((List<?>) value) : Collections.singletonList(value);
	}

	// The stored copy, and every copy handed out, refers to other instances as shells, and has its own lists, so
	// neither the caller's instance nor the stored one changes when the other does.
	private static SimpleInstance toStoredCopy(SimpleInstance instance) {
		return deepCopy(instance);
	}

	private static SimpleInstance deepCopy(SimpleInstance instance) {
		SimpleInstance copy = shellOf(instance);
		copy.setDefaultPersonId(instance.getDefaultPersonId());
		if (instance.getAttributes() != null) {
			Map<String, Object> attributes = new HashMap<>();
			instance.getAttributes().forEach((name, value) -> attributes.put(name, toShells(value)));
			copy.setAttributes(attributes);
		}
		return copy;
	}

	private static Object toShells(Object value) {
		if (value instanceof List) {
			return ((List<?>) value).stream().map(InMemoryCuratorToolAPI::toShells)
				.collect(Collectors.toCollection(ArrayList::new));
		}
		return value instanceof SimpleInstance ? shellOf((SimpleInstance) value) : value;
	}

	private static SimpleInstance shellOf(SimpleInstance instance) {
		SimpleInstance shell = new SimpleInstance();
		shell.setDbId(instance.getDbId());
		shell.setDisplayName(instance.getDisplayName());
		shell.setSchemaClassName(instance.getSchemaClassName());
		return shell;
	}

	// ---- Schema -------------------------------------------------------------------------------------------------

	private static List<SchemaAttribute> attributesOf(String className) {
		return SCHEMA.getOrDefault(className, Collections.emptyList());
	}

	private static Map<String, List<SchemaAttribute>> loadSchema() {
		try (InputStream schemaStream =
				 InMemoryCuratorToolAPI.class.getClassLoader().getResourceAsStream("curation_schema_attributes.json")) {
			if (schemaStream == null) {
				throw new IllegalStateException("curation_schema_attributes.json (from curator-tool-ws) not found");
			}
			return new ObjectMapper().readValue(schemaStream, new TypeReference<Map<String, List<SchemaAttribute>>>() {});
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
	}

	@com.fasterxml.jackson.annotation.JsonIgnoreProperties(ignoreUnknown = true)
	static class SchemaAttribute {
		private String name;
		private String definingType;

		public String getName() {
			return name;
		}

		public void setName(String name) {
			this.name = name;
		}

		public void setDefiningType(String definingType) {
			this.definingType = definingType;
		}

		boolean isDefining(String type) {
			return type.equals(definingType);
		}
	}
}
