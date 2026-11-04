package org.reactome.orthoinference;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.gk.model.ReactomeJavaConstants;
import org.json.simple.parser.ParseException;
import org.reactome.curation.model.SimpleInstance;
import org.reactome.server.graph.domain.model.DatabaseObject;
import org.reactome.server.graph.domain.model.PhysicalEntity;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.context.annotation.Bean;
import org.springframework.stereotype.Component;
import org.springframework.util.ReflectionUtils;
import org.springframework.util.StringUtils;

// GenerateInstance is meant to act as a catch-all for functions that are instance-oriented, such as creating, mocking,
// or identical-checking.
@Component
public class InstanceUtilities {

	private static final Logger logger = LogManager.getLogger();

	private final long DISEASE_PATHWAY_DB_ID = 1643685L;
	// Ends each value in a cache key, so that e.g. values "12","3" and "1","23" do not produce the same key.
	private static final String CACHE_KEY_SEPARATOR = "|";

	private ConfigProperties configProperties;
	private String targetSpeciesCode;
	private StableIdentifierGenerator stableIdentifierGenerator;
	private CuratorToolAPI curatorToolAPI;

	private SimpleInstance evidenceType;
	private SimpleInstance summation;

	private SimpleInstance speciesInst;
	private Map<String,SimpleInstance> mockedIdenticals = new HashMap<>();

	public InstanceUtilities(
		ConfigProperties configProperties,
		@Qualifier("targetSpeciesCode") String targetSpeciesCode,
		StableIdentifierGenerator stableIdentifierGenerator,
		CuratorToolAPI curatorToolAPI
	) {
		this.configProperties = configProperties;
		this.targetSpeciesCode = targetSpeciesCode;
		this.stableIdentifierGenerator = stableIdentifierGenerator;
		this.curatorToolAPI = curatorToolAPI;
	}

//	@Bean(name = "instanceEditInst")
//	public GKInstance getInstanceEdit() throws Exception {
//		if (instanceEdit == null) {
//			instanceEdit = InstanceEditUtils.createInstanceEdit(
//				getCurrentDBA(), getPersonId(), "org.reactome.orthoinference");
//		}
//		return instanceEdit;
//	}

	// Find the instance specific to this species
	@Bean(name = "speciesInst")
	public SimpleInstance getSpeciesInstance() throws Exception {
		if (speciesInst == null) {
			speciesInst = new SimpleInstance();
			speciesInst.setSchemaClassName(ReactomeJavaConstants.Species);
			speciesInst.setDefaultPersonId(getCuratorToolAPI().getPersonId());
			speciesInst.setAttribute(ReactomeJavaConstants.name, getSpeciesName());
			speciesInst.setDisplayName(getSpeciesName());
			speciesInst = checkForIdenticalInstances(speciesInst, null);
		}
		return speciesInst;
	}

	// Creates new instance that will be inferred based on the incoming instances class
	public SimpleInstance createNewInferredGKInstance(SimpleInstance instanceToBeInferred) throws Exception {
		String reactionClass = instanceToBeInferred.getSchemaClassName();
		if (reactionClass.matches(ReactomeJavaConstants.ReferenceIsoform)) {
			reactionClass = ReactomeJavaConstants.ReferenceGeneProduct;
		}
		SimpleInstance inferredInst = new SimpleInstance();
		inferredInst.setSchemaClassName(reactionClass);

		if (instanceToBeInferred.getAttribute(ReactomeJavaConstants.compartment) != null) {
			inferredInst.setAttribute(ReactomeJavaConstants.compartment, getCompartments(instanceToBeInferred));
		}
		if (instanceToBeInferred.getAttribute(ReactomeJavaConstants.species) != null) {
			inferredInst.setAttribute(ReactomeJavaConstants.species, speciesInst);
		}
		return inferredInst;
	}
	
	// Some 'Compartment' instances were actually 'GO_CellularComponent' instances. This meant that the instances that
	// were pulled from the original instance's Compartment attribute could not be added to the new instance, due to
	// them being a GO_CellularComponent. This function is the workaround, producing a Compartment instance that
	// contains all the same attribute values.
	public SimpleInstance createCompartmentInstance(SimpleInstance compartmentInst) throws Exception {
		logger.warn(compartmentInst + " is a " + compartmentInst.getSchemaClassName() + " instead of a Compartment" +
			" -- creating new Compartment instance");
		SimpleInstance newCompartmentInst = new SimpleInstance();
		newCompartmentInst.setSchemaClassName(ReactomeJavaConstants.Compartment);
		newCompartmentInst.setDisplayName(compartmentInst.getDisplayName());
		// Only attributes that Compartment also has are copied. The dbId is left unset so the new instance is treated
		// as new when it is committed, and the display name, like the dbId, is held outside the attributes map.
		for (String attributeName : getCuratorToolAPI().getAttributeNames(ReactomeJavaConstants.Compartment)) {
			if (attributeName.equals("dbId") || attributeName.equals("displayName")) {
				continue;
			}
			Object attributeValue = compartmentInst.getAttribute(attributeName);
			if (attributeValue != null) {
				// A multi-valued attribute's list is copied so the two instances do not share it.
				newCompartmentInst.setAttribute(attributeName,
					attributeValue instanceof List ? new ArrayList<>((List<?>) attributeValue) : attributeValue);
			}
		}
		return checkForIdenticalInstances(newCompartmentInst, null);
	}

	// Equivalent to create_ghost from Perl; Returns a mock homologue that is needed in cases where an inference is
	// rejected, but the component isn't essential for the inference to be completed.
	public SimpleInstance createMockGKInstance(SimpleInstance instanceToBeMocked) throws Exception {
		SimpleInstance mockedInst = new SimpleInstance();
		mockedInst.setSchemaClassName(ReactomeJavaConstants.GenomeEncodedEntity);
		mockedInst.setDefaultPersonId(getCuratorToolAPI().getPersonId());
		mockedInst.setDisplayName("Ghost homologue of " + instanceToBeMocked.getDisplayName());

		// A shell (dbId, displayName and schemaClass only) has no name or compartment to copy, so one is inflated
		// first. Only the copied values come from the inflated instance; the instance itself is still the one used.
		SimpleInstance sourceInst = instanceToBeMocked;
		if (sourceInst.getAttributes() == null && sourceInst.getDbId() != null && sourceInst.getDbId() > 0) {
			SimpleInstance inflatedInst = getCuratorToolAPI().inflate(sourceInst);
			if (inflatedInst != null) {
				sourceInst = inflatedInst;
			}
		}

		// name is multi-valued, so its value is a list; the first name is the one the ghost is named after.
		String mockedInstName = CuratorToolAPI.getFirstAttributeValueAsString(sourceInst, ReactomeJavaConstants.name);
		mockedInst.setAttribute(ReactomeJavaConstants.name, "Ghost homologue of " + mockedInstName);
		mockedInst.setAttribute(ReactomeJavaConstants.inferredFrom, instanceToBeMocked);
		mockedInst.setAttribute(ReactomeJavaConstants.species, speciesInst);
		mockedInst.setAttribute(ReactomeJavaConstants.compartment,
			sourceInst.getAttribute(ReactomeJavaConstants.compartment));

		// Caching based on an instance's defining attributes. This reduces the number of 'checkForIdenticalInstance'
		// calls, which is slow.
		String cacheKey = getCacheKey(mockedInst);
		if (mockedIdenticals.get(cacheKey) != null) {
			mockedInst = mockedIdenticals.get(cacheKey);
		} else {
			mockedInst = checkForIdenticalInstances(mockedInst, instanceToBeMocked);
			mockedIdenticals.put(cacheKey, mockedInst);
		}
		addInferredToValue(instanceToBeMocked, mockedInst);

		return mockedInst;
	}

	// Checks that equivalent instances don't already exist in the DB, substituting if they do
	public SimpleInstance checkForIdenticalInstances(SimpleInstance inferredInst, SimpleInstance originalInst)
		throws Exception {
		// An instance with a positive dbId is already stored -- e.g. the single inferred member that
		// createInfEntitySet substitutes for its EntitySet -- so it is its own identical instance. curator-tool-ws
		// does not always match an instance to itself, and committing it again here would send a copy whose
		// "modified" InstanceEdit is stale (later commits went through refreshed copies), which curator-tool-ws
		// rejects with an InstanceChangedException.
		if (inferredInst.getDbId() != null && inferredInst.getDbId() > 0) {
			return inferredInst;
		}
		List<SimpleInstance> identicalInstances = getCuratorToolAPI().fetchIdenticalInstances(inferredInst);
		if (!identicalInstances.isEmpty()) {
			// TODO: In future, could iterate through array of returned values and pull the 'most identical'.
			// For now, this mimics Perl.
			return identicalInstances.get(0);
		}
		if (isPhysicalEntity(inferredInst)) {
			SimpleInstance orthoStableIdentifierInst = stableIdentifierGenerator
				.generateOrthologousStableId(inferredInst, originalInst);
			// Null when inferredInst already has a stable identifier, which is then left as it is.
			if (orthoStableIdentifierInst != null) {
				inferredInst.setAttribute(ReactomeJavaConstants.stableIdentifier, orthoStableIdentifierInst);
			}
		}
		// Commit copies the dbId the database assigns back onto inferredInst.
		getCuratorToolAPI().commit(inferredInst);
		return inferredInst;
	}

	private boolean isPhysicalEntity(SimpleInstance instance) {
		return PhysicalEntity.class.isAssignableFrom(instance.getGraphModelClass());
	}

	// Adds the inferred instance to the source instance's inferredTo attribute and commits the source instance.
	public void addInferredToValue(SimpleInstance sourceInst, SimpleInstance inferredInst) {
		addAttributeValueAndCommit(sourceInst, inferredInst, ReactomeJavaConstants.inferredTo);
	}

	// Adds the source instance to the inferred instance's inferredFrom attribute and commits the inferred instance.
	public void addInferredFromValue(SimpleInstance inferredInst, SimpleInstance sourceInst) {
		addAttributeValueAndCommit(inferredInst, sourceInst, ReactomeJavaConstants.inferredFrom);
	}

	// Adds the orthologous instance to the event's orthologousEvent attribute and commits the event.
	public void addOrthologousEventValue(SimpleInstance eventInst, SimpleInstance orthologousEventInst) {
		addAttributeValueAndCommit(eventInst, orthologousEventInst, ReactomeJavaConstants.orthologousEvent);
	}

	// An attribute is valid if the graph model class for the instance's schema class has a getter for it
	public boolean isValidAttribute(SimpleInstance instance, String attributeName) {
		Class<?> graphModelClass = instance.getGraphModelClass();
		return graphModelClass != null &&
			ReflectionUtils.findMethod(graphModelClass, "get" + StringUtils.capitalize(attributeName)) != null;
	}

	// True if the entity's schema class is schemaClassType or any subclass of it. The graph model classes mirror the
	// graph database's schema hierarchy, so class assignability is equivalent to 'isa' on the schema class.
	public boolean isSchemaClassType(SimpleInstance entity, String schemaClassType) {
		Class<? extends DatabaseObject> entityClass = entity.getGraphModelClass();
		if (entityClass == null) {
			return false;
		}

		try {
			Class<?> schemaClass = Class.forName(DatabaseObject.class.getPackageName() + "." + schemaClassType);
			return schemaClass.isAssignableFrom(entityClass);
		} catch (ClassNotFoundException e) {
			throw new IllegalArgumentException("Unknown schema class: " + schemaClassType, e);
		}
	}

	// A commit replaces every attribute of the stored instance, so the value is added to an up-to-date copy: if this
	// run has already committed the instance, the copy the caller holds is stale and committing it would be rejected
	// or write its old values back over the stored ones. The caller's copy is left as it is.
	private void addAttributeValueAndCommit(SimpleInstance instance, SimpleInstance value, String attribute) {
		SimpleInstance currentInstance = getCuratorToolAPI().refresh(instance);
		// Null when this run has deleted the instance, and committing it would re-create it.
		if (currentInstance == null) {
			logger.warn(instance + " has been deleted -- not adding " + value + " to its " + attribute);
			return;
		}
		// Committing an unchanged instance would cost a full write and change nothing.
		if (hasAttributeValue(currentInstance, value, attribute)) {
			return;
		}
		addAttributeValueIfNecessary(currentInstance, value, attribute);
		getCuratorToolAPI().commit(currentInstance);
	}

	// True if the attribute already holds an instance with the value's dbId. A SimpleInstance holds an unset attribute
	// as null and a single value as itself rather than in a List; getAttributeValues handles both.
	private boolean hasAttributeValue(SimpleInstance instance, SimpleInstance value, String attribute) {
		Long valueDbId = value.getDbId();
		for (SimpleInstance attributeValue : getAttributeValues(instance, attribute)) {
			if (attributeValue.getDbId() != null && attributeValue.getDbId().equals(valueDbId)) {
				return true;
			}
		}
		return false;
	}

	// Returns a shell of the instance -- dbId, displayName and schemaClass, with no attributes -- which is all
	// curator-tool-ws needs to link to an instance that is already stored. An instance that is not stored yet is
	// returned as it is, since a commit needs its attributes to create it.
	private SimpleInstance createReference(SimpleInstance instance) {
		if (instance.getDbId() == null || instance.getDbId() < 0) {
			return instance;
		}
		SimpleInstance reference = new SimpleInstance();
		reference.setDbId(instance.getDbId());
		reference.setDisplayName(instance.getDisplayName());
		reference.setSchemaClassName(instance.getSchemaClassName());
		return reference;
	}

	// Checks if the instanceToCheck already contains the instanceToUse in the multi-value attribute
	@SuppressWarnings("unchecked")
	public SimpleInstance addAttributeValueIfNecessary(SimpleInstance instanceToBeCheckedForExistingAttribute,
														  SimpleInstance instanceContainingAttributeToBeChecked,
														  String attribute) {
		// Original version of this function had two checks: For 'multivalue attribute' and for 'instance-type object'. 
		// Now we know the only attributes being passed through here are inferredTo, inferredFrom, orthologousEvent,
		// and hasEvent, which are all multivalue attribute classes.
		// We also know that it will always be a GKInstance passed through here (see arguments), so we are able to
		// forego both checks.
		if (!hasAttributeValue(
			instanceToBeCheckedForExistingAttribute, instanceContainingAttributeToBeChecked, attribute)) {
			// The existing values are copied into a new List, which leaves the instance's own List unchanged until
			// it is set.
			List<SimpleInstance> attributeInstancesFromInferredInstance =
				getAttributeValues(instanceToBeCheckedForExistingAttribute, attribute);
			// Only a reference (dbId, displayName and schemaClass) is added, not the instance itself. The two
			// instances are often linked both ways (inferredFrom/inferredTo), and if each held the other with all its
			// attributes, committing either would hand curator-tool-ws a cycle that it walks until the stack overflows.
			attributeInstancesFromInferredInstance.add(createReference(instanceContainingAttributeToBeChecked));
			instanceToBeCheckedForExistingAttribute.setAttribute(attribute, attributeInstancesFromInferredInstance);
		}
		return instanceToBeCheckedForExistingAttribute;
	}
	
	// Returns the instance values of a multi-valued attribute as a new List. A SimpleInstance holds an unset attribute
	// as null and a single value as itself rather than in a List, so neither can be cast to a Collection.
	@SuppressWarnings("unchecked")
	public List<SimpleInstance> getAttributeValues(SimpleInstance instance, String attribute) {
		Object value = instance.getAttribute(attribute);
		if (value instanceof List) {
			return new ArrayList<>((List<SimpleInstance>) value);
		}
		return value != null ? new ArrayList<>(Collections.singletonList((SimpleInstance) value)) : new ArrayList<>();
	}

	// Caching. This function goes through each defining attribute of the incoming instance and produces a string of
	// the attribute values (DB IDs if the attribute is an instance).
	// This allows for identical instances held in memory to be used before trying to use fetchIdenticalInstances,
	// which is expensive.
	public String getCacheKey(SimpleInstance inferredInst) {
		StringBuilder key = new StringBuilder();
		for (String definingAttrName : getCuratorToolAPI().getDefiningAttributeNames(inferredInst.getSchemaClassName())) {
			Object attributeValue = inferredInst.getAttribute(definingAttrName);
			// A SimpleInstance holds a multi-valued attribute's values as a List and a single value as itself.
			if (attributeValue instanceof List) {
				List<?> multiValueAttributes = (List<?>) attributeValue;
				if (multiValueAttributes.isEmpty()) {
					key.append("null").append(CACHE_KEY_SEPARATOR);
				}
				for (Object attribute : multiValueAttributes) {
					key.append(getCacheKeyValue(attribute)).append(CACHE_KEY_SEPARATOR);
				}
			} else if (attributeValue != null) {
				key.append(getCacheKeyValue(attributeValue)).append(CACHE_KEY_SEPARATOR);
			} else {
				key.append("null").append(CACHE_KEY_SEPARATOR);
			}
		}
		return key.toString();
	}

	// Instance values are represented by their DB IDs, all others by their own value.
	private String getCacheKeyValue(Object attributeValue) {
		if (attributeValue instanceof SimpleInstance) {
			return String.valueOf(((SimpleInstance) attributeValue).getDbId());
		}
		return String.valueOf(attributeValue);
	}

	/**
	 * This method returns true if the only parent Pathway of the incoming Event instance is Disease.
	 * Instances with only Disease as a parent will not be inferred. Those that are a member of Disease AND
	 * another Pathway will be inferred, but the Disease Pathway (and sub-Pathway) inference will be skipped.
	 * @param eventInst -- GKInstance that will be checked for parent Pathways.
	 * @return boolean -- true if only parent is Disease Pathway, false if not.
	 * @throws Exception -- Thrown by MySQLAdaptor.
	 */
	public boolean onlyInDiseasePathway(SimpleInstance eventInst) throws Exception {
		Set<Long> topLevelPathwayDbIds = getTopLevelPathwayDbIds(eventInst);
		return topLevelPathwayDbIds.size() == 1 && topLevelPathwayDbIds.contains(DISEASE_PATHWAY_DB_ID);
	}

	/**
	 * Finds all 'hasEvent' referrals for the incoming Pathway instance. If there aren't any referrals via that
	 * attribute, then it is considered a TopLevelPathway. TopLevelPathway DbIds are added to a Set and returned.
	 * This method checks all Pathway referrals recursively up the Pathway hierarchy.
	 * @param pathway -- GKInstance that will be checked for 'hasEvent' referrals.
	 * @return -- Set<GKInstance> containing a list of all TopLevelPathway DbIds for incoming Event instance.
	 * @throws Exception -- Thrown by MySQLAdaptor.
	 */
	private Set<Long> getTopLevelPathwayDbIds(SimpleInstance pathway) throws Exception {
		List<SimpleInstance> parentPathways = safeList(getCuratorToolAPI().getReferrers(pathway, ReactomeJavaConstants.Event));
		if (parentPathways.isEmpty()) {
			return new HashSet<>(Arrays.asList(pathway.getDbId()));
		}

		Set<Long> topLevelPathwayDbIds = new HashSet<>();
		for (SimpleInstance parentPathway : parentPathways) {
			topLevelPathwayDbIds.addAll(getTopLevelPathwayDbIds(parentPathway));
		}

		return topLevelPathwayDbIds;
	}

	public SimpleInstance getEvidenceType() throws Exception {
		if (evidenceType == null) {
			String evidenceTypeText = "inferred by electronic annotation";

			evidenceType = new SimpleInstance();
			evidenceType.setSchemaClassName(ReactomeJavaConstants.EvidenceType);
			evidenceType.setDefaultPersonId(getCuratorToolAPI().getPersonId());
			evidenceType.setAttribute(ReactomeJavaConstants.name, Arrays.asList(evidenceTypeText, "IEA"));
			evidenceType.setDisplayName(evidenceTypeText);
			evidenceType = checkForIdenticalInstances(evidenceType, null);
		}
		return evidenceType;
	}

	public SimpleInstance getSummationInstance() throws Exception {
		if (summation == null) {
			summation = new SimpleInstance();
			summation.setSchemaClassName(ReactomeJavaConstants.Summation);
			summation.setDefaultPersonId(getCuratorToolAPI().getPersonId());
			String summationText = "This event has been computationally inferred from an event that has been" +
				" demonstrated in another species.<p>The inference is based on the homology mapping from PANTHER." +
				" Briefly, reactions for which all involved PhysicalEntities (in input, output and catalyst) have a" +
				" mapped orthologue/paralogue (for complexes at least 75% of components must have a mapping) are" +
				" inferred to the other species. High level events are also inferred for these events to allow for" +
				" easier navigation.<p><a href='/electronic_inference_compara.html' target = 'NEW'>More details and" +
				" caveats of the event inference in Reactome.</a> For details on PANTHER see also:" +
				" <a href='https://www.pantherdb.org/about.jsp' target='NEW'>https://www.pantherdb.org/about.jsp</a>";
			summation.setAttribute(ReactomeJavaConstants.text, summationText);
			summation.setDisplayName(summationText);
			summation = checkForIdenticalInstances(summation, null);
		}
		return summation;
	}

	public List<SimpleInstance> getCompartments(SimpleInstance instanceToBeInferred) throws Exception {
		List<SimpleInstance> compartments = new ArrayList<>();
		for (SimpleInstance compartmentInst : (List<SimpleInstance>) instanceToBeInferred.getAttribute(ReactomeJavaConstants.compartment)) {
			if (compartmentInst.getSchemaClassName().equals(ReactomeJavaConstants.Compartment)) {
				compartments.add(compartmentInst);
			} else {
				compartments.add(createCompartmentInstance(compartmentInst));
			}
		}
		return compartments;
	}

	public SimpleInstance getInflatedInstance(SimpleInstance instance, String attribute) {
		Object attributeValue = instance.getAttribute(attribute);
		return attributeValue != null ? curatorToolAPI.inflate((SimpleInstance) attributeValue) : null;
	}

	public List<SimpleInstance> getInflatedInstances(SimpleInstance instance, String attribute) {
		return safeList((Collection<SimpleInstance>) instance.getAttribute(attribute)).
			stream().map(i -> curatorToolAPI.inflate(i)).collect(Collectors.toList());
	}

	private String getSpeciesName() throws IOException, ParseException {
		return this.configProperties.getSpeciesConfig().getSpeciesName(this.targetSpeciesCode);
	}

	private long getPersonId() {
		return this.configProperties.getPersonId();
	}

	private CuratorToolAPI getCuratorToolAPI() {
		return this.curatorToolAPI;
	}

	private static <E> List<E> safeList(Collection<E> collection) {
		return Optional.ofNullable(collection).map(InstanceUtilities::convertToList).orElse(Collections.emptyList());
	}

	private static <E> List<E> convertToList(Collection<E> collection) {
		return new ArrayList<>(collection);
	}
}
