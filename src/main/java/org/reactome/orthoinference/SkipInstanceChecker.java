package org.reactome.orthoinference;

import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.ReactomeJavaConstants;
import org.reactome.curation.model.SimpleInstance;
import org.springframework.stereotype.Component;

import static org.reactome.util.general.CollectionUtils.safeList;

@Component
public class SkipInstanceChecker {

	private static final Logger logger = LogManager.getLogger();
	private InstanceUtilities instanceUtilities;
	private CuratorToolAPI curatorToolAPI;

	private Set<Long> skipList;
	private final long HIV_INFECTION_DB_ID = 162906L;
	private final long INFLUENZA_INFECTION_DB_ID = 168255L;
	private final long AMYLOID_FIBER_FORMATION_DB_ID = 977225L;

	public SkipInstanceChecker(InstanceUtilities instanceUtilities, CuratorToolAPI curatorToolAPI) {
		this.instanceUtilities = instanceUtilities;
		this.curatorToolAPI = curatorToolAPI;

        try {
            this.skipList = getAllReactionLikeEventDbIdsToSkip();
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

	public Set<Long> getSkipList() {
		return this.skipList;
	}

	// Skip orthoinference of this instance if:
	public boolean instanceShouldBeSkipped(SimpleInstance reactionInst) throws Exception {
		// it is found in skiplist array
		if (getSkipList().contains(reactionInst.getDbId())) {
			logger.info(reactionInst + " is in skipList -- skipping");
			return true;
		}

		// If the only TopLevelPathway of a Reaction is 'Disease', then it is skipped.
		// Otherwise, it is inferred, making sure in cases where a Reaction is also a part of 'Disease' that that
		// Pathway is not inferred.
		if (instanceUtilities.onlyInDiseasePathway(reactionInst)) {
			logger.info(reactionInst + " has only Disease TopLevelPathway -- skipping");
			return true;
		}
		// it is chimeric
		if (reactionInst.getAttribute(ReactomeJavaConstants.isChimeric) != null) {
			if ((boolean) reactionInst.getAttribute(ReactomeJavaConstants.isChimeric)) {
				logger.info(reactionInst + " is chimeric -- skipping");
				return true;
			}
		}
		// it has related species
		if (reactionInst.getAttribute("relatedSpecies") != null) {
			logger.info(reactionInst + " has related species -- skipping");
			return true;
		}
		// it is a disease reaction
		if (reactionInst.getAttribute(ReactomeJavaConstants.disease) != null) {
			logger.info(reactionInst + " is a disease reaction -- skipping");
			return true;
		}
		// it is manually inferred
		if (reactionInst.getAttribute(ReactomeJavaConstants.inferredFrom) != null) {
			logger.info(reactionInst + " is manually inferred -- skipping");
			return true;
		}
		// it contains multiple species
		Collection<SimpleInstance> speciesInstances = checkIfEntitiesContainMultipleSpecies(reactionInst);
		if (speciesInstances.size() > 1) {
			logger.info(reactionInst + " has multiple species -- skipping");
			return true;
		}
		return false;
	}

	// Skiplist was traditionally provided in a file, but since it's currently just 3 instances, I've just hard-coded
	// them here.
	private List<Long> getPathwayDbIdsToSkip() {
		return Arrays.asList(
			HIV_INFECTION_DB_ID,
			INFLUENZA_INFECTION_DB_ID,
			AMYLOID_FIBER_FORMATION_DB_ID
		);
	}

	private Set<Long> getAllReactionLikeEventDbIdsToSkip() throws Exception {
		Set<Long> eventDbIdsToSkip = new HashSet<>();
		for (long pathwayIdToSkip : getPathwayDbIdsToSkip()) {
			SimpleInstance pathwayInst = curatorToolAPI.findByDbId(pathwayIdToSkip);
			if (pathwayInst == null) {
				continue;
			}

			for (SimpleInstance reactionLikeEventToSkip : findReactionLikeEventInstancesToSkip(pathwayInst)) {
				eventDbIdsToSkip.add(reactionLikeEventToSkip.getDbId());
			}
		}
		return eventDbIdsToSkip;
	}

	private Set<SimpleInstance> findReactionLikeEventInstancesToSkip(SimpleInstance pathwayInst) {
		// Finds all ReactionLikeEvents associated with the skiplists Pathway and hasEvent attributes, and
		// adds them to skiplist.
		Map<Long, SimpleInstance> reactionLikeEventsToSkip = new HashMap<>();
		findReactionLikeEventInstancesToSkip(pathwayInst, reactionLikeEventsToSkip, new HashSet<>());
		return new HashSet<>(reactionLikeEventsToSkip.values());
	}

	// Follows the hasEvent attribute of the Pathway, and of every Pathway under it, collecting the ReactionlikeEvents
	// found along the way. Pathways already visited are not followed again, in case hasEvent loops back on itself.
	private void findReactionLikeEventInstancesToSkip(
		SimpleInstance pathwayInst, Map<Long, SimpleInstance> reactionLikeEventsToSkip, Set<Long> visitedPathwayDbIds) {

		if (!visitedPathwayDbIds.add(pathwayInst.getDbId())) {
			return;
		}

		for (SimpleInstance eventInst : getAttributeValues(pathwayInst, ReactomeJavaConstants.hasEvent)) {
			if (instanceUtilities.isSchemaClassType(eventInst, "ReactionLikeEvent")) {
				reactionLikeEventsToSkip.put(eventInst.getDbId(), eventInst);
			} else if (instanceUtilities.isSchemaClassType(eventInst, ReactomeJavaConstants.Pathway)) {
				// Referred instances may only be shells, so the Pathway is inflated before its hasEvent is read
				SimpleInstance inflatedPathwayInst = curatorToolAPI.inflate(eventInst);
				if (inflatedPathwayInst != null) {
					findReactionLikeEventInstancesToSkip(
						inflatedPathwayInst, reactionLikeEventsToSkip, visitedPathwayDbIds);
				}
			}
		}
	}

	// Goes through all input/output/catalystActivity/regulatedBy attribute instances, and captures all species
	// associates with them. Returns a collection of species instances.
	@SuppressWarnings("unchecked")
	private Collection<SimpleInstance> checkIfEntitiesContainMultipleSpecies(SimpleInstance reactionInst)
		throws Exception {
		List<SimpleInstance> physicalEntityInstances = new ArrayList<>();
		physicalEntityInstances.addAll(getInputs(reactionInst));
		physicalEntityInstances.addAll(getOutputs(reactionInst));
		physicalEntityInstances.addAll(getCatalystActivityPhysicalEntities(reactionInst));
		physicalEntityInstances.addAll(getRegulators(reactionInst));

		physicalEntityInstances = removeDuplicates(physicalEntityInstances);
		physicalEntityInstances = addConstituentPhysicalEntities(physicalEntityInstances);
		return findUniqueSpeciesInstances(physicalEntityInstances);
	}

	// Referred instances are shells (dbId, displayName and schemaClass only), so the inputs and outputs are inflated
	// before their species and constituents are read.
	private List<SimpleInstance> getInputs(SimpleInstance reactionLikeEvent) {
		return inflateAll(safeList((List<SimpleInstance>) reactionLikeEvent.getAttribute(ReactomeJavaConstants.input)));
	}

	private List<SimpleInstance> getOutputs(SimpleInstance reactionLikeEvent) {
		return inflateAll(safeList((List<SimpleInstance>) reactionLikeEvent.getAttribute(ReactomeJavaConstants.output)));
	}

	private List<SimpleInstance> inflateAll(List<SimpleInstance> instances) {
		List<SimpleInstance> inflatedInstances = new ArrayList<>();
		for (SimpleInstance instance : instances) {
			SimpleInstance inflatedInst = curatorToolAPI.inflate(instance);
			if (inflatedInst != null) {
				inflatedInstances.add(inflatedInst);
			}
		}
		return inflatedInstances;
	}

	private List<SimpleInstance> getCatalystActivityPhysicalEntities(SimpleInstance reactionLikeEvent) {

		List<SimpleInstance> catalystActivityPhysicalInstances = new ArrayList<>();
		for (SimpleInstance catalystActivityInst : instanceUtilities.getInflatedInstances(reactionLikeEvent, ReactomeJavaConstants.catalystActivity)) {
			SimpleInstance physicalEntityInst =
				instanceUtilities.getInflatedInstance(catalystActivityInst, ReactomeJavaConstants.physicalEntity);
			if (physicalEntityInst != null) {
				catalystActivityPhysicalInstances.add(physicalEntityInst);
			}
		}
		return catalystActivityPhysicalInstances;
	}

	private List<SimpleInstance> getRegulators(SimpleInstance reactionLikeEvent) throws Exception {
		List<SimpleInstance> regulatorInstances = new ArrayList<>();

		List<SimpleInstance> regulatedByInstances = instanceUtilities.getInflatedInstances(reactionLikeEvent, "regulatedBy");
		for (SimpleInstance regulatedByInst : regulatedByInstances) {
			SimpleInstance regulatorInst = instanceUtilities.getInflatedInstance(regulatedByInst, ReactomeJavaConstants.regulator);
			if (regulatorInst != null &&
				instanceUtilities.isSchemaClassType(regulatorInst, ReactomeJavaConstants.PhysicalEntity)) {
				regulatorInstances.add(regulatorInst);
			}
		}

		return regulatorInstances;
	}

	private List<SimpleInstance> removeDuplicates(List<SimpleInstance> physicalEntities) {
		Map<String, SimpleInstance> physicalEntityHash = new HashMap<>();
		// Remove duplicates using HashMap
		for (SimpleInstance physicalEntityInst : physicalEntities) {
			physicalEntityHash.put(physicalEntityInst.getDbId().toString(), physicalEntityInst);
		}

		return new ArrayList<>(physicalEntityHash.values());
	}

	private List<SimpleInstance> addConstituentPhysicalEntities(List<SimpleInstance> physicalEntities)
		throws Exception {

		Map<String, SimpleInstance> physicalEntitiesDbIdToInstance = new HashMap<>();
		for (SimpleInstance physicalEntityInst : physicalEntities) {
			physicalEntitiesDbIdToInstance.put(physicalEntityInst.getDbId().toString(), physicalEntityInst);
			for (SimpleInstance constituentInst : recursePhysicalEntityConstituentInstances(physicalEntityInst)) {
				physicalEntitiesDbIdToInstance.put(constituentInst.getDbId().toString(), constituentInst);
			}
		}
		return new ArrayList<>(physicalEntitiesDbIdToInstance.values());
	}

	private List<SimpleInstance> findUniqueSpeciesInstances(List<SimpleInstance> physicalEntities) {
		Map<Long, SimpleInstance> speciesHash = new HashMap<>();
		for (SimpleInstance physicalEntityInst : physicalEntities) {
			for (SimpleInstance speciesInst : getAttributeValues(physicalEntityInst, ReactomeJavaConstants.species)) {
				speciesHash.put(speciesInst.getDbId(), speciesInst);
			}
		}
		return new ArrayList<>(speciesHash.values());
	}

	// Looks at referrals of the constituent instances for the species attribute as well
	// The term 'constituent' is used as a catch-all for instances under the hasMember/hasComponent/repeatedUnit
	// attributes
	private Collection<SimpleInstance> recursePhysicalEntityConstituentInstances(SimpleInstance physicalEntity)
		throws Exception {
		List<String> attributes = Arrays.asList(
			ReactomeJavaConstants.hasMember, ReactomeJavaConstants.hasComponent,ReactomeJavaConstants.repeatedUnit);

		Map<Long, SimpleInstance> constituentInstances = new HashMap<>();
		for (String attribute : attributes) {
			// Constituents are shells too, so each is inflated before its species and own constituents are read
			for (SimpleInstance physicalEntityConstituent : inflateAll(getAttributeValues(physicalEntity, attribute))) {
				constituentInstances.put(physicalEntityConstituent.getDbId(), physicalEntityConstituent);
			}
		}

		Map<Long, SimpleInstance> finalConstituentInstancesMap = new HashMap<>();
		for (SimpleInstance constituentInst : constituentInstances.values()) {
			finalConstituentInstancesMap.put(constituentInst.getDbId(), constituentInst);
			if (isEntitySet(constituentInst) || isComplex(constituentInst) || isPolymer(constituentInst)) {
				for (SimpleInstance recursedConstituentInst : recursePhysicalEntityConstituentInstances(constituentInst)) {
					finalConstituentInstancesMap.put(recursedConstituentInst.getDbId(), recursedConstituentInst);
				}
			}
		}
		return finalConstituentInstancesMap.values();
	}

	// Returns the attribute's values as a list, whether the attribute is single- or multi-valued.
	// Returns an empty list if the attribute is not valid for the instance's schema class or has no value.
	private List<SimpleInstance> getAttributeValues(SimpleInstance instance, String attributeName) {
		List<SimpleInstance> attributeValues = new ArrayList<>();
		if (!instanceUtilities.isValidAttribute(instance, attributeName)) {
			return attributeValues;
		}

		Object attributeValue = instance.getAttribute(attributeName);
		if (attributeValue instanceof SimpleInstance) {
			attributeValues.add((SimpleInstance) attributeValue);
		} else if (attributeValue instanceof Collection) {
			for (Object attributeValueObj : (Collection<?>) attributeValue) {
				attributeValues.add((SimpleInstance) attributeValueObj);
			}
		}
		return attributeValues;
	}

	private boolean isEntitySet(SimpleInstance physicalEntity) {
		return instanceUtilities.isSchemaClassType(physicalEntity, ReactomeJavaConstants.EntitySet);
	}

	private boolean isComplex(SimpleInstance physicalEntity) {
		return instanceUtilities.isSchemaClassType(physicalEntity, ReactomeJavaConstants.Complex);
	}

	private boolean isPolymer(SimpleInstance physicalEntity) {
		return instanceUtilities.isSchemaClassType(physicalEntity, ReactomeJavaConstants.Polymer);
	}
}
