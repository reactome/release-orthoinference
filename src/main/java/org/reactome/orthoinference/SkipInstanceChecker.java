package org.reactome.orthoinference;

import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.ClassAttributeFollowingInstruction;
import org.gk.model.GKInstance;
import static org.gk.model.InstanceUtilities.followInstanceAttributes;

import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

@Component
public class SkipInstanceChecker {

	private static final Logger logger = LogManager.getLogger();
	private MySQLAdaptor dba;
	private InstanceUtilities instanceUtilities;
	private Set<Long> skipList;
	private final long HIV_INFECTION_DB_ID = 162906L;
	private final long INFLUENZA_INFECTION_DB_ID = 168255L;
	private final long AMYLOID_FIBER_FORMATION_DB_ID = 977225L;

	public SkipInstanceChecker(@Qualifier("currentDBA") MySQLAdaptor dba, InstanceUtilities instanceUtilities) {
		this.dba = dba;
		this.instanceUtilities = instanceUtilities;
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
	public boolean checkIfInstanceShouldBeSkipped(GKInstance reactionInst) throws Exception {
		// it is found in skiplist array
		if (getSkipList().contains(reactionInst.getDBID())) {
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
		if (reactionInst.getAttributeValue(ReactomeJavaConstants.isChimeric) != null) {
			if ((boolean) reactionInst.getAttributeValue(ReactomeJavaConstants.isChimeric)) {
				logger.info(reactionInst + " is chimeric -- skipping");
				return true;
			}
		}
		// it has related species
		if (reactionInst.getAttributeValue("relatedSpecies") != null) {
			logger.info(reactionInst + " has related species -- skipping");
			return true;
		}
		// it is a disease reaction
		if (reactionInst.getAttributeValue(ReactomeJavaConstants.disease) != null) {
			logger.info(reactionInst + " is a disease reaction -- skipping");
			return true;
		}
		// it is manually inferred
		if (reactionInst.getAttributeValue(ReactomeJavaConstants.inferredFrom) != null) {
			logger.info(reactionInst + " is manually inferred -- skipping");
			return true;
		}
		// it contains multiple species
		Collection<GKInstance> speciesInstances = checkIfEntitiesContainMultipleSpecies(reactionInst);
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
			GKInstance pathwayInst = dba.fetchInstance(pathwayIdToSkip);
			if (pathwayInst == null) {
				continue;
			}

			for (GKInstance reactionLikeEventToSkip : findReactionLikeEventInstancesToSkip(pathwayInst)) {
				eventDbIdsToSkip.add(reactionLikeEventToSkip.getDBID());
			}
		}
		return eventDbIdsToSkip;
	}

	private Set<GKInstance> findReactionLikeEventInstancesToSkip(GKInstance pathwayInst) throws Exception {
		// Finds all ReactionLikeEvents associated with the skiplists Pathway and hasEvent attributes, and
		// adds them to skiplist.
		List<ClassAttributeFollowingInstruction> classesToFollow = new ArrayList<>();
		classesToFollow.add(new ClassAttributeFollowingInstruction(
			ReactomeJavaConstants.Pathway, new String[]{ReactomeJavaConstants.hasEvent}, new String[]{}));
		String[] outClasses = new String[] {ReactomeJavaConstants.ReactionlikeEvent};
		return followInstanceAttributes(pathwayInst, classesToFollow, outClasses);
	}

	// Goes through all input/output/catalystActivity/regulatedBy attribute instances, and captures all species
	// associates with them. Returns a collection of species instances.
	@SuppressWarnings("unchecked")
	private Collection<GKInstance> checkIfEntitiesContainMultipleSpecies(GKInstance reactionInst)
		throws Exception {
		List<GKInstance> physicalEntityInstances = new ArrayList<>();
		physicalEntityInstances.addAll(getInputs(reactionInst));
		physicalEntityInstances.addAll(getOutputs(reactionInst));
		physicalEntityInstances.addAll(getCatalystActivityPhysicalEntities(reactionInst));
		physicalEntityInstances.addAll(getRegulators(reactionInst));

		physicalEntityInstances = removeDuplicates(physicalEntityInstances);
		physicalEntityInstances = addConstituentPhysicalEntities(physicalEntityInstances);
		return findUniqueSpeciesInstances(physicalEntityInstances);
	}

	private List<GKInstance> getInputs(GKInstance reactionLikeEvent) throws Exception {
		return reactionLikeEvent.getAttributeValuesList(ReactomeJavaConstants.input);
	}

	private List<GKInstance> getOutputs(GKInstance reactionLikeEvent) throws Exception {
		return reactionLikeEvent.getAttributeValuesList(ReactomeJavaConstants.output);
	}

	private List<GKInstance> getCatalystActivityPhysicalEntities(GKInstance reactionLikeEvent)
		throws Exception {

		List<GKInstance> catalystActivityPhysicalInstances = new ArrayList<>();
		for (GKInstance catalystActivityInst :
			(Collection<GKInstance>) reactionLikeEvent.getAttributeValuesList(ReactomeJavaConstants.catalystActivity)) {
			catalystActivityPhysicalInstances.addAll(catalystActivityInst.getAttributeValuesList(ReactomeJavaConstants.physicalEntity));
		}
		return catalystActivityPhysicalInstances;
	}

	private List<GKInstance> getRegulators(GKInstance reactionLikeEvent) throws Exception {
		List<GKInstance> regulatorInstances = new ArrayList<>();

		List<GKInstance> regulatedByInstances = reactionLikeEvent.getAttributeValuesList("regulatedBy");
		for (GKInstance regulatedByInst : regulatedByInstances) {
			for (GKInstance regulatorInst :
				(Collection<GKInstance>) regulatedByInst.getAttributeValuesList(ReactomeJavaConstants.regulator)) {
				if (regulatorInst.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity)) {
					regulatorInstances.add(regulatorInst);
				}
			}
		}

		return regulatorInstances;
	}

	private List<GKInstance> removeDuplicates(List<GKInstance> physicalEntities) {
		Map<String, GKInstance> physicalEntityHash = new HashMap<>();
		// Remove duplicates using HashMap
		for (GKInstance physicalEntityInst : physicalEntities) {
			physicalEntityHash.put(physicalEntityInst.getDBID().toString(), physicalEntityInst);
		}

		return new ArrayList<>(physicalEntityHash.values());
	}

	private List<GKInstance> addConstituentPhysicalEntities(List<GKInstance> physicalEntities)
		throws Exception {

		Map<String, GKInstance> physicalEntitiesDbIdToInstance = new HashMap<>();
		for (GKInstance physicalEntityInst : physicalEntities) {
			physicalEntitiesDbIdToInstance.put(physicalEntityInst.getDBID().toString(), physicalEntityInst);
			for (GKInstance constituentInst : recursePhysicalEntityConstituentInstances(physicalEntityInst)) {
				physicalEntitiesDbIdToInstance.put(constituentInst.getDBID().toString(), constituentInst);
			}
		}
		return new ArrayList<>(physicalEntitiesDbIdToInstance.values());
	}

	private List<GKInstance> findUniqueSpeciesInstances(List<GKInstance> physicalEntities) throws Exception {
		Map<String, GKInstance> speciesHash = new HashMap<>();
		for (GKInstance physicalEntityInst : physicalEntities) {
			if (physicalEntityInst.getSchemClass().isValidAttribute(ReactomeJavaConstants.species)) {
				for (GKInstance speciesInst :
					(Collection<GKInstance>) physicalEntityInst.getAttributeValuesList(ReactomeJavaConstants.species)) {
					speciesHash.put(speciesInst.getDBID().toString(), speciesInst);
				}
			}
		}
		return new ArrayList<>(speciesHash.values());
	}

	// Looks at referrals of the constituent instances for the species attribute as well
	// The term 'constituent' is used as a catch-all for instances under the hasMember/hasComponent/repeatedUnit
	// attributes
	@SuppressWarnings("unchecked")
	private Collection<GKInstance> recursePhysicalEntityConstituentInstances(GKInstance physicalEntity)
		throws Exception {
		List<String> attributes = Arrays.asList(
			ReactomeJavaConstants.hasMember, ReactomeJavaConstants.hasComponent,ReactomeJavaConstants.repeatedUnit);

		Map<Long, GKInstance> constituentInstances = new HashMap<>();
		for (String attribute : attributes) {
			if (physicalEntity.getSchemClass().isValidAttribute(attribute)) {
				for (GKInstance physicalEntityConstituent : (Collection<GKInstance>) physicalEntity.getAttributeValuesList(attribute)) {
					constituentInstances.put(physicalEntityConstituent.getDBID(), physicalEntityConstituent);
				}
			}
		}

		Map<Long, GKInstance> finalConstituentInstancesMap = new HashMap<>();
		for (GKInstance constituentInst : constituentInstances.values()) {
			finalConstituentInstancesMap.put(constituentInst.getDBID(), constituentInst);
			if (isEntitySet(constituentInst) || isComplex(constituentInst) || isPolymer(constituentInst)) {
				for (GKInstance recursedConstituentInst : recursePhysicalEntityConstituentInstances(constituentInst)) {
					finalConstituentInstancesMap.put(recursedConstituentInst.getDBID(), recursedConstituentInst);
				}
			}
		}
		return finalConstituentInstancesMap.values();
	}

	private boolean isEntitySet(GKInstance physicalEntity) {
		return isSchemaClassType(physicalEntity, ReactomeJavaConstants.EntitySet);
	}

	private boolean isComplex(GKInstance physicalEntity) {
		return isSchemaClassType(physicalEntity, ReactomeJavaConstants.Complex);
	}

	private boolean isPolymer(GKInstance physicalEntity) {
		return isSchemaClassType(physicalEntity, ReactomeJavaConstants.Polymer);
	}

	private boolean isSchemaClassType(GKInstance physicalEntity, String schemaClassType) {
		return physicalEntity.getSchemClass().isa(schemaClassType);
	}
}
