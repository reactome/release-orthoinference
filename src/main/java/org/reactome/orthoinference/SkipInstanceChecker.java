package org.reactome.orthoinference;

import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.ClassAttributeFollowingInstruction;
import org.gk.model.GKInstance;
import static org.gk.model.InstanceUtilities.followInstanceAttributes;
import static org.gk.model.ReactomeJavaConstants.*;

import org.gk.persistence.MySQLAdaptor;

public class SkipInstanceChecker {

	private static final Logger logger = LogManager.getLogger();
	private static MySQLAdaptor dba;
	private static Set<String> skipList = new HashSet<>();
	private static final long HIV_INFECTION_DB_ID = 162906L;
	private static final long INFLUENZA_INFECTION_DB_ID = 168255L;
	private static final long AMYLOID_FIBER_FORMATION_DB_ID = 977225L;

	// Skiplist was traditionally provided in a file, but since it's currently just 3 instances, I've just hard-coded
	// them here.
	public static void buildStaticSkipList() throws Exception {
		List<Long> pathwayIdsToSkip = Arrays.asList(
			HIV_INFECTION_DB_ID,
			INFLUENZA_INFECTION_DB_ID,
			AMYLOID_FIBER_FORMATION_DB_ID
		);
		for (long pathwayId : pathwayIdsToSkip) {
			GKInstance pathwayInst = dba.fetchInstance(pathwayId);
			if (pathwayInst != null) {
				// Finds all ReactionLikeEvents associated with the skiplists Pathway and hasEvent attributes, and
				// adds them to skiplist.
				List<ClassAttributeFollowingInstruction> classesToFollow = new ArrayList<>();
				classesToFollow.add(new ClassAttributeFollowingInstruction(
					Pathway, new String[]{hasEvent}, new String[]{}));
				String[] outClasses = new String[] {ReactionlikeEvent};
				@SuppressWarnings("unchecked")
				Collection<GKInstance> followedInstances = followInstanceAttributes(
					pathwayInst, classesToFollow, outClasses);

				for (GKInstance entityInst : followedInstances) {
					skipList.add(entityInst.getDBID().toString());
				}
			}
		}
	}

	// Skip orthoinference of this instance if:
	public static boolean checkIfInstanceShouldBeSkipped(GKInstance reactionInst) throws Exception {
		// it is found in skiplist array
		if (skipList.contains(reactionInst.getDBID().toString())) {
			logger.info(reactionInst + " is in skipList -- skipping");
			return true;
		}

		// If the only TopLevelPathway of a Reaction is 'Disease', then it is skipped.
		// Otherwise, it is inferred, making sure in cases where a Reaction is also a part of 'Disease' that that
		// Pathway is not inferred.
		if (InstanceUtilities.onlyInDiseasePathway(reactionInst)) {
			logger.info(reactionInst + " has only Disease TopLevelPathway -- skipping");
			return true;
		}
		// it is chimeric
		if (reactionInst.getAttributeValue(isChimeric) != null) {
			if ((boolean) reactionInst.getAttributeValue(isChimeric)) {
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
		if (reactionInst.getAttributeValue(disease) != null) {
			logger.info(reactionInst + " is a disease reaction -- skipping");
			return true;
		}
		// it is manually inferred
		if (reactionInst.getAttributeValue(inferredFrom) != null) {
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

	// Goes through all input/output/catalystActivity/regulatedBy attribute instances, and captures all species
	// associates with them. Returns a collection of species instances.
	@SuppressWarnings("unchecked")
	private static Collection<GKInstance> checkIfEntitiesContainMultipleSpecies(GKInstance reactionInst)
		throws Exception {
		List<GKInstance> physicalEntityInstances = new ArrayList<>();
		physicalEntityInstances.addAll(reactionInst.getAttributeValuesList(input));
		physicalEntityInstances.addAll(reactionInst.getAttributeValuesList(output));
		for (GKInstance catalystActivityInst :
			(Collection<GKInstance>) reactionInst.getAttributeValuesList(catalystActivity)) {
			physicalEntityInstances.addAll(catalystActivityInst.getAttributeValuesList(physicalEntity));
		}
		List<GKInstance> regulatedByInstances =
			(ArrayList<GKInstance>) reactionInst.getAttributeValuesList("regulatedBy");

		if (regulatedByInstances != null) {
			for (GKInstance regulatedByInst : regulatedByInstances) {
				for (GKInstance regulatorInst :
					(Collection<GKInstance>) regulatedByInst.getAttributeValuesList(regulator)) {
					if (regulatorInst.getSchemClass().isa(PhysicalEntity)) {
						physicalEntityInstances.add(regulatorInst);
					}
				}
			}
		}
		Map<String, GKInstance> physicalEntityHash = new HashMap<>();
		// Remove duplicates using HashMap
		for (GKInstance physicalEntityInst : physicalEntityInstances) {
			physicalEntityHash.put(physicalEntityInst.getDBID().toString(), physicalEntityInst);
		}
		Map<String, GKInstance> physicalEntitiesFinal = new HashMap<>();
		for (GKInstance physicalEntityInst : physicalEntityHash.values()) {
			physicalEntitiesFinal.put(physicalEntityInst.getDBID().toString(), physicalEntityInst);
			Collection<GKInstance> allConstituentInstances =
				recursePhysicalEntityConstituentInstances(physicalEntityInst);
			if (allConstituentInstances != null) {
				for (GKInstance constituentInst : allConstituentInstances) {
					physicalEntitiesFinal.put(constituentInst.getDBID().toString(), constituentInst);
				}
			}
		}
		Map<String, GKInstance> speciesHash = new HashMap<>();
		for (GKInstance physicalEntityInst : physicalEntitiesFinal.values()) {
			if (physicalEntityInst.getSchemClass().isValidAttribute(species)) {
				for (GKInstance speciesInst :
					(Collection<GKInstance>) physicalEntityInst.getAttributeValuesList(species)) {
					speciesHash.put(speciesInst.getDBID().toString(), speciesInst);
				}
			}
		}
		return speciesHash.values();
	}

	// Looks at referrals of the constituent instances for the species attribute as well
	// The term 'constituent' is used as a catch-all for instances under the hasMember/hasComponent/repeatedUnit
	// attributes
	@SuppressWarnings("unchecked")
	private static Collection<GKInstance> recursePhysicalEntityConstituentInstances(GKInstance physicalEntity)
		throws Exception {
		Map<String, GKInstance> constituentInstances = new HashMap<>();
		if (physicalEntity.getSchemClass().isValidAttribute(hasMember)) {
			for (GKInstance memberInst : (Collection<GKInstance>) physicalEntity.getAttributeValuesList(hasMember)) {
				constituentInstances.put(memberInst.getDBID().toString(), memberInst);
			}
		}
		if (physicalEntity.getSchemClass().isValidAttribute(hasComponent)) {
			for (GKInstance componentInst :
				(Collection<GKInstance>) physicalEntity.getAttributeValuesList(hasComponent)) {
				constituentInstances.put(componentInst.getDBID().toString(), componentInst);
			}
		}
		if (physicalEntity.getSchemClass().isValidAttribute(repeatedUnit)) {
			for (GKInstance repeatedUnitInst :
				(Collection<GKInstance>) physicalEntity.getAttributeValuesList(repeatedUnit)) {
				constituentInstances.put(repeatedUnitInst.getDBID().toString(), repeatedUnitInst);
			}
		}
		if (constituentInstances.size() > 0) {
			Map<String, GKInstance> finalConstituentInstancesMap = new HashMap<>();
			for (GKInstance constituentInst : constituentInstances.values()) {
				finalConstituentInstancesMap.put(constituentInst.getDBID().toString(), constituentInst);
				if (constituentInst.getSchemClass().isa(EntitySet) ||
					constituentInst.getSchemClass().isa(Complex) ||
					constituentInst.getSchemClass().isa(Polymer)) {
					Collection<GKInstance> recursedConstituentInstances =
						recursePhysicalEntityConstituentInstances(constituentInst);
					if (recursedConstituentInstances != null) {
						for (GKInstance recursedConstituentInst : recursedConstituentInstances) {
							finalConstituentInstancesMap.put(
								recursedConstituentInst.getDBID().toString(), recursedConstituentInst);
						}
					}
				} else {
					continue;
				}
			}
			return finalConstituentInstancesMap.values();
		}
		return null;
	}

	public static void setAdaptor(MySQLAdaptor dbAdaptor) {
		dba = dbAdaptor;
	}
}
