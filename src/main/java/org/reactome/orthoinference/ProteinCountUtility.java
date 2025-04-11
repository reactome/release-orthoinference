package org.reactome.orthoinference;

import java.util.*;

import org.gk.model.ClassAttributeFollowingInstruction;
import org.gk.model.GKInstance;
import org.gk.model.InstanceUtilities;
import org.gk.model.ReactomeJavaConstants;


public class ProteinCountUtility {
	
	private static Map<String, String[]> homologueMappings = new HashMap<>();
	
	/** This function is meant to emulate the count_distinct_proteins function found in infer_events.pl.
	 A crucial note is that the Perl version seems to be depend on the order by which instance groups are taken from
	 the DB. Often the DB IDs are ordered smallest to largest, but other times it is a consistent yet 'random' order.
	 What that means is that every time the Perl version will pull the instances out in the exact same order, but there
	 isn't a clear pattern (such as DB ID order) that is followed. Since this happens as well with the Java code,
	 sometimes the protein counts will differ from the Perl protein counts. The vast majority of the time this isn't
	 true, but this still suggests the protein count functionality should be re-written, for consistency's sake.
	 See the bottom of ProteinCount.checkCandidates for further elaboration. 
	*/
	
	public static List<Integer> getDistinctProteinCounts (GKInstance instanceToBeInferred) throws Exception {


		
		// With the output instances saved in followedInstances, begin the protein count process, which is based on
		// the homologue mappings (orthopairs) files.
		List<Integer> distinctProteinCounts = new ArrayList<>();
		int total = 0;
		int inferrable = 0;
		int max = 0;

		// Perform an AttributeQueryRequest with specified input attributes (ReactionlikeEvent, CatalystActivity,
		// Complex, Polymer, EWAS) and output attributes (ReferenceGeneProduct, EntitySet).
		Collection<GKInstance> sortedFollowedInstances = getSortedFollowedInstances(instanceToBeInferred);
		// If it is a ReferenceGene Product, the inferrable and max values are incremented depending on the number
		// of homologue mappings, while total is incremented for each entity.
		for (GKInstance entityInst : sortedFollowedInstances) {
			if (isReferenceGeneProduct(entityInst)) {
				int count = 0;
				String identifierName = entityInst.getAttributeValue(ReactomeJavaConstants.identifier).toString();
				if (homologueMappings.get(identifierName) != null) {
					count = homologueMappings.get(identifierName).length;
				}
				total++;
				if (count > max) {
					max = count;
				}
				if (count > 0) {
					inferrable++;
				}
			}
		}

		// For EntitySets, another AttributeQueryRequest is completed. This time the output classes are
		// Complex, Polymer, and ReferenceSequence.
		for (GKInstance entityInst : sortedFollowedInstances) {
			if (isEntitySet(entityInst)) {
				List<ClassAttributeFollowingInstruction> entitySetsInstancesToFollow = new ArrayList<ClassAttributeFollowingInstruction>();
				entitySetsInstancesToFollow.add(new ClassAttributeFollowingInstruction(
					ReactomeJavaConstants.DefinedSet, new String[]{ReactomeJavaConstants.hasMember}, new String[]{}));
				entitySetsInstancesToFollow.add(new ClassAttributeFollowingInstruction(
					ReactomeJavaConstants.CandidateSet, new String[]{ReactomeJavaConstants.hasMember}, new String[]{}));
				entitySetsInstancesToFollow.add(new ClassAttributeFollowingInstruction(
					ReactomeJavaConstants.EntityWithAccessionedSequence, new String[]{ReactomeJavaConstants.referenceEntity}, new String[]{}));

				String[] entitySetsOutClasses = new String[] {ReactomeJavaConstants.Complex, ReactomeJavaConstants.Polymer, ReactomeJavaConstants.ReferenceSequence};
				@SuppressWarnings("unchecked")
				Collection<GKInstance> entitySetsFollowedInstances = InstanceUtilities.followInstanceAttributes(
					entityInst, entitySetsInstancesToFollow, entitySetsOutClasses);
				Collection<GKInstance> entitySetsSortedInstances = sortInstancesByDbId(entitySetsFollowedInstances);
				
				if (entitySetsFollowedInstances.size() == 0 && entityInst.getSchemClass().isa(ReactomeJavaConstants.CandidateSet)) {
					// Protein counts are incremented depending on the number and types of candidates
					List<Integer> checkedCandidateInstances = getCandidateProteinCounts(
						entityInst, sortedFollowedInstances);
					if (checkedCandidateInstances.size() > 0) {
						total += checkedCandidateInstances.get(0);
						if (checkedCandidateInstances.size() > 1) {
							inferrable += checkedCandidateInstances.get(1);
						}
						if (checkedCandidateInstances.size() > 2) {
							max += checkedCandidateInstances.get(2);
						}
					}
					continue;
				}
				if (entitySetsFollowedInstances.size() > 0) {
					boolean uncountedInstances = false;
					// Little trick for breaking out of a nested loop
					outerloop:
					for (GKInstance physicalEntityInst : entitySetsSortedInstances) {
						for (GKInstance earlyFollowedInst : sortedFollowedInstances) {
							if (physicalEntityInst.getAttributeValue(ReactomeJavaConstants.DB_ID) ==
								earlyFollowedInst.getAttributeValue(ReactomeJavaConstants.DB_ID)) {
								continue outerloop;
							}
						}
						uncountedInstances = true;
					}
					if (!uncountedInstances) {
						continue;
					}
					// For Complexes and Polymers, the flag and flagInferred variables determine both if and by how
					// much the values are incremented.  These values are determined by the result of a recursive
					// ProteinCount call for each entity in the second followedInstances array.
					int flag = 0;
					int flagInferred = 0;
					for (GKInstance physicalEntityInst : entitySetsSortedInstances) {
						if (physicalEntityInst.getSchemClass().isa(ReactomeJavaConstants.Complex) ||
							physicalEntityInst.getSchemClass().isa(ReactomeJavaConstants.Polymer)) {
							List<Integer> complexProteinCounts = getDistinctProteinCounts(physicalEntityInst);
							if (complexProteinCounts != null) {
								int subTotal = complexProteinCounts.get(0);
								int subInferred = complexProteinCounts.get(1);
								int subMax = complexProteinCounts.get(2);
								if (subTotal > flag) {
									flag = subTotal;
								}
								if (subInferred > flagInferred) {
									flagInferred = subInferred;
								}
								if (subMax > max) {
									max = subMax;
								}
							}
						} else if (physicalEntityInst.getSchemClass().isa(ReactomeJavaConstants.ReferenceGeneProduct)) {
							flag = 1;
							String identifierName = physicalEntityInst.getAttributeValue(ReactomeJavaConstants.identifier).toString();
							int count = 0;
							if (homologueMappings.get(identifierName) != null) {
								count = homologueMappings.get(identifierName).length;
							}
							if (count > max) {
								max = count;
							}
							if (count > 0) {
								flagInferred = 1;
							}
						} 
					}
					// After going through the logic for Complexes/Polymers and ReferenceGeneProduct, the total and
					// inferrable values are incremented by their respective flag totals.
					total += flag;
					inferrable += flagInferred;
				}
			}
		}
		distinctProteinCounts.add(total);
		distinctProteinCounts.add(inferrable);
		distinctProteinCounts.add(max);
		return distinctProteinCounts;
	}

	private static Collection<GKInstance> getSortedFollowedInstances(GKInstance instanceToBeInferred)
		throws Exception {

		String[] outClasses = new String[] {ReactomeJavaConstants.ReferenceGeneProduct, ReactomeJavaConstants.EntitySet};
		@SuppressWarnings("unchecked")
		Collection<GKInstance> followedInstances = InstanceUtilities.followInstanceAttributes(
			instanceToBeInferred, getClassesToFollow(), outClasses);
		return sortInstancesByDbId(followedInstances);
	}

	// Function that determines protein counts of CandidateSets. Incoming arguments are the candidateSet of interest,
	// as well as the output array from the very first AttributeQueryRequest (AQR).
	// This 'output array from the first AQR' is used to prevent redundant counts, such as if a Candidate instance has
	// already undergone a protein count.
	private static List<Integer> getCandidateProteinCounts(
		GKInstance candidateSetInst, Collection<GKInstance> sortedFollowedInstances) throws Exception {
		List<Integer> checkedCandidateCounts = new ArrayList<>();
		if (hasCandidates(candidateSetInst)) {
			// AttributeQueryRequest for candidateSets, where the output instances are Complex, Polymer, and
			// ReferenceSequence
			int candidateTotal = 0;
			int candidateInferrable = 0;
			int candidateMax = 0;
			int flag = 0;
			List<ClassAttributeFollowingInstruction> candidateSetInstancesToFollow = new ArrayList<>();
			candidateSetInstancesToFollow.add(new ClassAttributeFollowingInstruction(
				ReactomeJavaConstants.CandidateSet, new String[]{ReactomeJavaConstants.hasCandidate}, new String[]{}));
			candidateSetInstancesToFollow.add(new ClassAttributeFollowingInstruction(
				ReactomeJavaConstants.EntityWithAccessionedSequence, new String[]{ReactomeJavaConstants.referenceEntity}, new String[]{}));
			String[] candidateSetOutClasses = new String[] {ReactomeJavaConstants.Complex, ReactomeJavaConstants.Polymer, ReactomeJavaConstants.ReferenceSequence};
			@SuppressWarnings("unchecked")
			Collection<GKInstance> candidateSetFollowedInstances = InstanceUtilities.followInstanceAttributes(
				candidateSetInst, candidateSetInstancesToFollow, candidateSetOutClasses);
			Collection<GKInstance> sortedCandidateSetFollowedInstances = sortInstancesByDbId(candidateSetFollowedInstances);

			boolean uncountedInstances = false;
			for (GKInstance physicalEntityInst : sortedCandidateSetFollowedInstances) {
				for (GKInstance earlyFollowedInst : sortedFollowedInstances) {
					if (physicalEntityInst.getAttributeValue(ReactomeJavaConstants.DB_ID) == earlyFollowedInst.getAttributeValue(ReactomeJavaConstants.DB_ID)) {
						continue;
					}
				}
				uncountedInstances = true;
			}
			if (!uncountedInstances) {
				return checkedCandidateCounts;
			}
			// For instances that are Complex or Polymer, the total, inferrable and max are incremented according to
			// the results of a recursive countDistinctProteins call.
			for (GKInstance physicalEntityInst : sortedCandidateSetFollowedInstances) {
				if (isComplex(physicalEntityInst) || isPolymer(physicalEntityInst)) {
					List<Integer> candidateComplexCounts = getDistinctProteinCounts(physicalEntityInst);
					if (candidateComplexCounts.size() > 0) {
						if (candidateComplexCounts.get(0) > 0 && candidateComplexCounts.get(1) == 0) {
							flag++;
						}
						if (candidateTotal > 0 && candidateComplexCounts.get(0) > candidateTotal) {
							candidateTotal = candidateComplexCounts.get(0);
						}
						if (candidateInferrable > 0 && candidateComplexCounts.get(1) > candidateInferrable) {
							candidateInferrable = candidateComplexCounts.get(1);
						}
						if (candidateMax > 0 && candidateComplexCounts.get(2) > candidateMax) {
							candidateMax = candidateComplexCounts.get(2);
						}
					}
					// ReferenceGeneProduct instances can only have an inferrable of 1 (So says the Perl version)
				} else if (physicalEntityInst.getSchemClass().isa(ReactomeJavaConstants.ReferenceGeneProduct)) {
					candidateTotal = 1;
					String identifierName = physicalEntityInst.getAttributeValue(ReactomeJavaConstants.identifier).toString();
					int count = 0;
					if (homologueMappings.get(identifierName) != null) {
						count = homologueMappings.get(identifierName).length;
					}
					if (count > 0) {
						candidateInferrable = 1;
					}
					if (candidateInferrable == 0) {
						flag++;
					}
				}
			} 
			// This was the tricky bit between Perl and Java. The flag is meant to drop any CandidateSet counts that
			// don't all have an inferrable protein.
			// In the Perl version, instance order makes a big difference here, since the flag value persists through
			// the for-loop.
			// Example: If a complex of 4 PEs has an inferrable in PE 1, 2, and 4, but not in 3, then the flag wouldn't
			// be raised (incorrect, I believe). Alternatively, if 1 doesn't have an inferrable, but 2, 3, 4 do, then
			// the flag would be raised (correctly). The second example given exemplifies how I believe the function is
			// meant to work, but the Perl version doesn't handle it. Currently this emulates the Perl version.
			if (flag > 0) {
				checkedCandidateCounts.add(candidateTotal);
				checkedCandidateCounts.add(0); // candidateInferred value is dropped, 0 is returned instead
				return checkedCandidateCounts;
			}
			checkedCandidateCounts.add(candidateTotal);
			checkedCandidateCounts.add(candidateInferrable);
			checkedCandidateCounts.add(candidateMax);
			return checkedCandidateCounts;
		}
		return checkedCandidateCounts;
	}

	private static List<ClassAttributeFollowingInstruction> getClassesToFollow() {
		Map<String, List<String>> classToAttributesToFollow = getClassToAttributesToFollow();

		List<ClassAttributeFollowingInstruction> classesToFollow = new ArrayList<>();
		for (String classToFollow : classToAttributesToFollow.keySet()) {
			List<String> attributesToFollow = classToAttributesToFollow.get(classToFollow);
			classesToFollow.add(new ClassAttributeFollowingInstruction(
				classToFollow, attributesToFollow.toArray(new String[0]), new String[]{}
			));
		}
		return classesToFollow;
	}

	private static Map<String, List<String>> getClassToAttributesToFollow() {
		Map<String, List<String>> classToAttributesToFollow = new LinkedHashMap<>();

		classToAttributesToFollow.put(
			ReactomeJavaConstants.ReactionlikeEvent,
			Arrays.asList(ReactomeJavaConstants.input, ReactomeJavaConstants.output, ReactomeJavaConstants.catalystActivity)
		);
		classToAttributesToFollow.put(
			ReactomeJavaConstants.CatalystActivity,
			Collections.singletonList(ReactomeJavaConstants.physicalEntity)
		);
		classToAttributesToFollow.put(
			ReactomeJavaConstants.Complex,
			Collections.singletonList(ReactomeJavaConstants.hasComponent)
		);
		classToAttributesToFollow.put(
			ReactomeJavaConstants.Polymer,
			Collections.singletonList(ReactomeJavaConstants.repeatedUnit)
		);
		classToAttributesToFollow.put(
			ReactomeJavaConstants.EntityWithAccessionedSequence,
			Collections.singletonList(ReactomeJavaConstants.referenceEntity)
		);
		return classToAttributesToFollow;
	}

	private static boolean isReferenceGeneProduct(GKInstance entity) {
		return entity.getSchemClass().isa(ReactomeJavaConstants.ReferenceGeneProduct);
	}

	private static boolean isEntitySet(GKInstance entity) {
		return entity.getSchemClass().isa(ReactomeJavaConstants.EntitySet);
	}

	private static boolean isComplex(GKInstance physicalEntity) {
		return physicalEntity.getSchemClass().isa(ReactomeJavaConstants.Complex);
	}

	private static boolean isPolymer(GKInstance physicalEntity) {
		return physicalEntity.getSchemClass().isa(ReactomeJavaConstants.Polymer);
	}

	private static boolean hasCandidates(GKInstance candidateSet) throws Exception {
		return candidateSet.getAttributeValue(ReactomeJavaConstants.hasCandidate) != null;
	}

	private static Collection<GKInstance> sortInstancesByDbId(Collection<GKInstance> unsortedInstances) {
		// Sort instances by DB ID
		List<Long> dbIds = new ArrayList<>();
		Collection<GKInstance> sortedInstances = new ArrayList<>();
		Map<Long,GKInstance> dbIdToInstance = new HashMap<>();
		for (GKInstance instance : unsortedInstances) {
			dbIds.add(instance.getDBID());
			dbIdToInstance.put(instance.getDBID(), instance);
		}
		Collections.sort(dbIds);
		for (Long id : dbIds) {
			sortedInstances.add(dbIdToInstance.get(id));
		}

		return sortedInstances;
	}

	public static void setHomologueMappingFile(Map<String, String[]> homologueMappingsCopy) {
		homologueMappings = homologueMappingsCopy;
	}
}
