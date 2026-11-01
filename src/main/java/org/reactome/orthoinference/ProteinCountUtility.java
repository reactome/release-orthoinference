package org.reactome.orthoinference;

import java.util.*;
import java.util.stream.Collectors;

import org.gk.model.ClassAttributeFollowingInstruction;
import org.gk.model.ReactomeJavaConstants;
import org.reactome.curation.model.SimpleInstance;
import org.reactome.server.graph.domain.model.DatabaseObject;
import org.springframework.stereotype.Component;

@Component
public class ProteinCountUtility {
	private Mappings mappings;
	private CuratorToolAPI curatorToolAPI;
	private InstanceUtilities instanceUtilities;

	public ProteinCountUtility(Mappings mappings, CuratorToolAPI curatorToolAPI, InstanceUtilities instanceUtilities) {
		this.mappings = mappings;
		this.curatorToolAPI = curatorToolAPI;
		this.instanceUtilities = instanceUtilities;
	}
	
	public List<Integer> getDistinctProteinCounts (SimpleInstance instanceToBeInferred) throws Exception {
		// With the output instances saved in followedInstances, begin the protein count process, which is based on
		// the homologue mappings (orthopairs) files.
		List<Integer> distinctProteinCounts = new ArrayList<>();
		int total = 0;
		int inferrable = 0;
		int max = 0;

		// Perform an AttributeQueryRequest with specified input attributes (ReactionlikeEvent, CatalystActivity,
		// Complex, Polymer, EWAS) and output attributes (ReferenceGeneProduct, EntitySet).
		Collection<SimpleInstance> sortedFollowedInstances = getSortedFollowedInstances(instanceToBeInferred);
		// If it is a ReferenceGene Product, the inferrable and max values are incremented depending on the number
		// of homologue mappings, while total is incremented for each entity.
		for (SimpleInstance entityInst : sortedFollowedInstances) {
			if (isReferenceGeneProduct(entityInst)) {
				int count = 0;
				String identifierName = entityInst.getAttribute(ReactomeJavaConstants.identifier).toString();
				if (getMappings().getHomologues(identifierName) != null) {
					count = getMappings().getHomologues(identifierName).length;
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
		for (SimpleInstance entityInst : sortedFollowedInstances) {
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
				Collection<SimpleInstance> entitySetsFollowedInstances = followInstanceAttributes(
					entityInst, entitySetsInstancesToFollow, entitySetsOutClasses);
				Collection<SimpleInstance> entitySetsSortedInstances = sortInstancesByDbId(entitySetsFollowedInstances);
				
				if (entitySetsFollowedInstances.size() == 0 && entityInst.getSchemaClassName().equals(ReactomeJavaConstants.CandidateSet)) {
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
					for (SimpleInstance physicalEntityInst : entitySetsSortedInstances) {
						for (SimpleInstance earlyFollowedInst : sortedFollowedInstances) {
							if (physicalEntityInst.getDbId().equals(earlyFollowedInst.getDbId())) {
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
					for (SimpleInstance physicalEntityInst : entitySetsSortedInstances) {
						if (physicalEntityInst.getSchemaClassName().equals(ReactomeJavaConstants.Complex) ||
							physicalEntityInst.getSchemaClassName().equals(ReactomeJavaConstants.Polymer)) {
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
						} else if (isReferenceGeneProduct(physicalEntityInst)) {
							flag = 1;
							String identifierName = physicalEntityInst.getAttribute(ReactomeJavaConstants.identifier).toString();
							int count = 0;
							if (getMappings().getHomologues(identifierName) != null) {
								count = getMappings().getHomologues(identifierName).length;
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

	private Collection<SimpleInstance> getSortedFollowedInstances(SimpleInstance instanceToBeInferred) {
		List<String> outClasses = Arrays.asList(ReactomeJavaConstants.ReferenceGeneProduct, ReactomeJavaConstants.EntitySet);
		Map<Long, SimpleInstance> followedInstances = new HashMap<>();
		followInstanceAttributes(
			instanceToBeInferred, getClassToAttributesToFollow(), outClasses, followedInstances, new HashSet<>());
		return followedInstances.values()
			.stream()
			.sorted(Comparator.comparing(SimpleInstance::getDbId))
			.collect(Collectors.toList());
	}

	// Takes the follow instructions and out classes in the form org.gk.model.InstanceUtilities.followInstanceAttributes
	// does, and returns the instances reached that are instances of one of the outClasses.
	@SuppressWarnings("unchecked")
	private Collection<SimpleInstance> followInstanceAttributes(
		SimpleInstance instance, List<ClassAttributeFollowingInstruction> classesToFollow, String[] outClasses) {

		Map<String, List<String>> classToAttributesToFollow = new LinkedHashMap<>();
		for (ClassAttributeFollowingInstruction classToFollow : classesToFollow) {
			classToAttributesToFollow
				.computeIfAbsent(classToFollow.getClassName(), className -> new ArrayList<>())
				.addAll((Collection<String>) classToFollow.getAttributes());
		}

		Map<Long, SimpleInstance> followedInstances = new HashMap<>();
		followInstanceAttributes(
			instance, classToAttributesToFollow, Arrays.asList(outClasses), followedInstances, new HashSet<>());
		return followedInstances.values();
	}

	// Equivalent of org.gk.model.InstanceUtilities.followInstanceAttributes: starting from the instance (which is
	// itself included), follows the attributes listed for each class it is an instance of, and collects every instance
	// reached that is an instance of one of the outClasses. Instances already visited are not followed again.
	private void followInstanceAttributes(
		SimpleInstance instance, Map<String, List<String>> classToAttributesToFollow, List<String> outClasses,
		Map<Long, SimpleInstance> followedInstances, Set<Long> visitedDbIds) {

		if (!visitedDbIds.add(instance.getDbId())) {
			return;
		}
		if (isAnySchemaClassType(instance, outClasses)) {
			followedInstances.put(instance.getDbId(), instance);
		}

		for (String classToFollow : classToAttributesToFollow.keySet()) {
			if (!getUtils().isSchemaClassType(instance, classToFollow)) {
				continue;
			}
			for (String attributeToFollow : classToAttributesToFollow.get(classToFollow)) {
				for (SimpleInstance referredInst : getAttributeValues(instance, attributeToFollow)) {
					// Only instances that are followed further or collected are of interest
					if (visitedDbIds.contains(referredInst.getDbId()) ||
						!(isAnySchemaClassType(referredInst, classToAttributesToFollow.keySet()) ||
						isAnySchemaClassType(referredInst, outClasses))) {
						continue;
					}
					// Referred instances are only shells, so they are inflated before their attributes are read
					SimpleInstance inflatedReferredInst = getCuratorToolAPI().inflate(referredInst);
					if (inflatedReferredInst != null) {
						followInstanceAttributes(
							inflatedReferredInst, classToAttributesToFollow, outClasses, followedInstances, visitedDbIds);
					}
				}
			}
		}
	}

	// Function that determines protein counts of CandidateSets. Incoming arguments are the candidateSet of interest,
	// as well as the output array from the very first AttributeQueryRequest (AQR).
	// This 'output array from the first AQR' is used to prevent redundant counts, such as if a Candidate instance has
	// already undergone a protein count.
	private List<Integer> getCandidateProteinCounts(
		SimpleInstance candidateSetInst, Collection<SimpleInstance> sortedFollowedInstances) throws Exception {
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
			Collection<SimpleInstance> candidateSetFollowedInstances = followInstanceAttributes(
				candidateSetInst, candidateSetInstancesToFollow, candidateSetOutClasses);
			Collection<SimpleInstance> sortedCandidateSetFollowedInstances = sortInstancesByDbId(candidateSetFollowedInstances);

			boolean uncountedInstances = false;
			for (SimpleInstance physicalEntityInst : sortedCandidateSetFollowedInstances) {
				for (SimpleInstance earlyFollowedInst : sortedFollowedInstances) {
					if (physicalEntityInst.getDbId().equals(earlyFollowedInst.getDbId())) {
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
			for (SimpleInstance physicalEntityInst : sortedCandidateSetFollowedInstances) {
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
				} else if (isReferenceGeneProduct(physicalEntityInst)) {
					candidateTotal = 1;
					String identifierName = physicalEntityInst.getAttribute(ReactomeJavaConstants.identifier).toString();
					int count = 0;
					if (getMappings().getHomologues(identifierName) != null) {
						count = getMappings().getHomologues(identifierName).length;
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

	private Map<String, List<String>> getClassToAttributesToFollow() {
		Map<String, List<String>> classToAttributesToFollow = new LinkedHashMap<>();

		classToAttributesToFollow.put(
			"ReactionLikeEvent",
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

	private boolean isReferenceGeneProduct(SimpleInstance entity) {
		return getUtils().isSchemaClassType(entity, ReactomeJavaConstants.ReferenceGeneProduct);
	}

	private boolean isEntitySet(SimpleInstance entity) {
		return getUtils().isSchemaClassType(entity, ReactomeJavaConstants.EntitySet);
	}

	private boolean isComplex(SimpleInstance physicalEntity) {
		return getUtils().isSchemaClassType(physicalEntity, ReactomeJavaConstants.Complex);
	}

	private boolean isPolymer(SimpleInstance physicalEntity) {
		return getUtils().isSchemaClassType(physicalEntity, ReactomeJavaConstants.Polymer);
	}

	private boolean hasCandidates(SimpleInstance candidateSet) throws Exception {
		return candidateSet.getAttribute(ReactomeJavaConstants.hasCandidate) != null;
	}

	private Collection<SimpleInstance> sortInstancesByDbId(Collection<SimpleInstance> unsortedInstances) {
		// Sort instances by DB ID
		List<Long> dbIds = new ArrayList<>();
		Collection<SimpleInstance> sortedInstances = new ArrayList<>();
		Map<Long,SimpleInstance> dbIdToInstance = new HashMap<>();
		for (SimpleInstance instance : unsortedInstances) {
			dbIds.add(instance.getDbId());
			dbIdToInstance.put(instance.getDbId(), instance);
		}
		Collections.sort(dbIds);
		for (Long id : dbIds) {
			sortedInstances.add(dbIdToInstance.get(id));
		}

		return sortedInstances;
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

	private boolean isAnySchemaClassType(SimpleInstance instance, Collection<String> schemaClassTypes) {
		return schemaClassTypes.stream().anyMatch(schemaClassType -> getUtils().isSchemaClassType(instance, schemaClassType));
	}

	private InstanceUtilities getUtils() {
		return this.instanceUtilities;
	}

	private Mappings getMappings() {
		return this.mappings;
	}

	private CuratorToolAPI getCuratorToolAPI() {
		return this.curatorToolAPI;
	}
}
