package org.reactome.orthoinference;

import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import static org.gk.model.ReactomeJavaConstants.*;

import org.gk.model.ReactomeJavaConstants;
import org.reactome.curation.model.SimpleInstance;
import org.springframework.stereotype.Component;

@Component
public class OrthologousEntityGenerator {
	
	private static final Logger logger = LogManager.getLogger();

	private EWASInferrer ewasInferrer;
	private InstanceUtilities instanceUtilities;
	private ProteinCountUtility proteinCountUtility;
	private SpeciesCheckUtility speciesCheckUtility;

	private SimpleInstance complexSummationInst;


	private static SimpleInstance nullInst = null;
	private static Map<SimpleInstance,SimpleInstance> orthologousEntityIdenticals = new HashMap<>();
	private static Map<SimpleInstance,SimpleInstance> homolEWASIdenticals = new HashMap<>();
	private static Map<SimpleInstance,SimpleInstance> complexPolymerIdenticals = new HashMap<>();
	private static Map<SimpleInstance,SimpleInstance> inferredEntitySetIdenticals = new HashMap<>();
	private static Map<String,SimpleInstance> definedSetIdenticals = new HashMap<>();
	private static Map<String,SimpleInstance> complexIdenticals = new HashMap<>();
	private static Map<String,SimpleInstance> entitySetIdenticals = new HashMap<>();



/** The heart of the OrthoInference process. This function takes PhysicalEntity (PE) instances and will infer those
 *  that are EWAS', Complexes/Polymers, or EntitySets.  The function's arguments are an incoming PE instance and an
 *  override attribute. Instances that are comprised of PE's will often recursively call this createOrthoEntity
 *  function on constituent PE's with the override attribute set to 'true'. This ensures that these PE's are inferred,
 *  despite the fact that they might not pass some filter criteria.  This is often handled using 'mock' instances
 *  (i.e. 'ghost instances' from Perl script), which allow a PE to be inferred without having to commit a 'real'
 *  instance to the DB.
*/
	public OrthologousEntityGenerator(
		EWASInferrer ewasInferrer,
		SpeciesCheckUtility speciesCheckUtility,
		InstanceUtilities instanceUtilities,
		ProteinCountUtility proteinCountUtility
	) throws Exception {
		this.ewasInferrer = ewasInferrer;
		this.speciesCheckUtility = speciesCheckUtility;
		this.instanceUtilities = instanceUtilities;
		this.proteinCountUtility = proteinCountUtility;
		this.complexSummationInst = getComplexSummationInstance();
	}

	public SimpleInstance createOrthoEntity(SimpleInstance entityInst, boolean override) throws Exception {
		logger.info("Attempting PE inference: " + entityInst);
		SimpleInstance infEntityInst = null;
		if (!isDrugOrOtherEntity(entityInst)) {
			// This used to have a conditional statement based on the returned value of the 'check_intracellular'
			// function.
			// That function doesn't exist anymore (only seemed to apply to the 'mtub' species, which hasn't been
			// inferred for a while).  Since the instance is species-agnostic, just returns the original instance.
			logger.info("Could not find valid species attribute, returning original instance: " + entityInst);
			return entityInst;
		}

		if (orthologousEntityIdenticals.get(entityInst) != null) {
			logger.info("Inferred PE instance already exists");
			return orthologousEntityIdenticals.get(entityInst);
		}

		// Checks that a species attribute exists in either the current instance or in constituent instances.
		if (!speciesCheckUtility.hasOrContainsSpeciesAttribute(entityInst)) {
			logger.info("No species attribute found in PE, using original instance");
			infEntityInst = entityInst;
		// Will either infer an EWAS or return a mock GEE instance if needed (i.e. if override is currently 'True')
		} else if (entityInst.getSchemaClassName().equals(ReactomeJavaConstants.GenomeEncodedEntity)) {
			// TODO: Try using 'isa' here instead of contains
			if (entityInst.getSchemaClassName().equals(ReactomeJavaConstants.EntityWithAccessionedSequence)) {
				infEntityInst = createInfEWAS(entityInst, override);
			} else {
				if (override) {
					logger.info("Mock GEE instance needed");
                    return instanceUtilities.createMockGKInstance(entityInst);
				}
			}
		// Infers Complex or Polymer instances -- Will recursively call createOrthoEntity with override on its
		// constituent PEs
		} else if (entityInst.getSchemaClassName().equals(ReactomeJavaConstants.Complex) ||
			entityInst.getSchemaClassName().equals(ReactomeJavaConstants.Polymer)) {
			infEntityInst = createInfComplexPolymer(entityInst, override);
		// Infers EntitySetInstances that themselves contain the species attribute (Not just constituent instances as
		// when hasSpecies is called above), returning the current instance if it doesn't.
		} else if (instanceUtilities.isSchemaClassType(entityInst, ReactomeJavaConstants.EntitySet)) {
			if (entityInst.getAttribute(ReactomeJavaConstants.species) != null) {
				infEntityInst = createInfEntitySet(entityInst, override);
			} else {
				logger.info("EntitySet has no species attribute, using original instance: " + entityInst);
				infEntityInst = entityInst;
			}
		// Handles SimpleEntities by returning the current instance. The idea behind this is that SimpleEntities
		// wouldn't need to be inferred since they wouldn't change between species {Note from infer_events.pl --
		// David Croft}.
		} else if (entityInst.getSchemaClassName().equals(ReactomeJavaConstants.SimpleEntity)) {
			logger.info("PE is a SimplyEntity, using original instance");
			infEntityInst = entityInst;
		} else {
			logger.warn("Unknown PhysicalEntity class: " + entityInst.getClass());
		}
		if (override) {
			return infEntityInst;
		}
		orthologousEntityIdenticals.put(entityInst, infEntityInst);
			logger.info("PE inference completed: " + entityInst);
			return infEntityInst;
	}
	
	// Function that first tries to infer any EWAS' associated with the instance. For those that have more than 1
	// returned EWAS instance, it's re-structured to a DefinedSet instance. If there is no EWAS instances inferred,
	// it will either return null or, if override is set, return a mock instance.
	private SimpleInstance createInfEWAS(SimpleInstance ewasInst, boolean override) throws Exception {
		if (homolEWASIdenticals.get(ewasInst) == null) {
			// Attempt to infer the EWAS 
			List<SimpleInstance> infEWASInstances = getEWASInferrer().inferEWAS(ewasInst);
			// If number of EWAS instances is greater than 1, then it is considered a DefinedSet. A new inferred
			// instance with definedSet class is created.
			if (infEWASInstances.size() > 1) {
				logger.info("Multiple EWAS homologues produced for single EWAS. Converting to DefinedSet");
				SimpleInstance infDefinedSetInst = new SimpleInstance();
				String definedSetName = "Homologues of " + ewasInst.getAttribute(name);
				infDefinedSetInst.setAttribute(name, definedSetName);
				infDefinedSetInst.setAttribute(ReactomeJavaConstants.compartment, instanceUtilities.getCompartments(ewasInst));
				
				infDefinedSetInst.setAttribute(ReactomeJavaConstants.species, getSpeciesInstance());
				infDefinedSetInst.setAttribute(ReactomeJavaConstants.hasMember, infEWASInstances);
				String definedSetDisplayName = infDefinedSetInst.getAttribute(ReactomeJavaConstants.name) +
					" [" +((SimpleInstance) ewasInst.getAttribute(compartment)).getDisplayName() + "]";
				infDefinedSetInst.setDisplayName(definedSetDisplayName);
				// Caching based on an instance's defining attributes. This reduces the number of
				// 'checkForIdenticalInstance' calls, which is slow.
				String cacheKey = instanceUtilities.getCacheKey(infDefinedSetInst);
				if (definedSetIdenticals.get(cacheKey) != null) {
					infDefinedSetInst = definedSetIdenticals.get(cacheKey);
				} else {
					infDefinedSetInst = instanceUtilities.checkForIdenticalInstances(infDefinedSetInst, ewasInst);
					definedSetIdenticals.put(cacheKey, infDefinedSetInst);
				}
				instanceUtilities.addInferredFromValue(infDefinedSetInst, ewasInst);
				instanceUtilities.addInferredToValue(ewasInst, infDefinedSetInst);
				homolEWASIdenticals.put(ewasInst, infDefinedSetInst);
				logger.info("Successfully converted to DefinedSet");
			} else if (infEWASInstances.size() == 1) {
				homolEWASIdenticals.put(ewasInst, infEWASInstances.get(0));
			} else {
				if (override) {
					logger.info("Mock EWAS instance needed");
					return instanceUtilities.createMockGKInstance(ewasInst);
				} else {
					return nullInst;
				}
			}
		} else {
			logger.info("Inferred EWAS already exists");
		}
		return homolEWASIdenticals.get(ewasInst);
	}

	private boolean isDrugOrOtherEntity(SimpleInstance entityInst) {
		String schemaClassName = entityInst.getSchemaClassName();
		return schemaClassName.contains(ReactomeJavaConstants.Drug) ||
			schemaClassName.equals(ReactomeJavaConstants.OtherEntity);
	}

	// Infers Complex or Polymer instances. These instances are generally comprised of more than 1 PhysicalEntity,
	// and calls 'createOrthoEntity' for each one. Complex/Polymer instances are also subject to the
	// 'countDistinctProteins' function. The result from this needs to have at least 75% of total proteins to be
	// inferrable for inference to continue.
	private SimpleInstance createInfComplexPolymer(SimpleInstance complexInst, boolean override)
		throws Exception {
		if (complexPolymerIdenticals.get(complexInst) == null) {
			List<Integer> complexProteinCounts = getProteinCountUtility().getDistinctProteinCounts(complexInst);
			int complexTotalProteinCounts = complexProteinCounts.get(0);
			int complexInferrableProteinCounts = complexProteinCounts.get(1);
			// int complexMax = complexProteinCounts.get(2); // Doesn't get used, since MaxHomologue isn't a valid
			// attribute anymore.

			// Filtering based on results of ProteinCounts and threshold (currently hard-coded at 75%).
			int percent = 0;
			if (complexTotalProteinCounts > 0) {
				percent = (complexInferrableProteinCounts * 100)/complexTotalProteinCounts;
			}
			if (!override) {
				if ((complexTotalProteinCounts > 0 && complexInferrableProteinCounts == 0) || percent < 75)
				{
					logger.info("Complex/Polymer protein count is below 75% threshold (" + percent + "%) -- " +
						"terminating inference");
					return nullInst;
				}
			}
			logger.info("Complex protein counts. Total: " + complexTotalProteinCounts +
				"  Inferrable: " + complexInferrableProteinCounts);
			SimpleInstance infComplexInst = instanceUtilities.createNewInferredGKInstance(complexInst);
			infComplexInst.setAttribute(ReactomeJavaConstants.summation, complexSummationInst);
			infComplexInst.setAttribute(ReactomeJavaConstants.name, complexInst.getAttribute(ReactomeJavaConstants.name));
			List<SimpleInstance> infComponentInstances = new ArrayList<>();
			// Inference handling is different depending on if it is a Complex or a Polymer. Complexes will infer all
			// 'components' while Polymers will infer all 'repeatedUnits'.
			// TODO: Log the ratio of inferred complex/polyer from human?
			if (complexInst.getSchemaClassName().equals(ReactomeJavaConstants.Complex)) {
				Collection<SimpleInstance> componentInstances = (Collection<SimpleInstance>) complexInst.getAttribute(ReactomeJavaConstants.hasComponent);
				logger.info("Complex components: " + componentInstances);
				for (SimpleInstance componentInst : componentInstances) {
					infComponentInstances.add(createOrthoEntity(componentInst, true));
				}
				infComplexInst.setAttribute(ReactomeJavaConstants.hasComponent, infComponentInstances);
			} else if (complexInst.getSchemaClassName().equals(ReactomeJavaConstants.Polymer)) {
				Collection<SimpleInstance> repeatedUnitInstances = (Collection<SimpleInstance>) complexInst.getAttribute(ReactomeJavaConstants.repeatedUnit);
				logger.info("Polymer repeated units: " + repeatedUnitInstances);
				for (SimpleInstance repeatedUnitInst : repeatedUnitInstances) {
					infComponentInstances.add(createOrthoEntity(repeatedUnitInst, true));
				}
				infComplexInst.setAttribute(repeatedUnit, infComponentInstances);
			} else {
				logger.warn(complexInst + " is not a Complex or a Polymer");
				return nullInst;
			}
			infComplexInst.setDisplayName(complexInst.getDisplayName());
			
			// Caching based on an instance's defining attributes. This reduces the number of
			// 'checkForIdenticalInstance' calls, which is slow.
			String cacheKey = instanceUtilities.getCacheKey(infComplexInst);
			if (complexIdenticals.get(cacheKey) != null) {
				infComplexInst = complexIdenticals.get(cacheKey);
			} else {
				infComplexInst = instanceUtilities.checkForIdenticalInstances(infComplexInst, complexInst);
				complexIdenticals.put(cacheKey, infComplexInst);
			}

			instanceUtilities.addInferredFromValue(infComplexInst, complexInst);
			instanceUtilities.addInferredToValue(complexInst, infComplexInst);
			
			if (override) {
				return infComplexInst;
			} 
			complexPolymerIdenticals.put(complexInst, infComplexInst);
		} else {
			logger.info("Inferred Complex/Polymer already exists");
		}
		return complexPolymerIdenticals.get(complexInst);
	}
	
	// EntitySet inference function. This function will initially call createOrthoEntity on all 'members' before
	// filtering by the type of EntitySet (Open, Candidate, or Defined Sets) and completing a specific inference.
	// Important to note is that while there are multiple cases where createOrthoEntity is called (for members and
	// candidates) in createInfEntitySet, the override functionality is not used here.  Presumably, this is because
	// the instances aren't a constituent part of a single instance (as in Complexes), but rather are stand-alone ones
	// that also happen to be included in a Set.  This means they should be subject  to the stringency of a typical
	// instance, rather then using override to create mock instances that allow an instance to be inferred more easily.
	@SuppressWarnings("unchecked")
	private SimpleInstance createInfEntitySet(SimpleInstance entitySetInst, boolean override)
		throws Exception {
		if (inferredEntitySetIdenticals.get(entitySetInst) == null) {
			// Equivalent to infer_members function in infer_events.pl
			Set<String> existingMemberInstances = new HashSet<>();
			List<SimpleInstance> infMembersList = new ArrayList<>();
			Collection<SimpleInstance> memberInstances =
				(Collection<SimpleInstance>) entitySetInst.getAttribute(ReactomeJavaConstants.hasMember);
			if (!entitySetInst.getSchemaClassName().equals(ReactomeJavaConstants.CandidateSet)) {
				logger.info("Total member instances: " + memberInstances.size());
				logger.info("Member instances: " + memberInstances);
			}
			for (SimpleInstance memberInst : memberInstances) {
				SimpleInstance infMemberInst = createOrthoEntity(memberInst, false);
				if (infMemberInst != null &&
					!existingMemberInstances.contains(infMemberInst.getAttribute(ReactomeJavaConstants.name).toString())) {

					existingMemberInstances.add(infMemberInst.getAttribute(ReactomeJavaConstants.name).toString());
					infMembersList.add(infMemberInst);
				}
			}
			if (!entitySetInst.getSchemaClassName().equals(ReactomeJavaConstants.CandidateSet)) {
				logger.info("Total number of inferred members: " +
					infMembersList.size() + "/" + memberInstances.size());
			}

			// Begin inference of EntitySet
			SimpleInstance infEntitySetInst = instanceUtilities.createNewInferredGKInstance(entitySetInst);
			infEntitySetInst.setAttribute(ReactomeJavaConstants.name, entitySetInst.getAttribute(ReactomeJavaConstants.name));
			infEntitySetInst.setAttribute(ReactomeJavaConstants.hasMember, infMembersList);

			// Begin specific inference process for each type of DefinedSet entity.
			List<Integer> entitySetProteinCounts = getProteinCountUtility().getDistinctProteinCounts(entitySetInst);
			int entitySetTotalCount = entitySetProteinCounts.get(0);
			int entitySetInferrableCount = entitySetProteinCounts.get(1);
			// int entitySetMax = entitySetProteinCounts.get(2);  // Doesn't get used, since MaxHomologue isn't a
			// valid attribute anymore
			
			// Filtering based on ProteinCount results
			if (!override && entitySetTotalCount > 0 && entitySetInferrableCount == 0) {
				logger.info("No distinct proteins found in EntitySet -- terminating inference");
				return nullInst;
			}
			
			if (instanceUtilities.isSchemaClassType(entitySetInst,ReactomeJavaConstants.CandidateSet)) {
				Set<String> existingCandidateInstances = new HashSet<>();
				List<SimpleInstance> infCandidatesList = new ArrayList<>();
				// Equivalent to infer_members function in infer_events.pl
				Collection<SimpleInstance> candidateInstances = (Collection<SimpleInstance>)
					entitySetInst.getAttribute(ReactomeJavaConstants.hasCandidate);
				logger.info("Total candidate instances: " + candidateInstances.size());
				logger.info("Candidate instances: " + candidateInstances);
				for (SimpleInstance candidateInst : candidateInstances) {
					SimpleInstance infCandidateInst = createOrthoEntity(candidateInst, false);
					if (infCandidateInst != null &&
						!existingMemberInstances.contains(infCandidateInst.getAttribute(ReactomeJavaConstants.name).toString()) &&
						!existingCandidateInstances.contains(infCandidateInst.getAttribute(ReactomeJavaConstants.name).toString())) {

						existingCandidateInstances.add(infCandidateInst.getAttribute(ReactomeJavaConstants.name).toString());
						infCandidatesList.add(infCandidateInst);
					}
				}
				logger.info("Total number of inferred candidates: " +
					infCandidatesList.size() + "/" + candidateInstances.size());
				// Handling of CandidateSets
				if (infCandidatesList.size() > 0) {
					infEntitySetInst.setAttribute(ReactomeJavaConstants.hasCandidate, infCandidatesList);
				} else {
					if (infMembersList.size() != 0) {
						if (infMembersList.size() == 1) {
							infEntitySetInst = infMembersList.get(0);
						} else {
							logger.info("No candidates inferred, but there are inferred members. " +
								"Converting to DefinedSet");
							//SchemaClass definedSetClass = dba.getSchema().getClassByName(DefinedSet);
							SimpleInstance infDefinedSetInst = new SimpleInstance();
							infDefinedSetInst.setSchemaClassName(ReactomeJavaConstants.DefinedSet);
							infDefinedSetInst.setAttribute(ReactomeJavaConstants.name, infEntitySetInst.getAttribute(ReactomeJavaConstants.name));
							infDefinedSetInst.setAttribute(ReactomeJavaConstants.hasMember, infMembersList);
							if (instanceUtilities.isValidAttribute(entitySetInst, ReactomeJavaConstants.compartment) &&
								entitySetInst.getAttribute(ReactomeJavaConstants.compartment) != null) {

								infDefinedSetInst.setAttribute(ReactomeJavaConstants.compartment, instanceUtilities.getCompartments(entitySetInst));

//								for (Object compartmentInst : entitySetInst.getAttribute(ReactomeJavaConstants.compartment)) {
//									GKInstance compartmentInstGk = (GKInstance) compartmentInst;
//									if (compartmentInstGk.getSchemClass().isa(Compartment)) {
//										infDefinedSetInst.addAttributeValue(compartment, compartmentInstGk);
//									} else {
//										GKInstance newCompartmentInst =
//											instanceUtilities.createCompartmentInstance(compartmentInstGk);
//										infDefinedSetInst.addAttributeValue(compartment, newCompartmentInst);
//									}
//								}
							}
							infDefinedSetInst.setAttribute(species, getSpeciesInstance());
							infEntitySetInst = infDefinedSetInst;
							logger.info("Successfully converted to DefinedSet");
						}
					} else {
						if (override) {
							logger.info("Mock CandidateSet instance needed");
							infEntitySetInst = instanceUtilities.createMockGKInstance(entitySetInst);
						} else {
							return nullInst;
						}
					}
				}	
			} else if (instanceUtilities.isSchemaClassType(entitySetInst, ReactomeJavaConstants.DefinedSet)) {
				if (infMembersList.size() == 0) {
					if (override) {
						logger.info("Mock DefinedSet instance needed");
						return instanceUtilities.createMockGKInstance(entitySetInst);
					} else {
						logger.info("No member instances found -- terminating inference");
						return nullInst;
					}
				} else if (infMembersList.size() == 1) {
					logger.info("Only 1 member from EntitySet was inferred, converting to PE: " +
						infMembersList.get(0));
					infEntitySetInst = infMembersList.get(0);
				}
				// If it has more than 1 member (which is the logic that would theoretically go here), nothing happens; 
				// All members are stored in this inferred instances 'hasMember' attribute near the beginning of this
				// function.
			}
			infEntitySetInst.setDisplayName(entitySetInst.getDisplayName());
			// Caching based on an instance's defining attributes. This reduces the number of
			// 'checkForIdenticalInstance' calls, which is slow.
			String cacheKey = instanceUtilities.getCacheKey(infEntitySetInst);
			if (entitySetIdenticals.get(cacheKey) != null) {
				infEntitySetInst = entitySetIdenticals.get(cacheKey);
			} else {
				infEntitySetInst = instanceUtilities.checkForIdenticalInstances(infEntitySetInst, entitySetInst);
				entitySetIdenticals.put(cacheKey, infEntitySetInst);
			}
			if (instanceUtilities.isValidAttribute(infEntitySetInst, ReactomeJavaConstants.species) &&
				entitySetInst.getAttribute(ReactomeJavaConstants.species) != null) {
				instanceUtilities.addInferredFromValue(infEntitySetInst, entitySetInst);
				instanceUtilities.addInferredToValue(entitySetInst, infEntitySetInst);
			}
			if (override) {
				return infEntitySetInst;
			}
			inferredEntitySetIdenticals.put(entitySetInst, infEntitySetInst);
		} else {
			logger.info("Inferred EntitySet already exists");
		}
		return inferredEntitySetIdenticals.get(entitySetInst);
	}

	public SimpleInstance getComplexSummationInstance() throws Exception {
		SimpleInstance complexSummationInst = new SimpleInstance();
		complexSummationInst.setSchemaClassName(ReactomeJavaConstants.Summation);
		String complexSummationText = "This complex/polymer has been computationally inferred (based on PANTHER) " +
			"from a complex/polymer involved in an event that has been demonstrated in another species.";
		complexSummationInst.setAttribute(text, complexSummationText);
		complexSummationInst.setDisplayName(complexSummationText);
		complexSummationInst = instanceUtilities.checkForIdenticalInstances(complexSummationInst, null);
		return complexSummationInst;
	}

	private EWASInferrer getEWASInferrer() {
		return this.ewasInferrer;
	}

	private ProteinCountUtility getProteinCountUtility() {
		return this.proteinCountUtility;
	}

	private SimpleInstance getSpeciesInstance() throws Exception {
		return this.instanceUtilities.getSpeciesInstance();
	}
}
