package org.reactome.orthoinference;

import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;
import static org.gk.model.ReactomeJavaConstants.*;

import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.*;

public class OrthologousEntityGenerator {
	
	private static final Logger logger = LogManager.getLogger();
	private static MySQLAdaptor dba;
	private static GKInstance instanceEditInst;
	private static GKInstance complexSummationInst;
	private static GKInstance speciesInst;
	private static GKInstance nullInst = null;
	private static Map<GKInstance, GKInstance> orthologousEntityIdenticals = new HashMap<>();
	private static Map<GKInstance, GKInstance> homolEWASIdenticals = new HashMap<>();
	private static Map<GKInstance, GKInstance> complexPolymerIdenticals = new HashMap<>();
	private static Map<GKInstance, GKInstance> inferredEntitySetIdenticals = new HashMap<>();
	private static Map<String,GKInstance> definedSetIdenticals = new HashMap<>();
	private static Map<String,GKInstance> complexIdenticals = new HashMap<>();
	private static Map<String,GKInstance> entitySetIdenticals = new HashMap<>();
	private static Map<GKInstance, Set<GKInstance>> nonHumanParticpants = new HashMap<>();
	private static Map<GKInstance, GKInstance> inferredSARSIdenticals = new HashMap<>();
	private static Map<GKInstance, GKInstance> humanComplexIdenticals = new HashMap<>();

	/** The heart of the OrthoInference process. This function takes PhysicalEntity (PE) instances and will infer those that are EWAS', Complexes/Polymers, or EntitySets.
	 The function's arguments are an incoming PE instance and an override attribute. Instances that are comprised of PE's will often recursively call this createOrthoEntity function
	 on constituent PE's with the override attribute set to 'true'. This ensures that these PE's are inferred, despite the fact that they might not pass some filter criteria.
	 This is often handled using 'mock' instances (i.e. 'ghost instances' from Perl script), which allow a PE to be inferred without having to commit a 'real' instance to the DB.
*/
	public static GKInstance createOrthoEntity(GKInstance entityInst, boolean override) throws Exception
	{
		logger.info("Attempting PE inference: " + entityInst);
		GKInstance infEntityInst = null;
		if (!entityInst.getSchemClass().isValidAttribute(species)) {
			// This used to have a conditional statement based on the returned value of the 'check_intracellular' function.
			// That function doesn't exist anymore (only seemed to apply to the 'mtub' species, which hasn't been inferred for a while).
			// Since the instance is species-agnostic, just returns the original instance.
			logger.info("Could not find valid species attribute, returning original instance: " + entityInst);
			return entityInst;
		}

		if (orthologousEntityIdenticals.get(entityInst) != null) {
			logger.info("Inferred PE instance already exists");
			return orthologousEntityIdenticals.get(entityInst);
		}

		GKInstance entitySpeciesInst = (GKInstance) entityInst.getAttributeValue(species);
		if (entitySpeciesInst != null && entitySpeciesInst.getDBID().equals(48887L)) {
			inferSARSParticipants(entityInst);
			return entityInst;
		}

		// Checks that a species attribute exists in either the current instance or in constituent instances.
		if (!SpeciesCheckUtility.checkForSpeciesAttribute(entityInst))
		{
			logger.info("No species attribute found in PE, using original instance");
			infEntityInst = entityInst;
		// Will either infer an EWAS or return a mock GEE instance if needed (i.e. if override is currently 'True')
		} else if (entityInst.getSchemClass().isa(GenomeEncodedEntity))
		{
			// TODO: Try using 'isa' here instead of contains
			if (entityInst.getSchemClass().toString().contains(EntityWithAccessionedSequence))
			{
				infEntityInst = createInfEWAS(entityInst, override);
			} else {
				if (override)
				{
					logger.info("Mock GEE instance needed");
					GKInstance mockedInst = InstanceUtilities.createMockGKInstance(entityInst);
					return mockedInst;
				}
			}
		// Infers Complex or Polymer instances -- Will recursively call createOrthoEntity with override on its constituent PEs
		} else if (entityInst.getSchemClass().isa(Complex) || entityInst.getSchemClass().isa(Polymer))
		{
			infEntityInst = createInfComplexPolymer(entityInst, override);
		// Infers EntitySetInstances that themselves contain the species attribute (Not just constituent instances as when hasSpecies is called above),
		// returning the current instance if it doesn't.
		} else if (entityInst.getSchemClass().isa(EntitySet))
		{
			if (entityInst.getAttributeValue(species) != null)
			{
				infEntityInst = createInfEntitySet(entityInst, override);
			} else {
				logger.info("EntitySet has no species attribute, using original instance: " + entityInst);
				infEntityInst = entityInst;
			}
		// Handles SimpleEntities by returning the current instance. The idea behind this is that SimpleEntities wouldn't need
		// to be inferred since they wouldn't change between species {Note from infer_events.pl -- David Croft}.
		} else if (entityInst.getSchemClass().isa(SimpleEntity))
		{
			logger.info("PE is a SimplyEntity, using original instance");
			infEntityInst = entityInst;
		} else {
			logger.warn("Unknown PhysicalEntity class: " + entityInst.getClass());
		}
		if (override)
		{
			return infEntityInst;
		}
		orthologousEntityIdenticals.put(entityInst, infEntityInst);
			logger.info("PE inference completed: " + entityInst);
			return infEntityInst;
	}

	private static GKInstance inferSARSParticipants(GKInstance entityInst) throws Exception {

		if (humanComplexIdenticals.get(entityInst) == null) {
			Set<GKInstance> containedInstances = getComplexEntitySetContainedInstances(entityInst);

			boolean hasContainedSARSInstance = false;
			for (GKInstance containedInst : containedInstances) {
				if (hasSARSSpecies(containedInst)) {
					hasContainedSARSInstance = true;
					if (inferredSARSIdenticals.get(containedInst) == null) {
						GKInstance inferredSARSEntityInst = createOrthoEntity(containedInst, false);
						inferredSARSIdenticals.put(containedInst, inferredSARSEntityInst);
					}
				}
			}

			if (hasContainedSARSInstance) {
				GKInstance copiedHumanComplex = InstanceUtilities.createNewInferredGKInstance(entityInst);
				for (SchemaAttribute complexAttr : (Collection<SchemaAttribute>) entityInst.getSchemClass().getAttributes()) {
					if (!complexAttr.getName().equals(authored)
							&& !complexAttr.getName().equals(created)
							&& !complexAttr.getName().equals(modified)
							&& !complexAttr.getName().equals(relatedSpecies)
							&& !complexAttr.getName().equals(disease)
							&& !complexAttr.getName().equals(reviewed)
							&& !complexAttr.getName().equals(inferredFrom)
							&& !complexAttr.getName().equals(inferredTo)
							&& !complexAttr.getName().equals(DB_ID)
							&& !complexAttr.getName().equals(stableIdentifier)
							&& !complexAttr.getName().equals(revised)
							&& !complexAttr.getName().equals(edited)
							&& !complexAttr.getName().equals(compartment)
							&& !complexAttr.getName().equals(species)) {

						if (entityInst.getAttributeValuesList(complexAttr).size() > 0) {
							for (Object attrValue : entityInst.getAttributeValuesList(complexAttr)) {
								copiedHumanComplex.addAttributeValue(complexAttr, attrValue);
							}
						}
					}
				}

				List<GKInstance> components = (List<GKInstance>) copiedHumanComplex.getAttributeValuesList(hasComponent);
				List<GKInstance> updatedComponents = new ArrayList<>();
				for (GKInstance component : components) {
					if (hasSARSSpecies(component)) {
						updatedComponents.add(inferredSARSIdenticals.get(component));
					} else {
						updatedComponents.add(component);
					}
				}
				copiedHumanComplex.setAttributeValue(hasComponent, updatedComponents);

				copiedHumanComplex = InstanceUtilities.checkForIdenticalInstances(copiedHumanComplex, entityInst);

				copiedHumanComplex = InstanceUtilities.addAttributeValueIfNecessary(copiedHumanComplex, entityInst, inferredFrom);
				dba.updateInstanceAttribute(copiedHumanComplex, inferredFrom);
				entityInst = InstanceUtilities.addAttributeValueIfNecessary(entityInst, copiedHumanComplex, inferredTo);
				dba.updateInstanceAttribute(entityInst, inferredTo);

				humanComplexIdenticals.put(entityInst, copiedHumanComplex);


				/////// This code was used for troubleshooting and to see how far down 'multi-species' instances went in the Complex/EntitySet hierarchy
//			for (String attr : complexAttrs) {
//				System.out.println(attr);
//			}
//			for (GKInstance containedInst : containedInstances) {

//				if (hasSARSSpecies(containedInst)) {
//////					System.out.println("\t" + containedInst);
////					Set<GKInstance> subContainedInstances = getComplexEntitySetContainedInstances(containedInst);
////					for (GKInstance subContainedInst : subContainedInstances) {
////						if (hasSARSSpecies(subContainedInst)) {
//////							System.out.println("\t\t" + subContainedInst);
//////							Set<GKInstance> subSubContainedInstances = getComplexEntitySetContainedInstances(subContainedInst);
//////							for (GKInstance subSubContainedInst : subSubContainedInstances) {
//////								if (hasSARSSpecies(subSubContainedInst)) {
////////									System.out.println("\t\t\t" + subSubContainedInst);
//////									Set<GKInstance> subSubSubContainedInstances = getComplexEntitySetContainedInstances(subSubContainedInst);
//////									for (GKInstance subSubSubContainedInst : subSubSubContainedInstances) {
//////										if (hasSARSSpecies(subSubSubContainedInst)) {
////////											System.out.println("\t\t\t\t" + subSubSubContainedInst);
//////										} else if (hasContainedSARSInstance(subSubSubContainedInst)) {
//////
//////										} else {
//////											System.out.println(subSubSubContainedInst.getAttributeValue(species) + "\t\t" + subSubSubContainedInst);
//////										}
//////									}
//////								} else if (hasContainedSARSInstance(subSubContainedInst)) {
//////
//////								}
//////							}
////						} else if (hasContainedSARSInstance(subContainedInst)) {
////
////						} else {
////							System.out.println(subContainedInst.getAttributeValue(species) + "\t\t" + subContainedInst);
////						}
////					}
//				} else if (hasContainedSARSInstance(containedInst)) {
////					Set<GKInstance> subContainedInstances = getComplexEntitySetContainedInstances(containedInst);
////					System.out.println(subContainedInstances.size());
////					for (GKInstance subContainedInst : subContainedInstances) {
////						System.out.println(subContainedInst);
////						if (hasSARSSpecies(subContainedInst)) {
////							System.out.println("\t\tTWOO: " + subContainedInst);
////						} else if (hasContainedSARSInstance(subContainedInst)) {
////							System.out.println("\t\t\tTEE: " + subContainedInst);
////						} else {
////							System.out.println("\t\t\t\t\t\t\tDUDDD: " + subContainedInst);
////						}
////					}
//				} else {
////					System.out.println("\t\t\t\t\tDUD: " + containedInst);
////					System.out.println(containedInst.getAttributeValue(species) + "\t\t" + containedInst);
//				}
//			}
				/////////////

			}

		}
		return humanComplexIdenticals.get(entityInst);
	}

	private static boolean hasSARSSpecies(GKInstance entityInst) throws Exception {
		if (entityInst.getSchemClass().isValidAttribute(species)) {
			GKInstance speciesInst = (GKInstance) entityInst.getAttributeValue(species);
			return speciesInst != null && speciesInst.getDBID().equals(9678119L);
		}
		return false;
	}

	private static boolean hasContainedSARSInstance(GKInstance subEntityInst) throws Exception {
		boolean hasContainedSARSInstance = false;
		for (GKInstance subContainedInst : getComplexEntitySetContainedInstances(subEntityInst)) {
			if (hasSARSSpecies(subContainedInst)) {
				hasContainedSARSInstance = true;
			}
		}
		return hasContainedSARSInstance;
	}

	private static Set<GKInstance> getComplexEntitySetContainedInstances(GKInstance entityInst) throws Exception {
		return org.gk.model.InstanceUtilities.getContainedInstances(entityInst,
				ReactomeJavaConstants.hasMember,
				ReactomeJavaConstants.hasCandidate,
				ReactomeJavaConstants.hasComponent,
				ReactomeJavaConstants.repeatedUnit
		);
	}

	// Function that first tries to infer any EWAS' associated with the instance. For those that have more than 1 returned EWAS instance, 
	// it's re-structured to a DefinedSet instance. If there is no EWAS instances inferred, it will either return null or, if override is set, return a mock instance. 
	private static GKInstance createInfEWAS(GKInstance ewasInst, boolean override) throws InvalidAttributeException, Exception
	{
		if (homolEWASIdenticals.get(ewasInst) == null)
		{
			// Attempt to infer the EWAS 
			List<GKInstance> infEWASInstances = EWASInferrer.inferEWAS(ewasInst);
			// If number of EWAS instances is greater than 1, then it is considered a DefinedSet. A new inferred instance with definedSet class is created.
			if (infEWASInstances.size() > 1)
			{	
				logger.info("Multiple EWAS homologues produced for single EWAS. Converting to DefinedSet");
				SchemaClass definedSetClass = dba.getSchema().getClassByName(DefinedSet);
				GKInstance infDefinedSetInst = new GKInstance(definedSetClass);
				infDefinedSetInst.setDbAdaptor(dba);
				infDefinedSetInst.addAttributeValue(created, instanceEditInst);
				String definedSetName = "Homologues of " + ewasInst.getAttributeValue(name);
				infDefinedSetInst.addAttributeValue(name, definedSetName);
				
				GKInstance compartmentInstGk = (GKInstance) ewasInst.getAttributeValue(compartment);
				if (compartmentInstGk.getSchemClass().isa(Compartment)) {
					infDefinedSetInst.addAttributeValue(compartment, ewasInst.getAttributeValue(compartment));
				} else {
					GKInstance newCompartmentInst = InstanceUtilities.createCompartmentInstance(compartmentInstGk);
					infDefinedSetInst.addAttributeValue(compartment, newCompartmentInst);
				}
				
				infDefinedSetInst.addAttributeValue(species, speciesInst);
				infDefinedSetInst.addAttributeValue(hasMember, infEWASInstances);
				String definedSetDisplayName = (String) infDefinedSetInst.getAttributeValue(name) + " [" +((GKInstance) ewasInst.getAttributeValue(compartment)).getDisplayName() + "]";
				infDefinedSetInst.setAttributeValue(_displayName, definedSetDisplayName);
				// Caching based on an instance's defining attributes. This reduces the number of 'checkForIdenticalInstance' calls, which is slow.
				String cacheKey = InstanceUtilities.getCacheKey((GKSchemaClass) infDefinedSetInst.getSchemClass(), infDefinedSetInst);
				if (definedSetIdenticals.get(cacheKey) != null)
				{
					infDefinedSetInst = definedSetIdenticals.get(cacheKey);
				} else {
					infDefinedSetInst = InstanceUtilities.checkForIdenticalInstances(infDefinedSetInst, ewasInst);
					definedSetIdenticals.put(cacheKey, infDefinedSetInst);
				}
				infDefinedSetInst = InstanceUtilities.addAttributeValueIfNecessary(infDefinedSetInst, ewasInst, inferredFrom);
				dba.updateInstanceAttribute(infDefinedSetInst, inferredFrom);
				ewasInst = InstanceUtilities.addAttributeValueIfNecessary(ewasInst, infDefinedSetInst, inferredTo);
				dba.updateInstanceAttribute(ewasInst, inferredTo);
				homolEWASIdenticals.put(ewasInst, infDefinedSetInst);
				logger.info("Successfully converted to DefinedSet");
			} else if (infEWASInstances.size() == 1)
			{
				homolEWASIdenticals.put(ewasInst, infEWASInstances.get(0));
			} else {
				if (override) 
				{
					logger.info("Mock EWAS instance needed");
					return InstanceUtilities.createMockGKInstance(ewasInst);
				} else {
					return nullInst;
				}
			}
		} else {
			logger.info("Inferred EWAS already exists");
		}
		return homolEWASIdenticals.get(ewasInst);
	}
	// Infers Complex or Polymer instances. These instances are generally comprised of more than 1 PhysicalEntity, and calls 'createOrthoEntity' for each one. Complex/Polymer instances
	// are also subject to the 'countDistinctProteins' function. The result from this needs to have at least 75% of total proteins to be inferrable for inference to continue. 
	private static GKInstance createInfComplexPolymer(GKInstance complexInst, boolean override) throws InvalidAttributeException, InvalidAttributeValueException, Exception
	{
		if (complexPolymerIdenticals.get(complexInst) == null)
		{
			List<Integer> complexProteinCounts = ProteinCountUtility.getDistinctProteinCounts(complexInst);
			int complexTotalProteinCounts = complexProteinCounts.get(0);
			int complexInferrableProteinCounts = complexProteinCounts.get(1);
//			int complexMax = complexProteinCounts.get(2); // Doesn't get used, since MaxHomologue isn't a valid attribute anymore.
			
			// Filtering based on results of ProteinCounts and threshold (currently hard-coded at 75%).
			int percent = 0;
			if (complexTotalProteinCounts > 0)
			{
				percent = (complexInferrableProteinCounts * 100)/complexTotalProteinCounts;
			}
			if (!override)
			{
				if ((complexTotalProteinCounts > 0 && complexInferrableProteinCounts == 0) || percent < 75)
				{
					logger.info("Complex/Polymer protein count is below 75% threshold (" + percent + "%) -- terminating inference");
					return nullInst;
				}
			}
			logger.info("Complex protein counts. Total: " + complexTotalProteinCounts + "  Inferrable: " + complexInferrableProteinCounts);
			GKInstance infComplexInst = InstanceUtilities.createNewInferredGKInstance(complexInst);
//			infComplexInst.addAttributeValue(summation, complexSummationInst);
			infComplexInst.addAttributeValue(name, complexInst.getAttributeValue(name));
			List<GKInstance> infComponentInstances = new ArrayList<>();
			// Inference handling is different depending on if it is a Complex or a Polymer. Complexes will infer all 'components' while Polymers will infer all 'repeatedUnits'.
			// TODO: Log the ratio of inferred complex/polyer from human?
			if (complexInst.getSchemClass().isa(Complex))
			{
				Collection<GKInstance> componentInstances = complexInst.getAttributeValuesList(hasComponent);
				logger.info("Complex components: " + componentInstances);
				for (GKInstance componentInst : componentInstances)
				{	
					infComponentInstances.add(createOrthoEntity(componentInst, true));
				}
				infComplexInst.addAttributeValue(hasComponent, infComponentInstances);
			} else  if (complexInst.getSchemClass().isa(Polymer))
			{
				Collection<GKInstance> repeatedUnitInstances = complexInst.getAttributeValuesList(repeatedUnit);
				logger.info("Polymer repeated units: " + repeatedUnitInstances);
				for (GKInstance repeatedUnitInst : repeatedUnitInstances)
				{		
					infComponentInstances.add(createOrthoEntity(repeatedUnitInst, true));
				}
				infComplexInst.addAttributeValue(repeatedUnit, infComponentInstances);
			} else {
				logger.warn(complexInst + " is not a Complex or a Polymer");
				return nullInst;
			}
			infComplexInst.setAttributeValue(_displayName, complexInst.getAttributeValue(_displayName));
			
			// Caching based on an instance's defining attributes. This reduces the number of 'checkForIdenticalInstance' calls, which is slow.
			String cacheKey = InstanceUtilities.getCacheKey((GKSchemaClass) infComplexInst.getSchemClass(), infComplexInst);
			if (complexIdenticals.get(cacheKey) != null)
			{
				infComplexInst = complexIdenticals.get(cacheKey);
			} else {
				infComplexInst = InstanceUtilities.checkForIdenticalInstances(infComplexInst, complexInst);
				complexIdenticals.put(cacheKey, infComplexInst);
			}

			infComplexInst = InstanceUtilities.addAttributeValueIfNecessary(infComplexInst, complexInst, inferredFrom);
			dba.updateInstanceAttribute(infComplexInst, inferredFrom);
			complexInst = InstanceUtilities.addAttributeValueIfNecessary(complexInst, infComplexInst, inferredTo);
			dba.updateInstanceAttribute(complexInst, inferredTo);
			
			if (override)
			{
				return infComplexInst;
			} 
			complexPolymerIdenticals.put(complexInst, infComplexInst);
		} else {
			logger.info("Inferred Complex/Polymer already exists");
		}
		return complexPolymerIdenticals.get(complexInst);
	}
	
	// EntitySet inference function. This function will initially call createOrthoEntity on all 'members' before filtering by the type of EntitySet (Open, Candidate, or Defined Sets) and completing a specific inference.
	// Important to note is that while there are multiple cases where createOrthoEntity is called (for members and candidates) in createInfEntitySet, the override functionality is not used here. 
	// Presumably, this is because the instances aren't a constituent part of a single instance (as in Complexes), but rather are stand-alone ones that also happen to be included in a Set. 
	// This means they should be subject  to the stringency of a typical instance, rather then using override to create mock instances that allow an instance to be inferred more easily.
	@SuppressWarnings("unchecked")
	private static GKInstance createInfEntitySet(GKInstance entitySetInst, boolean override) throws InvalidAttributeException, Exception
	{
		if (inferredEntitySetIdenticals.get(entitySetInst) == null)
		{
			// Equivalent to infer_members function in infer_events.pl
			Set<String> existingMemberInstances = new HashSet<>();
			List<GKInstance> infMembersList = new ArrayList<>();
			Collection<GKInstance> memberInstances = (Collection<GKInstance>) entitySetInst.getAttributeValuesList(hasMember);
			if (!entitySetInst.getSchemClass().isa(CandidateSet)) {
				logger.info("Total member instances: " + memberInstances.size());
				logger.info("Member instances: " + memberInstances);
			}
			for (GKInstance memberInst : memberInstances)
			{
				GKInstance infMemberInst = createOrthoEntity(memberInst, false);
				if (infMemberInst != null && !existingMemberInstances.contains(infMemberInst.getAttributeValue(name).toString()))
				{
					existingMemberInstances.add(infMemberInst.getAttributeValue(name).toString());
					infMembersList.add(infMemberInst);
				}
			}
			if (!entitySetInst.getSchemClass().isa(CandidateSet)) {
				logger.info("Total number of inferred members: " + infMembersList.size() + "/" + memberInstances.size());
			}

			// Begin inference of EntitySet
			GKInstance infEntitySetInst = InstanceUtilities.createNewInferredGKInstance(entitySetInst);
			infEntitySetInst.addAttributeValue(name, entitySetInst.getAttributeValuesList(name));
			infEntitySetInst.addAttributeValue(hasMember, infMembersList);

			// Begin specific inference process for each type of DefinedSet entity.
			List<Integer> entitySetProteinCounts = ProteinCountUtility.getDistinctProteinCounts(entitySetInst);
			int entitySetTotalCount = entitySetProteinCounts.get(0);
			int entitySetInferrableCount = entitySetProteinCounts.get(1);
//				int entitySetMax = entitySetProteinCounts.get(2);  // Doesn't get used, since MaxHomologue isn't a valid attribute anymore
			
			// Filtering based on ProteinCount results
			if (!override && entitySetTotalCount > 0 && entitySetInferrableCount == 0)
			{
				logger.info("No distinct proteins found in EntitySet -- terminating inference");
				return nullInst;
			}
			
			if (entitySetInst.getSchemClass().isa(CandidateSet))
			{
				Set<String> existingCandidateInstances = new HashSet<>();
				List<GKInstance> infCandidatesList = new ArrayList<>();
				// Equivalent to infer_members function in infer_events.pl
				Collection<GKInstance> candidateInstances = (Collection<GKInstance>) entitySetInst.getAttributeValuesList(hasCandidate);
				logger.info("Total candidate instances: " + candidateInstances.size());
				logger.info("Candidate instances: " + candidateInstances);
				for (GKInstance candidateInst : candidateInstances)
				{
					GKInstance infCandidateInst = createOrthoEntity(candidateInst, false);
					if (infCandidateInst != null && !existingMemberInstances.contains(infCandidateInst.getAttributeValue(name).toString()) && !existingCandidateInstances.contains(infCandidateInst.getAttributeValue(name).toString()))
					{
						existingCandidateInstances.add(infCandidateInst.getAttributeValue(name).toString());
						infCandidatesList.add(infCandidateInst);
					}
				}
				logger.info("Total number of inferred candidates: " + infCandidatesList.size() + "/" + candidateInstances.size());
				// Handling of CandidateSets
				if (infCandidatesList.size() > 0)
				{
					infEntitySetInst.addAttributeValue(hasCandidate, infCandidatesList);
				} else {
					if (infMembersList.size() != 0)
					{
						if (infMembersList.size() == 1)
						{
							infEntitySetInst = infMembersList.get(0);
						} else {
							logger.info("No candidates inferred, but there are inferred members. Converting to DefinedSet");
							SchemaClass definedSetClass = dba.getSchema().getClassByName(DefinedSet);
							GKInstance infDefinedSetInst = new GKInstance(definedSetClass);
							infDefinedSetInst.setDbAdaptor(dba);
							infDefinedSetInst.addAttributeValue(created, instanceEditInst);
							infDefinedSetInst.setAttributeValue(name, infEntitySetInst.getAttributeValuesList(name));
							infDefinedSetInst.setAttributeValue(hasMember, infMembersList);
							if (entitySetInst.getSchemClass().isValidAttribute(compartment) && entitySetInst.getAttributeValue(compartment) != null) 
							{
								for (Object compartmentInst : entitySetInst.getAttributeValuesList(compartment)) {
									GKInstance compartmentInstGk = (GKInstance) compartmentInst;
									if (compartmentInstGk.getSchemClass().isa(Compartment)) 
									{
										infDefinedSetInst.addAttributeValue(compartment, compartmentInstGk);
									} else {
										GKInstance newCompartmentInst = InstanceUtilities.createCompartmentInstance(compartmentInstGk);
										infDefinedSetInst.addAttributeValue(compartment, newCompartmentInst);
									}
								}
							}
							infDefinedSetInst.addAttributeValue(species, speciesInst);
							infEntitySetInst = infDefinedSetInst;
							logger.info("Successfully converted to DefinedSet");
						}
					} else {
						if (override)
						{
							logger.info("Mock CandidateSet instance needed");
							infEntitySetInst = InstanceUtilities.createMockGKInstance(entitySetInst);
						} else {
							return nullInst;
						}
					}
				}	
			} else if (entitySetInst.getSchemClass().isa(DefinedSet))
			{
				if (infMembersList.size() == 0)
				{
					if (override)
					{
						logger.info("Mock DefinedSet instance needed");
						return InstanceUtilities.createMockGKInstance(entitySetInst);
					} else {
						logger.info("No member instances found -- terminating inference");
						return nullInst;
					}
				} else if (infMembersList.size() == 1)
				{
					logger.info("Only 1 member from EntitySet was inferred, converting to PE: " + infMembersList.get(0));
					infEntitySetInst = infMembersList.get(0);
				}
				// If it has more than 1 member (which is the logic that would theoretically go here), nothing happens; 
				// All members are stored in this inferred instances 'hasMember' attribute near the beginning of this function.
			}
			infEntitySetInst.setAttributeValue(_displayName, entitySetInst.getAttributeValue(_displayName));
			// Caching based on an instance's defining attributes. This reduces the number of 'checkForIdenticalInstance' calls, which is slow.
			String cacheKey = InstanceUtilities.getCacheKey((GKSchemaClass) infEntitySetInst.getSchemClass(), infEntitySetInst);
			if (entitySetIdenticals.get(cacheKey) != null)
			{
				infEntitySetInst = entitySetIdenticals.get(cacheKey);
			} else {
				infEntitySetInst = InstanceUtilities.checkForIdenticalInstances(infEntitySetInst, entitySetInst);
				entitySetIdenticals.put(cacheKey, infEntitySetInst);
			}
			if (infEntitySetInst.getSchemClass().isValidAttribute(species) && entitySetInst.getAttributeValue(species) != null)
			{
				infEntitySetInst = InstanceUtilities.addAttributeValueIfNecessary(infEntitySetInst, entitySetInst, inferredFrom);
				dba.updateInstanceAttribute(infEntitySetInst, inferredFrom);
				entitySetInst = InstanceUtilities.addAttributeValueIfNecessary(entitySetInst, infEntitySetInst, inferredTo);
				dba.updateInstanceAttribute(entitySetInst, inferredTo);
			}
			if (override)
			{
			return infEntitySetInst;
			}
			inferredEntitySetIdenticals.put(entitySetInst, infEntitySetInst);
		} else {
			logger.info("Inferred EntitySet already exists");
		}
		return inferredEntitySetIdenticals.get(entitySetInst);
	}
	
	public static void setAdaptor(MySQLAdaptor dbAdaptor)
	{
		dba = dbAdaptor;
	}
	
	public static void setSpeciesInstance(GKInstance speciesInstCopy)
	{
		speciesInst = speciesInstCopy;
	}
	
	public static void setInstanceEdit(GKInstance instanceEditCopy) 
	{
		instanceEditInst = instanceEditCopy;
	}
	
	public static void setComplexSummationInstance() throws Exception
	{
		complexSummationInst = new GKInstance(dba.getSchema().getClassByName(Summation));
		complexSummationInst.setDbAdaptor(dba);
		complexSummationInst.addAttributeValue(created, instanceEditInst);
		String complexSummationText = "This complex/polymer has been computationally inferred (based on PANTHER) from a complex/polymer involved in an event that has been demonstrated in another species.";
		complexSummationInst.addAttributeValue(text, complexSummationText);
		complexSummationInst.setAttributeValue(_displayName, complexSummationText);
		complexSummationInst = InstanceUtilities.checkForIdenticalInstances(complexSummationInst, null);
	}

	public static Map<GKInstance, Set<GKInstance>> getNonHumanParticipants() {
		return nonHumanParticpants;
	}
}
