package org.reactome.orthoinference;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;
import static org.gk.model.ReactomeJavaConstants.*;
import static org.reactome.orthoinference.InstanceUtilities.inferDengueNameToZika;

import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.*;

public class OrthologousEntityGenerator {
	
	private static final Logger logger = LogManager.getLogger();
	private static final Pattern dengueNamePattern = Pattern.compile("dengue|denv[^a-zA-Z]?", Pattern.CASE_INSENSITIVE);

	private static MySQLAdaptor dba;
	private static GKInstance instanceEditInst;
	private static GKInstance complexSummationInst;
	private static GKInstance speciesInst;
	private static GKInstance nullInst = null;
	private static Map<GKInstance, GKInstance> orthologousEntityIdenticals = new HashMap<>();
	private static Map<GKInstance, GKInstance> homolEWASIdenticals = new HashMap<>();
	private static Map<GKInstance, GKInstance> complexPolymerIdenticals = new HashMap<>();
	private static Map<GKInstance, GKInstance> inferredEntitySetIdenticals = new HashMap<>();
	private static Map<GKInstance, GKInstance> inferredOtherEntityIdenticals = new HashMap<>();
	private static Map<String,GKInstance> definedSetIdenticals = new HashMap<>();
	private static Map<String,GKInstance> complexIdenticals = new HashMap<>();
	private static Map<String,GKInstance> entitySetIdenticals = new HashMap<>();
	private static Map<GKInstance, Set<GKInstance>> nonHumanParticpants = new HashMap<>();
	private static Map<GKInstance, GKInstance> inferredZikaIdenticals = new HashMap<>();
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
			if (dengueSpecificName(entityInst) && entityInst.getSchemClass().isa(OtherEntity)) {
				return createInfOtherEntity(entityInst);
			}

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
			if (entityInst.getSchemClass().isa(EntityWithAccessionedSequence) ||
				entityInst.getSchemClass().isa(GenomeEncodedEntity) ||
				!containsDengue(entityInst)) {
				return entityInst;
			}
			return inferZikaParticipants(entityInst);
		}

		// Checks that a species attribute exists in either the current instance or in constituent instances.
		if (!SpeciesCheckUtility.checkForSpeciesAttribute(entityInst) && !dengueSpecificName(entityInst))
		{
			logger.info("No species attribute found in " + entityInst + " - using original instance");
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
			if (entityInst.getAttributeValue(species) != null || dengueSpecificName(entityInst))
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

	private static boolean containsDengue(GKInstance entityInst) throws Exception {
		if (hasConstituents(entityInst)) {
			List<GKInstance> constituents = getConstituents(entityInst);
			return constituents.stream().anyMatch(constituent -> hasDengueSpecies(constituent));
		}
		return hasDengueSpecies(entityInst);
	}

	private static boolean dengueSpecificName(GKInstance entityInst) {
		String entityDisplayName = entityInst.getDisplayName().toLowerCase();

		Matcher matcher = dengueNamePattern.matcher(entityDisplayName);

		return matcher.find();
	}

	private static GKInstance inferZikaParticipants(GKInstance entityInst) throws Exception {

		if (humanComplexIdenticals.get(entityInst) == null) {
			Set<GKInstance> containedInstances = getComplexEntitySetContainedInstances(entityInst);

			boolean hasContainedDengueInstance = false;
			for (GKInstance containedInst : containedInstances) {
				if (containsDengue(containedInst)) {
					hasContainedDengueInstance = true;
					if (inferredZikaIdenticals.get(containedInst) == null) {
						GKInstance inferredZikaEntityInst = createOrthoEntity(containedInst, false);
						inferredZikaIdenticals.put(containedInst, inferredZikaEntityInst);
					}
				}
			}

			if (hasContainedDengueInstance) {
				// Outputs Human Complexes/EntitySets that contain Dengue instances.
//				System.out.println(entityInst);
				GKInstance copiedHumanComplexOrSet = InstanceUtilities.createNewInferredGKInstance(entityInst);
				for (SchemaAttribute complexOrSetAttr : (Collection<SchemaAttribute>) entityInst.getSchemClass().getAttributes()) {
					if (!complexOrSetAttr.getName().equals(authored)
							&& !complexOrSetAttr.getName().equals(created)
							&& !complexOrSetAttr.getName().equals(modified)
							&& !complexOrSetAttr.getName().equals(relatedSpecies)
							&& !complexOrSetAttr.getName().equals(disease)
							&& !complexOrSetAttr.getName().equals(reviewed)
							&& !complexOrSetAttr.getName().equals(inferredFrom)
							&& !complexOrSetAttr.getName().equals(inferredTo)
							&& !complexOrSetAttr.getName().equals(DB_ID)
							&& !complexOrSetAttr.getName().equals(stableIdentifier)
							&& !complexOrSetAttr.getName().equals(revised)
							&& !complexOrSetAttr.getName().equals(edited)
							&& !complexOrSetAttr.getName().equals(compartment)
							&& !complexOrSetAttr.getName().equals(species)) {

						if (entityInst.getAttributeValuesList(complexOrSetAttr).size() > 0) {
							for (Object attrValue : entityInst.getAttributeValuesList(complexOrSetAttr)) {
								if (complexOrSetAttr.getName().equals(name) || complexOrSetAttr.getName().equals(_displayName)) {
									attrValue = inferDengueNameToZika((String) attrValue);
								}
								copiedHumanComplexOrSet.addAttributeValue(complexOrSetAttr, attrValue);
							}
						}
					}
				}

				updateConstituents(copiedHumanComplexOrSet);

				copiedHumanComplexOrSet = InstanceUtilities.checkForIdenticalInstances(copiedHumanComplexOrSet, entityInst);

				copiedHumanComplexOrSet = InstanceUtilities.addAttributeValueIfNecessary(copiedHumanComplexOrSet, entityInst, inferredFrom);
				dba.updateInstanceAttribute(copiedHumanComplexOrSet, inferredFrom);
				entityInst = InstanceUtilities.addAttributeValueIfNecessary(entityInst, copiedHumanComplexOrSet, inferredTo);
				dba.updateInstanceAttribute(entityInst, inferredTo);

				humanComplexIdenticals.put(entityInst, copiedHumanComplexOrSet);


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

	public static boolean hasDengueSpecies(GKInstance entityInst) {
		final long dengueVirusType2SpeciesDbId = 3244621L;
		final long dengueVirusType2ThailandStrainSpeciesDbId = 9918331L;

		if (entityInst.getSchemClass().isValidAttribute(species)) {
			GKInstance speciesInst;
			try {
				speciesInst = (GKInstance) entityInst.getAttributeValue(species);
			} catch (Exception e) {
				throw new RuntimeException("Unable to get species instances from " + entityInst);
			}
			return speciesInst != null &&
				(speciesInst.getDBID().equals(dengueVirusType2SpeciesDbId) ||
				 speciesInst.getDBID().equals(dengueVirusType2ThailandStrainSpeciesDbId));
		}
		return false;
	}

	private static void updateConstituents(GKInstance complexOrSetInstance) throws Exception {
		if (complexOrSetInstance.getSchemClass().isa(ReactomeJavaConstants.Complex)) {
			updateComponents(complexOrSetInstance);
		} else if (complexOrSetInstance.getSchemClass().isa(ReactomeJavaConstants.EntitySet)) {
			updateMembers(complexOrSetInstance);
			if (complexOrSetInstance.getSchemClass().isa(ReactomeJavaConstants.CandidateSet)) {
				updateCandidates(complexOrSetInstance);
			}
		} else {
			throw new RuntimeException(complexOrSetInstance + " is not a complex or set");
		}
	}

	private static void updateComponents(GKInstance complex) throws Exception {
		List<GKInstance> components = getComponents(complex);
		List<GKInstance> updatedComponents = getUpdatedConstituents(components);
		complex.setAttributeValue(hasComponent, updatedComponents);
	}

	private static void updateMembers(GKInstance entitySet) throws Exception {
		List<GKInstance> members = getMembers(entitySet);
		List<GKInstance> updatedMembers = getUpdatedConstituents(members);
		entitySet.setAttributeValue(hasMember, updatedMembers);
	}

	private static void updateCandidates(GKInstance candidateSet) throws Exception {
		List<GKInstance> candidates = getCandidates(candidateSet);
		List<GKInstance> updatedCandidates = getUpdatedConstituents(candidates);
		candidateSet.setAttributeValue(hasCandidate, updatedCandidates);
	}

	private static boolean hasConstituents(GKInstance entityInst) {
		return entityInst.getSchemClass().isa(ReactomeJavaConstants.Complex) ||
			entityInst.getSchemClass().isa(ReactomeJavaConstants.EntitySet) ||
			entityInst.getSchemClass().isa(ReactomeJavaConstants.Polymer);
	}

	private static List<GKInstance> getConstituents(GKInstance entityInst) throws Exception {

		if (!hasConstituents(entityInst)) {
			return Collections.singletonList(entityInst);
		}

		List<GKInstance> constituents = new ArrayList<>();
		if (entityInst.getSchemClass().isa(Complex)) {
			for (GKInstance component : getComponents(entityInst)) {
				constituents.addAll(getConstituents(component));
			}
		} else if (entityInst.getSchemClass().isa(EntitySet)) {
			for (GKInstance member : getMembers(entityInst)) {
				constituents.addAll(getConstituents(member));
			}
			if (entityInst.getSchemClass().isa(CandidateSet)) {
				for (GKInstance candidate : getCandidates(entityInst)) {
					constituents.addAll(getConstituents(candidate));
				}
			}
		} else if (entityInst.getSchemClass().isa(Polymer)) {
			for (GKInstance repeatedUnit : getRepeatedUnits(entityInst)) {
				constituents.addAll(getConstituents(repeatedUnit));
			}
		}
		return constituents;
	}


	private static List<GKInstance> getComponents(GKInstance complex) throws Exception {
		return (List<GKInstance>) complex.getAttributeValuesList(hasComponent);
	}

	private static List<GKInstance> getMembers(GKInstance entitySet) throws Exception {
		return (List<GKInstance>) entitySet.getAttributeValuesList(hasMember);
	}

	private static List<GKInstance> getCandidates(GKInstance candidateSet) throws Exception {
		return (List<GKInstance>) candidateSet.getAttributeValuesList(hasCandidate);
	}

	private static List<GKInstance> getRepeatedUnits(GKInstance polymer) throws Exception {
		return (List<GKInstance>) polymer.getAttributeValuesList(repeatedUnit);
	}

	private static List<GKInstance> getUpdatedConstituents(List<GKInstance> constituents) throws Exception {
		List<GKInstance> updatedConstituents = new ArrayList<>();
		for (GKInstance constituent : constituents) {
			if (containsDengue(constituent)) {
				updatedConstituents.add(inferredZikaIdenticals.get(constituent));
			} else {
				updatedConstituents.add(constituent);
			}
		}
		return updatedConstituents;
	}

	private static boolean hasContainedDengueInstance(GKInstance subEntityInst) throws Exception {
		boolean hasContainedDengueInstance = false;
		for (GKInstance subContainedInst : getComplexEntitySetContainedInstances(subEntityInst)) {
			if (hasDengueSpecies(subContainedInst)) {
				hasContainedDengueInstance = true;
			}
		}
		return hasContainedDengueInstance;
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
				infDefinedSetInst.addAttributeValue(name, inferDengueNameToZika(definedSetName));
				
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
				infDefinedSetInst.setAttributeValue(_displayName, inferDengueNameToZika(definedSetDisplayName));
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
			logger.info("Inferred EWAS for " + ewasInst + " already exists");
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
			infComplexInst.addAttributeValue(name, inferDengueNameToZika((String) complexInst.getAttributeValue(name)));
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
			infComplexInst.setAttributeValue(_displayName, inferDengueNameToZika((String) complexInst.getAttributeValue(_displayName)));
			
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
			infEntitySetInst.addAttributeValue(name, inferDengueNameToZika(entitySetInst.getAttributeValuesList(name)));
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
							infDefinedSetInst.setAttributeValue(name, inferDengueNameToZika(infEntitySetInst.getAttributeValuesList(name)));
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
			infEntitySetInst.setAttributeValue(_displayName, inferDengueNameToZika((String) entitySetInst.getAttributeValue(_displayName)));
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

	private static GKInstance createInfOtherEntity(GKInstance entityInst) throws Exception
	{
		if (inferredOtherEntityIdenticals.containsKey(entityInst)) {
			return inferredOtherEntityIdenticals.get(entityInst);
		}

		GKInstance infOtherEntityInst = InstanceUtilities.createNewInferredGKInstance(entityInst);
		infOtherEntityInst.setDbAdaptor(dba);
		infOtherEntityInst.addAttributeValue(created, instanceEditInst);
		infOtherEntityInst.addAttributeValue(name, inferDengueNameToZika(entityInst.getAttributeValuesList(name)));
		infOtherEntityInst.setAttributeValue(_displayName, InstanceDisplayNameGenerator.generateDisplayName(infOtherEntityInst));

		if (infOtherEntityInst.getSchemClass().isValidAttribute(compartment) && infOtherEntityInst.getAttributeValue(compartment) != null)
		{
			for (Object compartmentInst : entityInst.getAttributeValuesList(compartment)) {
				GKInstance compartmentInstGk = (GKInstance) compartmentInst;
				if (compartmentInstGk.getSchemClass().isa(Compartment))
				{
					infOtherEntityInst.addAttributeValue(compartment, compartmentInstGk);
				} else {
					GKInstance newCompartmentInst = InstanceUtilities.createCompartmentInstance(compartmentInstGk);
					infOtherEntityInst.addAttributeValue(compartment, newCompartmentInst);
				}
			}
		}

		InstanceUtilities.addAttributeValueIfNecessary(infOtherEntityInst, entityInst, inferredFrom);
		//dba.updateInstanceAttribute(infOtherEntityInst, inferredFrom);
		entityInst = InstanceUtilities.addAttributeValueIfNecessary(entityInst, infOtherEntityInst, inferredTo);
		dba.updateInstanceAttribute(entityInst, inferredTo);



		inferredOtherEntityIdenticals.put(entityInst, infOtherEntityInst);
		return infOtherEntityInst;
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
