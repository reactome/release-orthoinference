package org.reactome.orthoinference;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.GKSchemaClass;
import org.gk.schema.SchemaClass;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

@Component
public class OrthologousEntityGenerator implements OrthoEntityInferrer {
	
	private static final Logger logger = LogManager.getLogger();

	private EntitySetInferrer entitySetInferrer;
	private ComplexPolymerInferrer complexPolymerInferrer;
	private EWASInferrer ewasInferrer;
	private InstanceUtilities instanceUtilities;
	private MySQLAdaptor dba;

	//private GKInstance complexSummationInst;


	private static GKInstance nullInst = null;
	private static Map<GKInstance, GKInstance> orthologousEntityIdenticals = new HashMap<>();
	private static Map<GKInstance, GKInstance> homolEWASIdenticals = new HashMap<>();
	private static Map<String,GKInstance> definedSetIdenticals = new HashMap<>();



/** The heart of the OrthoInference process. This function takes PhysicalEntity (PE) instances and will infer those
 *  that are EWAS', Complexes/Polymers, or EntitySets.  The function's arguments are an incoming PE instance and an
 *  override attribute. Instances that are comprised of PE's will often recursively call this createOrthoEntity
 *  function on constituent PE's with the override attribute set to 'true'. This ensures that these PE's are inferred,
 *  despite the fact that they might not pass some filter criteria.  This is often handled using 'mock' instances
 *  (i.e. 'ghost instances' from Perl script), which allow a PE to be inferred without having to commit a 'real'
 *  instance to the DB.
*/
	public OrthologousEntityGenerator(
		@Qualifier("currentDBA") MySQLAdaptor dba,
		ComplexPolymerInferrer complexPolymerInferrer,
		EntitySetInferrer entitySetInferrer,
		EWASInferrer ewasInferrer,
		InstanceUtilities instanceUtilities
	) throws Exception {
		this.dba = dba;
		this.complexPolymerInferrer = complexPolymerInferrer;
		this.entitySetInferrer = entitySetInferrer;
		this.ewasInferrer = ewasInferrer;
		this.instanceUtilities = instanceUtilities;
	}

	@Override
	public GKInstance createOrthoEntity(GKInstance entityInst, boolean override) throws Exception {
		logger.info("Attempting PE inference: " + entityInst);

		if (!hasSpeciesAttribute(entityInst)) {
			return entityInst;
		}

		if (orthologousEntityIdenticals.get(entityInst) != null) {
			logger.info("Inferred PE instance already exists");
			return orthologousEntityIdenticals.get(entityInst);
		}

		GKInstance infEntityInst = inferEntityBasedOnType(entityInst, override);

		if (override) {
			return infEntityInst;
		}

		orthologousEntityIdenticals.put(entityInst, infEntityInst);
		logger.info("PE inference completed: " + entityInst);
		return infEntityInst;
	}

	private GKInstance inferEntityBasedOnType(GKInstance entityInst, boolean override) throws Exception {
		// No species attribute case
		if (!SpeciesCheckUtility.hasOrContainsSpeciesAttribute(entityInst)) {
			logger.info("No species attribute found in PE, using original instance");
			return entityInst;
		}

		// Handle different entity types
		if (entityInst.getSchemClass().isa(ReactomeJavaConstants.GenomeEncodedEntity)) {
			return handleGenomeEncodedEntity(entityInst, override);
		}

		if (entityInst.getSchemClass().isa(ReactomeJavaConstants.Complex) ||
			entityInst.getSchemClass().isa(ReactomeJavaConstants.Polymer)) {
			return complexPolymerInferrer.createInfComplexPolymer(entityInst, override);
		}

		if (entityInst.getSchemClass().isa(ReactomeJavaConstants.EntitySet)) {
			return handleEntitySet(entityInst, override);
		}

		if (entityInst.getSchemClass().isa(ReactomeJavaConstants.SimpleEntity)) {
			logger.info("PE is a SimplyEntity, using original instance");
			return entityInst;
		}

		logger.warn("Unknown PhysicalEntity class: " + entityInst.getClass());
		return null;
	}

	private GKInstance handleGenomeEncodedEntity(GKInstance entityInst, boolean override) throws Exception {
		if (entityInst.getSchemClass().isa(ReactomeJavaConstants.EntityWithAccessionedSequence)) {
			return createInfEWAS(entityInst, override);
		}

		if (override) {
			logger.info("Mock GEE instance needed");
			return instanceUtilities.createMockGKInstance(entityInst);
		}

		return null;
	}

	private GKInstance handleEntitySet(GKInstance entityInst, boolean override) throws Exception {
		if (entityInst.getAttributeValue(ReactomeJavaConstants.species) != null) {
			return entitySetInferrer.createInfEntitySet(entityInst, override);
		}

		logger.info("EntitySet has no species attribute, using original instance: " + entityInst);
		return entityInst;
	}

	private boolean hasSpeciesAttribute(GKInstance instance) {
		return instance.getSchemClass().isValidAttribute(ReactomeJavaConstants.species);
	}
	
	// Function that first tries to infer any EWAS' associated with the instance. For those that have more than 1
	// returned EWAS instance, it's re-structured to a DefinedSet instance. If there is no EWAS instances inferred,
	// it will either return null or, if override is set, return a mock instance.
	private GKInstance createInfEWAS(GKInstance ewasInst, boolean override) throws Exception {
		if (homolEWASIdenticals.get(ewasInst) == null) {
			// Attempt to infer the EWAS
			List<GKInstance> infEWASInstances = getEWASInferrer().inferEWAS(ewasInst);

			if (handleInferredEWASInstances(ewasInst, infEWASInstances, override)) {
				return homolEWASIdenticals.get(ewasInst);
			}

			if (override) {
				logger.info("Mock EWAS instance needed");
				return instanceUtilities.createMockGKInstance(ewasInst);
			}
			return nullInst;
		}

		logger.info("Inferred EWAS already exists");
		return homolEWASIdenticals.get(ewasInst);
	}

	private boolean handleInferredEWASInstances(GKInstance ewasInst, List<GKInstance> infEWASInstances, boolean override)
		throws Exception {

		if (infEWASInstances.size() > 1) {
			handleMultipleEWASHomologues(ewasInst, infEWASInstances);
			return true;
		} else if (infEWASInstances.size() == 1) {
			homolEWASIdenticals.put(ewasInst, infEWASInstances.get(0));
			return true;
		}
		return false;
	}

	private void handleMultipleEWASHomologues(GKInstance ewasInst, List<GKInstance> infEWASInstances)
		throws Exception {

		logger.info("Multiple EWAS homologues produced for single EWAS. Converting to DefinedSet");
		GKInstance infDefinedSetInst = createDefinedSetForHomologues(ewasInst, infEWASInstances);
		updateDefinedSetReferences(infDefinedSetInst, ewasInst);
		homolEWASIdenticals.put(ewasInst, infDefinedSetInst);
		logger.info("Successfully converted to DefinedSet");
	}

	private GKInstance createDefinedSetForHomologues(GKInstance ewasInst, List<GKInstance> infEWASInstances)
			throws Exception {
		SchemaClass definedSetClass = dba.getSchema().getClassByName(ReactomeJavaConstants.DefinedSet);
		GKInstance infDefinedSetInst = new GKInstance(definedSetClass);
		infDefinedSetInst.setDbAdaptor(dba);

		// Set basic attributes
		infDefinedSetInst.addAttributeValue(ReactomeJavaConstants.created, getInstanceEdit());
		String definedSetName = "Homologues of " + ewasInst.getAttributeValue(ReactomeJavaConstants.name);
		infDefinedSetInst.addAttributeValue(ReactomeJavaConstants.name, definedSetName);

		// Handle compartment
		setDefinedSetCompartment(infDefinedSetInst, ewasInst);

		// Set species and members
		infDefinedSetInst.addAttributeValue(ReactomeJavaConstants.species, getSpeciesInstance());
		infDefinedSetInst.addAttributeValue(ReactomeJavaConstants.hasMember, infEWASInstances);

		// Set display name
		String definedSetDisplayName = buildDefinedSetDisplayName(infDefinedSetInst, ewasInst);
		infDefinedSetInst.setAttributeValue(ReactomeJavaConstants._displayName, definedSetDisplayName);

		return handleDefinedSetCaching(infDefinedSetInst, ewasInst);
	}

	private void setDefinedSetCompartment(GKInstance infDefinedSetInst, GKInstance ewasInst) throws Exception {
		GKInstance compartmentInstGk = (GKInstance) ewasInst.getAttributeValue(ReactomeJavaConstants.compartment);
		if (compartmentInstGk.getSchemClass().isa(ReactomeJavaConstants.Compartment)) {
			infDefinedSetInst.addAttributeValue(ReactomeJavaConstants.compartment, compartmentInstGk);
		} else {
			GKInstance newCompartmentInst = instanceUtilities.createCompartmentInstance(compartmentInstGk);
			infDefinedSetInst.addAttributeValue(ReactomeJavaConstants.compartment, newCompartmentInst);
		}
	}

	private String buildDefinedSetDisplayName(GKInstance infDefinedSetInst, GKInstance ewasInst) throws Exception {
		return infDefinedSetInst.getAttributeValue(ReactomeJavaConstants.name) +
				" [" + ((GKInstance) ewasInst.getAttributeValue(ReactomeJavaConstants.compartment)).getDisplayName() + "]";
	}

	private GKInstance handleDefinedSetCaching(GKInstance infDefinedSetInst, GKInstance ewasInst) throws Exception {
		String cacheKey = instanceUtilities.getCacheKey(
				(GKSchemaClass) infDefinedSetInst.getSchemClass(), infDefinedSetInst);
		if (definedSetIdenticals.get(cacheKey) != null) {
			return definedSetIdenticals.get(cacheKey);
		}
		infDefinedSetInst = instanceUtilities.checkForIdenticalInstances(infDefinedSetInst, ewasInst);
		definedSetIdenticals.put(cacheKey, infDefinedSetInst);
		return infDefinedSetInst;
	}

	private void updateDefinedSetReferences(GKInstance infDefinedSetInst, GKInstance ewasInst) throws Exception {
		infDefinedSetInst = instanceUtilities.addAttributeValueIfNecessary(
				infDefinedSetInst, ewasInst, ReactomeJavaConstants.inferredFrom);
		dba.updateInstanceAttribute(infDefinedSetInst, ReactomeJavaConstants.inferredFrom);

		ewasInst = instanceUtilities.addAttributeValueIfNecessary(
				ewasInst, infDefinedSetInst, ReactomeJavaConstants.inferredTo);
		dba.updateInstanceAttribute(ewasInst, ReactomeJavaConstants.inferredTo);
	}

	private EWASInferrer getEWASInferrer() {
		return this.ewasInferrer;
	}

	private GKInstance getInstanceEdit() throws Exception {
		return this.instanceUtilities.getInstanceEdit();
	}

	private ProteinCountUtility getProteinCountUtility() {
		return this.instanceUtilities.getProteinCountUtility();
	}

	private GKInstance getSpeciesInstance() throws Exception {
		return this.instanceUtilities.getSpeciesInstance();
	}
}
