package org.reactome.orthoinference;

import java.io.IOException;
import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;
import static org.reactome.util.general.CollectionUtils.safeList;

import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.GKSchemaAttribute;
import org.gk.schema.GKSchemaClass;
import org.gk.schema.SchemaClass;
import org.json.simple.parser.ParseException;
import org.reactome.release.common.database.InstanceEditUtils;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.context.annotation.Bean;
import org.springframework.stereotype.Component;

// GenerateInstance is meant to act as a catch-all for functions that are instance-oriented, such as creating, mocking,
// or identical-checking.
@Component
public class InstanceUtilities {

	private static final Logger logger = LogManager.getLogger();

	private final long DISEASE_PATHWAY_DB_ID = 1643685L;

	private MySQLAdaptor dba;
	private long personId;
	private String targetSpeciesCode;
	private SpeciesConfig speciesConfig;

	private GKInstance speciesInst;
	private GKInstance instanceEdit;
	private Map<String,GKInstance> mockedIdenticals = new HashMap<>();

	public InstanceUtilities(
		@Qualifier("currentDBA") MySQLAdaptor dba,
		@Qualifier("personId") long personId,
		@Qualifier("targetSpeciesCode") String targetSpeciesCode,
		SpeciesConfig speciesConfig
	) {
		this.dba = dba;
		this.personId = personId;
		this.targetSpeciesCode = targetSpeciesCode;
		this.speciesConfig = speciesConfig;
	}

	@Bean(name = "instanceEditInst")
	public GKInstance getInstanceEdit() throws Exception {
		if (instanceEdit == null) {
			instanceEdit = InstanceEditUtils.createInstanceEdit(
				dba, personId, "org.reactome.orthoinference");
		}
		return instanceEdit;
	}

	// Find the instance specific to this species
	@Bean(name = "speciesInst")
	public GKInstance getSpeciesInstance() throws Exception {
		if (speciesInst == null) {
			SchemaClass referenceDb = dba.getSchema().getClassByName(ReactomeJavaConstants.Species);
			speciesInst = new GKInstance(referenceDb);
			speciesInst.setDbAdaptor(dba);
			speciesInst.addAttributeValue(ReactomeJavaConstants.created, getInstanceEdit());
			speciesInst.addAttributeValue(ReactomeJavaConstants.name, getSpeciesName());
			speciesInst.addAttributeValue(ReactomeJavaConstants._displayName, getSpeciesName());
			speciesInst = checkForIdenticalInstances(speciesInst, null);
		}
		return speciesInst;
	}

	// Creates new instance that will be inferred based on the incoming instances class
	public GKInstance createNewInferredGKInstance(GKInstance instanceToBeInferred) throws Exception {
		String reactionClass = instanceToBeInferred.getSchemClass().getName();
		if (reactionClass.matches(ReactomeJavaConstants.ReferenceIsoform)) {
			reactionClass = ReactomeJavaConstants.ReferenceGeneProduct;
		}
		SchemaClass instanceClass = dba.getSchema().getClassByName(reactionClass);
		GKInstance inferredInst = new GKInstance(instanceClass);
		inferredInst.setDbAdaptor(dba);
		inferredInst.addAttributeValue(ReactomeJavaConstants.created, instanceEdit);
		if (instanceToBeInferred.getSchemClass().isValidAttribute(ReactomeJavaConstants.compartment) &&
			instanceToBeInferred.getAttributeValue(ReactomeJavaConstants.compartment) != null) {
			for (Object compartmentInst : instanceToBeInferred.getAttributeValuesList(ReactomeJavaConstants.compartment)) {
				GKInstance compartmentInstGk = (GKInstance) compartmentInst;
				if (compartmentInstGk.getSchemClass().isa(ReactomeJavaConstants.Compartment)) {
					inferredInst.addAttributeValue(ReactomeJavaConstants.compartment, compartmentInstGk);
				} else {
					GKInstance newCompartmentInst = createCompartmentInstance(compartmentInstGk);
					inferredInst.addAttributeValue(ReactomeJavaConstants.compartment, newCompartmentInst);
				}
			}
		}
		if (instanceToBeInferred.getSchemClass().isValidAttribute(ReactomeJavaConstants.species) &&
			instanceToBeInferred.getAttributeValue(ReactomeJavaConstants.species) != null) {
			inferredInst.addAttributeValue(ReactomeJavaConstants.species, speciesInst);
		}
		return inferredInst;
	}
	
	// Some 'Compartment' instances were actually 'GO_CellularComponent' instances. This meant that the instances that
	// were pulled from the original instance's Compartment attribute could not be added to the new instance, due to
	// them being a GO_CellularComponent. This function is the workaround, producing a Compartment instance that
	// contains all the same attribute values.
	@SuppressWarnings("unchecked")
	public GKInstance createCompartmentInstance(GKInstance compartmentInstGk) throws Exception {
		logger.warn(compartmentInstGk + " is a " + compartmentInstGk.getSchemClass() + " instead of a Compartment" +
			" -- creating new Compartment instance");
		SchemaClass compartmentClass = dba.getSchema().getClassByName(ReactomeJavaConstants.Compartment);
		GKInstance newCompartmentInst = new GKInstance(compartmentClass);
		newCompartmentInst.setDbAdaptor(dba);
		Collection<GKSchemaAttribute> compartmentAttributes = compartmentClass.getAttributes();
		for (GKSchemaAttribute compartmentAttribute : compartmentAttributes) {
			if (!compartmentAttribute.getName().matches("DB_ID") &&
				compartmentInstGk.getAttributeValue(compartmentAttribute.getName()) != null) {
				for (Object attribute : compartmentInstGk.getAttributeValuesList(compartmentAttribute.getName())) {
					newCompartmentInst.addAttributeValue(compartmentAttribute.getName(), attribute);
				}
			}
		}
		newCompartmentInst = checkForIdenticalInstances(newCompartmentInst, null);
		return newCompartmentInst;
	}

	// Equivalent to create_ghost from Perl; Returns a mock homologue that is needed in cases where an inference is
	// rejected, but the component isn't essential for the inference to be completed.
	public GKInstance createMockGKInstance(GKInstance instanceToBeMocked) throws Exception {
		SchemaClass genomeEncodedEntityClass = dba.getSchema().getClassByName(ReactomeJavaConstants.GenomeEncodedEntity);
		GKInstance mockedInst = new GKInstance(genomeEncodedEntityClass);
		mockedInst.setDbAdaptor(dba);
		mockedInst.addAttributeValue(ReactomeJavaConstants.created, instanceEdit);
		String mockedInstName = (String) instanceToBeMocked.getAttributeValue(ReactomeJavaConstants.name);
		mockedInst.addAttributeValue(ReactomeJavaConstants.name, "Ghost homologue of " + mockedInstName);
		mockedInst.addAttributeValue(ReactomeJavaConstants._displayName, "Ghost homologue of " +
			instanceToBeMocked.getAttributeValue(ReactomeJavaConstants._displayName));
		mockedInst.addAttributeValue(ReactomeJavaConstants.inferredFrom, instanceToBeMocked);
		mockedInst.addAttributeValue(ReactomeJavaConstants.species, speciesInst);
		mockedInst.addAttributeValue(ReactomeJavaConstants.compartment, instanceToBeMocked.getAttributeValue(ReactomeJavaConstants.compartment));

		// Caching based on an instance's defining attributes. This reduces the number of 'checkForIdenticalInstance'
		// calls, which is slow.
		String cacheKey = getCacheKey((GKSchemaClass) mockedInst.getSchemClass(), mockedInst);
		if (mockedIdenticals.get(cacheKey) != null) {
			mockedInst = mockedIdenticals.get(cacheKey);
		} else {
			mockedInst = checkForIdenticalInstances(mockedInst, instanceToBeMocked);
			mockedIdenticals.put(cacheKey, mockedInst);
		}
		instanceToBeMocked = addAttributeValueIfNecessary(instanceToBeMocked, mockedInst, ReactomeJavaConstants.inferredTo);
		dba.updateInstanceAttribute(instanceToBeMocked, ReactomeJavaConstants.inferredTo);

		return mockedInst;
	}

	// Checks that equivalent instances don't already exist in the DB, substituting if they do
	public GKInstance checkForIdenticalInstances(GKInstance inferredInst, GKInstance originalInst)
		throws Exception {
		@SuppressWarnings("unchecked")
		Collection<GKInstance> identicalInstances = dba.fetchIdenticalInstances(inferredInst);
		if (identicalInstances != null) {
			if (identicalInstances.size() == 1) {
				return identicalInstances.iterator().next();
			} else {
				// TODO: In future, could iterate through array of returned values and pull the 'most identical'.
				// For now, this mimics Perl.
				return identicalInstances.iterator().next();
			}
		} else {
			if (inferredInst.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity)) {
//				GKInstance orthoStableIdentifierInst = getUtils().getStableIdentifierGenerator()
//					.generateOrthologousStableId(inferredInst, originalInst);
//				inferredInst.addAttributeValue(ReactomeJavaConstants.stableIdentifier, orthoStableIdentifierInst);
			}
			dba.storeInstance(inferredInst);
			return inferredInst;
		}
	}
	// Checks if the instanceToCheck already contains the instanceToUse in the multi-value attribute
	@SuppressWarnings("unchecked")
	public GKInstance addAttributeValueIfNecessary(GKInstance instanceToBeCheckedForExistingAttribute,
														  GKInstance instanceContainingAttributeToBeChecked,
														  String attribute) throws Exception {
		// Original version of this function had two checks: For 'multivalue attribute' and for 'instance-type object'. 
		// Now we know the only attributes being passed through here are inferredTo, inferredFrom, orthologousEvent,
		// and hasEvent, which are all multivalue attribute classes.
		// We also know that it will always be a GKInstance passed through here (see arguments), so we are able to
		// forego both checks.
		Collection<GKInstance> attributeInstancesFromInferredInstance =
			instanceToBeCheckedForExistingAttribute.getAttributeValuesList(attribute);
		Set<Long> dbIdsFromAttributeInstances = new HashSet<>();
		for (GKInstance attributeInstance : attributeInstancesFromInferredInstance) {
			dbIdsFromAttributeInstances.add(attributeInstance.getDBID());
		}
		boolean attributeExists = false;
		for (Long attributeInstanceDbId : dbIdsFromAttributeInstances) {
			if (attributeInstanceDbId == instanceContainingAttributeToBeChecked.getDBID()) {
				attributeExists = true;
			}
		}
		if (!attributeExists) {
			instanceToBeCheckedForExistingAttribute.addAttributeValue(
				attribute, instanceContainingAttributeToBeChecked);
		}
		return instanceToBeCheckedForExistingAttribute;
	}
	
	// Caching. This function goes through each defining attribute of the incoming instance and produces a string of
	// the attribute values (DB IDs if the attribute is an instance).
	// This allows for identical instances held in memory to be used before trying to use fetchIdenticalInstances,
	// which is expensive.
	@SuppressWarnings("unchecked")
	public String getCacheKey(GKSchemaClass instanceClass, GKInstance inferredInst) throws Exception {
		String key = "";
		for (GKSchemaAttribute definingAttr : (Collection<GKSchemaAttribute>) instanceClass.getDefiningAttributes()) {
			if (definingAttr.isMultiple()) {
				Collection<Object> multiValueAttributes = inferredInst.getAttributeValuesList(definingAttr.getName());
				if (multiValueAttributes.size() > 0) {
					for (Object attribute : multiValueAttributes) {
						if (attribute.getClass().getSimpleName().equals("GKInstance")) {
							GKInstance gkInstance = (GKInstance) attribute;
							key += gkInstance.getDBID().toString();
						} else {
							key += (String) attribute;
						}
					}
				} else {
					key += "null";
				}
			} else {
				if (definingAttr.isInstanceTypeAttribute() &&
					inferredInst.getAttributeValue(definingAttr.getName()) != null) {
					key += ((GKInstance) inferredInst.getAttributeValue(definingAttr.getName())).getDBID().toString();
				} else if (inferredInst.getAttributeValue(definingAttr.getName()) != null) {
					key += inferredInst.getAttributeValue(definingAttr.getName());
				} else {
					key += "null";
				}
			}
		}
		return key;
	}

	/**
	 * This method returns true if the only parent Pathway of the incoming Event instance is Disease.
	 * Instances with only Disease as a parent will not be inferred. Those that are a member of Disease AND
	 * another Pathway will be inferred, but the Disease Pathway (and sub-Pathway) inference will be skipped.
	 * @param eventInst -- GKInstance that will be checked for parent Pathways.
	 * @return boolean -- true if only parent is Disease Pathway, false if not.
	 * @throws Exception -- Thrown by MySQLAdaptor.
	 */
	public boolean onlyInDiseasePathway(GKInstance eventInst) throws Exception {
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
	private Set<Long> getTopLevelPathwayDbIds(GKInstance pathway) throws Exception {
		List<GKInstance> parentPathways = safeList(pathway.getReferers(ReactomeJavaConstants.Event));
		if (parentPathways.isEmpty()) {
			return new HashSet<>(Arrays.asList(pathway.getDBID()));
		}

		Set<Long> topLevelPathwayDbIds = new HashSet<>();
		for (GKInstance parentPathway : parentPathways) {
			topLevelPathwayDbIds.addAll(getTopLevelPathwayDbIds(parentPathway));
		}

		return topLevelPathwayDbIds;
	}
	
	public void setAdaptor(MySQLAdaptor dbAdaptor) {
		dba = dbAdaptor;
	}

	public void setSpeciesInstance(GKInstance speciesInstCopy) {
		speciesInst = speciesInstCopy;
	}

	public void setInstanceEdit(GKInstance instanceEditCopy) {
		instanceEdit = instanceEditCopy;
	}

	public long getDiseasePathwayDbId() {
		return DISEASE_PATHWAY_DB_ID;
	}

	private String getSpeciesName() throws IOException, ParseException {
		return this.speciesConfig.getSpeciesName(this.targetSpeciesCode);
	}
}
