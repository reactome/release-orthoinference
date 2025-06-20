package org.reactome.orthoinference;

import java.sql.SQLException;
import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;

import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.GKSchemaClass;
import org.gk.schema.SchemaClass;
import org.springframework.stereotype.Component;

@Component
public class EWASInferrer {

	private static final Logger logger = LogManager.getLogger();

	private ConfigProperties configProperties;
	private OrthologyReferenceDatabase orthologyReferenceDatabase;
	private Mappings mappings;
	private InstanceUtilities instanceUtilities;

	private static Map<String, GKInstance> referenceGeneProductIdenticals = new HashMap<>();
	private static Map<String,GKInstance> ewasIdenticals = new HashMap<>();
	private static Map<String,GKInstance> residueIdenticals = new HashMap<>();

	public EWASInferrer(
		ConfigProperties configProperties,
		OrthologyReferenceDatabase orthologyReferenceDatabase,
		Mappings mappings,
		InstanceUtilities instanceUtilities
	) {
		this.configProperties = configProperties;
		this.orthologyReferenceDatabase = orthologyReferenceDatabase;
		this.mappings = mappings;
		this.instanceUtilities = instanceUtilities;
	}

	// Creates an array of inferred EWAS instances from the homologue mappings file (hsap_species_mapping.txt).
	@SuppressWarnings("unchecked")
	public List<GKInstance> inferEWAS(GKInstance ewasInst) throws Exception {
		String rgpIdentifier = getReferenceEntityIdentifier(ewasInst);
		String[] homologues = getMappings().getHomologues(rgpIdentifier);
		if (homologues == null) {
			logger.info("Could not infer EWAS, unable to find homologue for " + rgpIdentifier);
			return new ArrayList<>();
		}
		List<GKInstance> infEWASInstances = new ArrayList<>();

		// Iterate through the array of homologue mappings, attempting to infer EWAS instances for each.
		logger.info("EWAS homologue(s): " + Arrays.toString(homologues));
		for (String homologue : homologues) {
			if (!hasEnsEMBLGeneMapping(getHomologueId(homologue))) {
				logger.info("Gene ID corresponding to " + homologue + " not found in gene_protein_mapping file" +
					" -- skipping EWAS inference");
				continue;
			}

			infEWASInstances.add(inferEWASHomologue(homologue, ewasInst));
		}
		logger.info("Total orthologous EWAS' created: " + infEWASInstances.size());
		return infEWASInstances;
	}

	private GKInstance inferEWASHomologue(String homologue, GKInstance ewasInst) throws Exception {
		logger.info("Homologue:" + homologue + "  Source:" + getReferenceEntityIdentifier(ewasInst));

		String homologueId = getHomologueId(homologue);
		if (referenceGeneProductIdenticals.get(homologueId) == null) {
			logger.info("Creating ReferenceGeneProduct for " + homologue);

			referenceGeneProductIdenticals.put(homologueId, inferReferenceGeneProduct(ewasInst, homologue));
		}
		GKInstance infReferenceGeneProductInst = referenceGeneProductIdenticals.get(homologueId);

		// Creating inferred EWAS
		GKInstance infEWASInst = instanceUtilities.createNewInferredGKInstance(ewasInst);
		infEWASInst.addAttributeValue(ReactomeJavaConstants.referenceEntity, infReferenceGeneProductInst);

		inferStartAndEndCoordinates(ewasInst, infEWASInst);
		infEWASInst.addAttributeValue(ReactomeJavaConstants.name, homologueId);

		// If the species-specific gene name was retrieved from UniProt's mapping service, it is added
		// as the primary name for the EWAS.
		if (getMappings().homologueHasGeneName(homologueId)) {
			List<String> ewasNames = infEWASInst.getAttributeValuesList(ReactomeJavaConstants.name);
			infEWASInst.setAttributeValue(ReactomeJavaConstants.name, getMappings().getGeneName(homologueId));
			for (String ewasName : ewasNames) {
				infEWASInst.addAttributeValue(ReactomeJavaConstants.name, ewasName);
			}

		}
		// New display name is generated using the updated 'name' attribute
		infEWASInst.setAttributeValue(
			ReactomeJavaConstants._displayName, InstanceDisplayNameGenerator.generateDisplayName(infEWASInst));

		// Infer residue modifications. This was another step where the name of an EWAS can change.
		// For this, it is based on the existence of the string 'phospho' in the name of the psiMod
		// attribute.
		// If true, 'phospho-' is prepended to the EWAS' name attribute.
		List<GKInstance> infModifiedResidueInstances = new ArrayList<>();
		boolean phosFlag = true;
		for (GKInstance modifiedResidueInst : getModifiedResidues(ewasInst)) {
			logger.info("Inferring ModifiedResidue: " + modifiedResidueInst);
			String infModifiedResidueDisplayName = "";
			GKInstance infModifiedResidueInst = instanceUtilities.createNewInferredGKInstance(modifiedResidueInst);
			infModifiedResidueInst.addAttributeValue(ReactomeJavaConstants.referenceSequence, infReferenceGeneProductInst);
			infModifiedResidueDisplayName += infReferenceGeneProductInst.getDisplayName();
			for (int coordinateValue :
				(Collection<Integer>) modifiedResidueInst.getAttributeValuesList(ReactomeJavaConstants.coordinate)) {
				infModifiedResidueInst.addAttributeValue(ReactomeJavaConstants.coordinate, coordinateValue);
			}
			if (infModifiedResidueInst.getSchemClass().isValidAttribute(ReactomeJavaConstants.modification)) {
				for (GKInstance modifiedInst :
					(Collection<GKInstance>) modifiedResidueInst.getAttributeValuesList(ReactomeJavaConstants.modification)) {
					infModifiedResidueInst.addAttributeValue(ReactomeJavaConstants.modification, modifiedInst);
				}
				if (infModifiedResidueInst.getAttributeValue(ReactomeJavaConstants.modification) != null) {
					infModifiedResidueDisplayName += " " +
						((GKInstance) infModifiedResidueInst.getAttributeValue(ReactomeJavaConstants.modification)).getDisplayName();
				}
			}
			// Update name depending on the presence of 'phospho' in the Psimod's name attribute
			GKInstance firstPsiModInst = (GKInstance) modifiedResidueInst.getAttributeValue(ReactomeJavaConstants.psiMod);
			if (phosFlag && firstPsiModInst.getAttributeValue(ReactomeJavaConstants.name).toString().contains("phospho")) {
				String phosphoName = "phospho-" + infEWASInst.getAttributeValue(ReactomeJavaConstants.name);
				List<String> ewasNames = (ArrayList<String>) infEWASInst.getAttributeValuesList(ReactomeJavaConstants.name);
				String originalName = ewasNames.remove(0);
				infEWASInst.setAttributeValue(ReactomeJavaConstants.name, phosphoName);
				// In the Perl version, this code block modifies the 'name' attribute to include
				// 'phosopho-', but in the process it drops the other names contained. I believe this is
				// unintentional.
				// This would mean attributes without the 'phospho- ' addition would retain their array of
				// names, while attributes containing 'phospho-' would only contain a single name attribute.
				// I've assumed this is incorrect for the rewrite -- Instances that modify the name
				// attribute to prepend 'phospho-' retain their name array. (Justin Cook 2018)
				infEWASInst.addAttributeValue(ReactomeJavaConstants.name, ewasNames);
				String phosphoDisplayName = phosphoName +
					" [" + ((GKInstance) ewasInst.getAttributeValue(ReactomeJavaConstants.compartment)).getDisplayName() + "]";
				infEWASInst.setAttributeValue(ReactomeJavaConstants._displayName, phosphoDisplayName);
				// This flag ensures the 'phospho-' is only prepended once.
				logger.info("Updated EWAS name to reflect phosphorylation. Original: " + originalName +
					". Updated: " + phosphoName);
				phosFlag = false;
			}
			for (GKInstance psiModInst :
				(Collection<GKInstance>) modifiedResidueInst.getAttributeValuesList(ReactomeJavaConstants.psiMod)) {
				infModifiedResidueInst.addAttributeValue(ReactomeJavaConstants.psiMod, psiModInst);
			}
			if (infModifiedResidueInst.getAttributeValue(ReactomeJavaConstants.psiMod) != null) {
				infModifiedResidueDisplayName += " " +
					((GKInstance) infModifiedResidueInst.getAttributeValue(ReactomeJavaConstants.psiMod)).getDisplayName();
			}
			infModifiedResidueInst.setAttributeValue(
				ReactomeJavaConstants._displayName, modifiedResidueInst.getAttributeValue(ReactomeJavaConstants._displayName));
			// Update name to reflect that coordinate values are taken from humans. This takes place after
			// cache retrieval, since the name from DB won't contain updated name.
			if (modifiedResidueInst.getAttributeValue(ReactomeJavaConstants.coordinate) != null) {
				String newModifiedResidueDisplayName =
					modifiedResidueInst.getAttributeValue(ReactomeJavaConstants._displayName).toString() + " (in Homo sapiens)";
				infModifiedResidueInst.setAttributeValue(ReactomeJavaConstants._displayName, newModifiedResidueDisplayName);

			} else {
				if (infModifiedResidueInst.getSchemClass().isa(ReactomeJavaConstants.InterChainCrosslinkedResidue)) {
					infModifiedResidueInst.setDisplayName(infModifiedResidueDisplayName);
				}
			}
			// Database-checker gave errors related to missing 'secondReferenceSequence' and
			// 'equivalentTo' attributes in InterChainCrosslinkedResidues
			// This was because they were never populated. This block is the fix.
			if (infModifiedResidueInst.getSchemClass().isa(ReactomeJavaConstants.InterChainCrosslinkedResidue)) {
				if (modifiedResidueInst.getAttributeValue(ReactomeJavaConstants.secondReferenceSequence) != null) {
					for (GKInstance secondRefSequenceInst : (Collection<GKInstance>)
						modifiedResidueInst.getAttributeValuesList(ReactomeJavaConstants.secondReferenceSequence)) {

						infModifiedResidueInst.addAttributeValue(
							ReactomeJavaConstants.secondReferenceSequence, secondRefSequenceInst);
					}
				}
				if (modifiedResidueInst.getAttributeValue("equivalentTo") != null) {
					for (GKInstance equivalentToInst : (Collection<GKInstance>)
						modifiedResidueInst.getAttributeValuesList("equivalentTo")) {

						infModifiedResidueInst.addAttributeValue("equivalentTo", equivalentToInst);
					}
				}
			}
			// Caching based on an instance's defining attributes. This reduces the number of
			// 'checkForIdenticalInstance' calls, which slows things.
			String cacheKey = instanceUtilities.getCacheKey(
				(GKSchemaClass) infModifiedResidueInst.getSchemClass(), infModifiedResidueInst);
			if (residueIdenticals.get(cacheKey) != null) {
				infModifiedResidueInst = residueIdenticals.get(cacheKey);
			} else {
				infModifiedResidueInst = instanceUtilities.checkForIdenticalInstances(
					infModifiedResidueInst, null);
				residueIdenticals.put(cacheKey, infModifiedResidueInst);
			}
			logger.info("Successfully inferred ModifiedResidue");
			infModifiedResidueInstances.add(infModifiedResidueInst);
		}
		infEWASInst.addAttributeValue(ReactomeJavaConstants.hasModifiedResidue, infModifiedResidueInstances);
		// Caching based on an instance's defining attributes. This reduces the number of
		// 'checkForIdenticalInstance' calls, which slows things.
		String cacheKey = instanceUtilities.getCacheKey(
			(GKSchemaClass) infEWASInst.getSchemClass(), infEWASInst);
		if (ewasIdenticals.get(cacheKey) != null) {
			infEWASInst = ewasIdenticals.get(cacheKey);
		} else {
			infEWASInst = instanceUtilities.checkForIdenticalInstances(infEWASInst, ewasInst);
			ewasIdenticals.put(cacheKey, infEWASInst);
		}

		infEWASInst = instanceUtilities.addAttributeValueIfNecessary(infEWASInst, ewasInst, ReactomeJavaConstants.inferredFrom);
		getCurrentDBA().updateInstanceAttribute(infEWASInst, ReactomeJavaConstants.inferredFrom);
		ewasInst = instanceUtilities.addAttributeValueIfNecessary(ewasInst, infEWASInst, ReactomeJavaConstants.inferredTo);
		getCurrentDBA().updateInstanceAttribute(ewasInst, ReactomeJavaConstants.inferredTo);
		logger.info("Successfully inferred EWAS instance for " + homologue + " homologue");
		return infEWASInst;
	}

	private List<GKInstance> getModifiedResidues(GKInstance ewasInst) throws Exception {
		return (List<GKInstance>) ewasInst.getAttributeValuesList(ReactomeJavaConstants.hasModifiedResidue);
	}

//	private static GKInstance inferModifiedResidue(GKInstance modifiedResidueInst, GKInstance ewasInst, GKInstance infEWASInst, GKInstance infReferenceGeneProductInst, boolean phosFlag) throws Exception {
//
//	}

	private String getModifiedResidueDisplayName() {
		return "";
	}

	private String getReferenceEntityIdentifier(GKInstance ewasInst) throws Exception {
		return ((GKInstance) ewasInst
			.getAttributeValue(ReactomeJavaConstants.referenceEntity))
			.getAttributeValue(ReactomeJavaConstants.identifier).toString();
	}

	private GKInstance inferReferenceGeneProduct(GKInstance ewasInst, String homologue) throws Exception {
		String homologueId = getHomologueId(homologue);

		GKInstance infReferenceGeneProductInst = instanceUtilities.createNewInferredGKInstance(
			(GKInstance) ewasInst.getAttributeValue(ReactomeJavaConstants.referenceEntity));
		infReferenceGeneProductInst.addAttributeValue(ReactomeJavaConstants.identifier, homologueId);
		// Reference DB can differ between homologue mappings, but can be differentiated by the
		// 'homologueSource' found in each mapping.
		// With PANTHER data, the Protein IDs are exclusively UniProt
		GKInstance referenceDatabaseInst = getReferenceDatabase().getReferenceDatabase(getHomologueSource(homologue));
		infReferenceGeneProductInst.addAttributeValue(ReactomeJavaConstants.referenceDatabase, referenceDatabaseInst);

		// Creates ReferenceDNASequence instance from ReferenceEntity
		List<GKInstance> inferredReferenceDNAInstances = createReferenceDNASequences(homologueId);
		infReferenceGeneProductInst.addAttributeValue(ReactomeJavaConstants.referenceGene, inferredReferenceDNAInstances);

		infReferenceGeneProductInst.addAttributeValue(ReactomeJavaConstants.species, getSpeciesInstance());
		String referenceGeneProductSource = getReferenceDatabase().getReferenceDatabaseSourceName(getHomologueSource(homologue));
		infReferenceGeneProductInst.setAttributeValue(
			ReactomeJavaConstants._displayName, referenceGeneProductSource + ":" + homologueId);

		// GeneName value comes from UniProt's identifier mapping service.
		if (getMappings().homologueHasGeneName(homologueId)) {
			infReferenceGeneProductInst.addAttributeValue(ReactomeJavaConstants.geneName, getMappings().getGeneName(homologueId));
		}

		logger.info("ReferenceGeneProduct instance created");
		infReferenceGeneProductInst = instanceUtilities.checkForIdenticalInstances(
			infReferenceGeneProductInst, null);
		return infReferenceGeneProductInst;
	}

	private String getHomologueSource(String homologue) {
		return homologue.contains(":") ? homologue.split(":")[0] : "";
	}

	private String getHomologueId(String homologue) {
		return homologue.contains(":") ? homologue.split(":")[1] : homologue;
	}

	private void inferStartAndEndCoordinates(GKInstance ewasInst, GKInstance infEWASInst) throws Exception {
		// Method for adding start/end coordinates. It is convoluted due to a quirk with assigning the
		// name differently based on coordinate value (see infer_events.pl lines 1190-1192).
		// The name of the entity needs to be at the front of the 'name' array if the coordinate is over 1,
		// and rearranging arrays in Java for this was a bit tricky.
		for (int startCoord : (Collection<Integer>) ewasInst.getAttributeValuesList(ReactomeJavaConstants.startCoordinate)) {
			infEWASInst.addAttributeValue(ReactomeJavaConstants.startCoordinate, startCoord);
		}
		for (int endCoord : (Collection<Integer>) ewasInst.getAttributeValuesList(ReactomeJavaConstants.endCoordinate)) {
			infEWASInst.addAttributeValue(ReactomeJavaConstants.endCoordinate, endCoord);
		}
		if (infEWASInst.getAttributeValue(ReactomeJavaConstants.startCoordinate) != null &&
			(int) infEWASInst.getAttributeValue(ReactomeJavaConstants.startCoordinate) > 1 ||
			infEWASInst.getAttributeValue(ReactomeJavaConstants.endCoordinate) != null &&
				(int) infEWASInst.getAttributeValue(ReactomeJavaConstants.endCoordinate) > 1) {
			List<String> infEWASInstNames = (ArrayList<String>) (ewasInst).getAttributeValuesList(ReactomeJavaConstants.name);
			infEWASInst.addAttributeValue(ReactomeJavaConstants.name, infEWASInstNames.get(0));
		}
	}

	// Homologous Protein IDs can exist in ${source}_${target}_mapping.txt but the corresponding Gene ID might not
	// exist in ${target}_gene_protein_mapping.txt.
	// This is different from when we built Orthopairs files using Compara, since the homology mapping file was
	// generated using IDs from the gene-protein file.
	// This function prevents a Null Exception from killing the entire Reaction's inference, rather than just the
	// EWAS inference.
	private boolean hasEnsEMBLGeneMapping(String homologueId) {
		return getMappings().homologueHasEnsEMBLGeneMapping(homologueId);
	}

	// Creates ReferenceGeneSequence instance based on ENSG identifier mapped to protein.
	// Creates an instance for the primary database and an alternate, if it exists.
	private List<GKInstance> createReferenceDNASequences(String homologueId) throws Exception {
		List<GKInstance> referenceDNAInstances = new ArrayList<>();
		List<String> ensgIds = getMappings().getEnsEMBLGeneIdentifiers(homologueId);
		logger.info("Gene ID(s): " + ensgIds);
		for (String ensgId : ensgIds) {
			referenceDNAInstances.add(createReferenceDNASequence(ensgId));
			if (getReferenceDatabase().alternateReferenceDatabaseExists()) {
				referenceDNAInstances.add(createReferenceDNASequenceWithAlternateReferenceDatabase(ensgId));
			}
		}
		logger.info("Total ReferenceDNASequence instance(s) created: " + referenceDNAInstances.size());
		return referenceDNAInstances;
	}

	private GKInstance createReferenceDNASequence(String ensgId) throws Exception {
		logger.info("Creating ReferenceDNASequence for " + ensgId);
		SchemaClass referenceDNAClass = getCurrentDBA().getSchema().getClassByName(ReactomeJavaConstants.ReferenceDNASequence);
		GKInstance referenceDNAInst = new GKInstance(referenceDNAClass);
		referenceDNAInst.setDbAdaptor(getCurrentDBA());
		referenceDNAInst.addAttributeValue(ReactomeJavaConstants.created, getInstanceEdit());
		referenceDNAInst.addAttributeValue(ReactomeJavaConstants.identifier, ensgId);
		referenceDNAInst.addAttributeValue(ReactomeJavaConstants.referenceDatabase, getReferenceDatabase().getEnsemblDBInst());
		referenceDNAInst.addAttributeValue(ReactomeJavaConstants.species, getSpeciesInstance());
		referenceDNAInst.setAttributeValue(ReactomeJavaConstants._displayName, "ENSEMBL:" + ensgId);
		referenceDNAInst = instanceUtilities.checkForIdenticalInstances(referenceDNAInst, null);
		return referenceDNAInst;
	}

	private GKInstance createReferenceDNASequenceWithAlternateReferenceDatabase(String ensgId) throws Exception {
		logger.info("Creating ReferenceDNASequence for " + ensgId + " using alternate reference database");
		SchemaClass referenceDNAClass = getCurrentDBA().getSchema().getClassByName(ReactomeJavaConstants.ReferenceDNASequence);
		GKInstance alternateRefDNAInst = new GKInstance(referenceDNAClass);
		alternateRefDNAInst.setDbAdaptor(getCurrentDBA());
		String altDbIdentifier = ensgId;
		if (getReferenceDatabase().getAlternateReferenceDbId() != null) {
			altDbIdentifier = altDbIdentifier.replaceAll(getReferenceDatabase().getAlternateReferenceDbId(), "");
		}
		alternateRefDNAInst.addAttributeValue(ReactomeJavaConstants.created, getInstanceEdit());
		alternateRefDNAInst.addAttributeValue(ReactomeJavaConstants.identifier, altDbIdentifier);
		alternateRefDNAInst.addAttributeValue(ReactomeJavaConstants.referenceDatabase, getReferenceDatabase().getAlternateReferenceDatabase());
		alternateRefDNAInst.addAttributeValue(ReactomeJavaConstants.species, getSpeciesInstance());
		alternateRefDNAInst.setAttributeValue(
			ReactomeJavaConstants._displayName, getReferenceDatabase().getAlternateReferenceDatabase().getAttributeValue(ReactomeJavaConstants.name) + ":" + ensgId);
		alternateRefDNAInst = instanceUtilities.checkForIdenticalInstances(alternateRefDNAInst, null);
		return alternateRefDNAInst;
	}

	private Mappings getMappings() {
		return this.mappings;
	}

	private OrthologyReferenceDatabase getReferenceDatabase() {
		return this.orthologyReferenceDatabase;
	}

	private MySQLAdaptor getCurrentDBA() throws SQLException {
		return this.configProperties.getCurrentDBA();
	}

	private GKInstance getSpeciesInstance() throws Exception {
		return this.instanceUtilities.getSpeciesInstance();
	}

	private GKInstance getInstanceEdit() throws Exception {
		return this.instanceUtilities.getInstanceEdit();
	}
}
