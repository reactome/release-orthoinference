package org.reactome.orthoinference;

import java.sql.SQLException;
import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;

import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.GKSchemaClass;
import org.gk.schema.SchemaClass;
import org.reactome.curation.model.SimpleInstance;
import org.springframework.stereotype.Component;

@Component
public class EWASInferrer {

	private static final Logger logger = LogManager.getLogger();

	private ConfigProperties configProperties;
	private OrthologyReferenceDatabase orthologyReferenceDatabase;
	private Mappings mappings;
	private InstanceUtilities instanceUtilities;

	private static Map<String, SimpleInstance> referenceGeneProductIdenticals = new HashMap<>();
	private static Map<String, SimpleInstance> ewasIdenticals = new HashMap<>();
	private static Map<String, SimpleInstance> residueIdenticals = new HashMap<>();

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
	public List<SimpleInstance> inferEWAS(SimpleInstance ewasInst) throws Exception {
		String rgpIdentifier = getReferenceEntityIdentifier(ewasInst);
		String[] homologues = getMappings().getHomologues(rgpIdentifier);
		if (homologues == null) {
			logger.info("Could not infer EWAS, unable to find homologue for " + rgpIdentifier);
			return new ArrayList<>();
		}
		List<SimpleInstance> infEWASInstances = new ArrayList<>();

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

	private SimpleInstance inferEWASHomologue(String homologue, SimpleInstance ewasInst) throws Exception {
		logger.info("Homologue:" + homologue + "  Source:" + getReferenceEntityIdentifier(ewasInst));

		String homologueId = getHomologueId(homologue);
		if (referenceGeneProductIdenticals.get(homologueId) == null) {
			logger.info("Creating ReferenceGeneProduct for " + homologue);

			referenceGeneProductIdenticals.put(homologueId, inferReferenceGeneProduct(ewasInst, homologue));
		}
		SimpleInstance infReferenceGeneProductInst = referenceGeneProductIdenticals.get(homologueId);

		// Creating inferred EWAS
		SimpleInstance infEWASInst = instanceUtilities.createNewInferredGKInstance(ewasInst);
		infEWASInst.setAttribute(ReactomeJavaConstants.referenceEntity, infReferenceGeneProductInst);

		inferStartAndEndCoordinates(ewasInst, infEWASInst);
		infEWASInst.setAttribute(ReactomeJavaConstants.name, homologueId);

		// If the species-specific gene name was retrieved from UniProt's mapping service, it is added
		// as the primary name for the EWAS.
		if (getMappings().homologueHasGeneName(homologueId)) {
			List<String> ewasNames = (List<String>) infEWASInst.getAttribute(ReactomeJavaConstants.name);

			ewasNames.add(0, getMappings().getGeneName(homologueId));
			infEWASInst.setAttribute(ReactomeJavaConstants.name, ewasNames);
		}
		// New display name is generated using the updated 'name' attribute
		infEWASInst.setDisplayName(generateEWASDisplayName(infEWASInst));

		// Infer residue modifications. This was another step where the name of an EWAS can change.
		// For this, it is based on the existence of the string 'phospho' in the name of the psiMod
		// attribute.
		// If true, 'phospho-' is prepended to the EWAS' name attribute.
		List<SimpleInstance> infModifiedResidueInstances = new ArrayList<>();
		boolean phosFlag = true;
		for (SimpleInstance modifiedResidueInst : getModifiedResidues(ewasInst)) {
			logger.info("Inferring ModifiedResidue: " + modifiedResidueInst);
			String infModifiedResidueDisplayName = "";
			SimpleInstance infModifiedResidueInst = instanceUtilities.createNewInferredGKInstance(modifiedResidueInst);
			infModifiedResidueInst.setAttribute(ReactomeJavaConstants.referenceSequence, infReferenceGeneProductInst);
			infModifiedResidueDisplayName += infReferenceGeneProductInst.getDisplayName();
			for (int coordinateValue :
				(Collection<Integer>) modifiedResidueInst.getAttribute(ReactomeJavaConstants.coordinate)) {
				infModifiedResidueInst.setAttribute(ReactomeJavaConstants.coordinate, coordinateValue);
			}
			if (instanceUtilities.isValidAttribute(infModifiedResidueInst, ReactomeJavaConstants.modification)) {
				for (SimpleInstance modifiedInst :
					(Collection<SimpleInstance>) modifiedResidueInst.getAttribute(ReactomeJavaConstants.modification)) {
					infModifiedResidueInst.setAttribute(ReactomeJavaConstants.modification, modifiedInst);
				}
				if (infModifiedResidueInst.getAttribute(ReactomeJavaConstants.modification) != null) {
					infModifiedResidueDisplayName += " " +
						((SimpleInstance) infModifiedResidueInst.getAttribute(ReactomeJavaConstants.modification)).getDisplayName();
				}
			}
			// Update name depending on the presence of 'phospho' in the Psimod's name attribute
			SimpleInstance firstPsiModInst = (SimpleInstance) modifiedResidueInst.getAttribute(ReactomeJavaConstants.psiMod);
			if (phosFlag && firstPsiModInst.getAttribute(ReactomeJavaConstants.name).toString().contains("phospho")) {
				String phosphoName = "phospho-" + infEWASInst.getAttribute(ReactomeJavaConstants.name);
				List<String> ewasNames = (ArrayList<String>) infEWASInst.getAttribute(ReactomeJavaConstants.name);
				String originalName = ewasNames.remove(0);

				ewasNames.add(0, phosphoName);
				// In the Perl version, this code block modifies the 'name' attribute to include
				// 'phosopho-', but in the process it drops the other names contained. I believe this is
				// unintentional.
				// This would mean attributes without the 'phospho- ' addition would retain their array of
				// names, while attributes containing 'phospho-' would only contain a single name attribute.
				// I've assumed this is incorrect for the rewrite -- Instances that modify the name
				// attribute to prepend 'phospho-' retain their name array. (Justin Cook 2018)
				infEWASInst.setAttribute(ReactomeJavaConstants.name, ewasNames);

				String phosphoDisplayName = phosphoName +
					" [" + ((GKInstance) ewasInst.getAttribute(ReactomeJavaConstants.compartment)).getDisplayName() + "]";
				infEWASInst.setDisplayName(phosphoDisplayName);
				// This flag ensures the 'phospho-' is only prepended once.
				logger.info("Updated EWAS name to reflect phosphorylation. Original: " + originalName +
					". Updated: " + phosphoName);
				phosFlag = false;
			}
			for (GKInstance psiModInst :
				(Collection<GKInstance>) modifiedResidueInst.getAttribute(ReactomeJavaConstants.psiMod)) {
				infModifiedResidueInst.setAttribute(ReactomeJavaConstants.psiMod, psiModInst);
			}
			if (infModifiedResidueInst.getAttribute(ReactomeJavaConstants.psiMod) != null) {
				infModifiedResidueDisplayName += " " +
					((GKInstance) infModifiedResidueInst.getAttribute(ReactomeJavaConstants.psiMod)).getDisplayName();
			}
			infModifiedResidueInst.setDisplayName(modifiedResidueInst.getDisplayName());
			// Update name to reflect that coordinate values are taken from humans. This takes place after
			// cache retrieval, since the name from DB won't contain updated name.
			if (modifiedResidueInst.getAttribute(ReactomeJavaConstants.coordinate) != null) {
				String newModifiedResidueDisplayName = modifiedResidueInst.getDisplayName() + " (in Homo sapiens)";
				infModifiedResidueInst.setDisplayName(newModifiedResidueDisplayName);

			} else {
				if (instanceUtilities.isSchemaClassType(infModifiedResidueInst, ReactomeJavaConstants.InterChainCrosslinkedResidue)) {
					infModifiedResidueInst.setDisplayName(infModifiedResidueDisplayName);
				}
			}
			// Database-checker gave errors related to missing 'secondReferenceSequence' and
			// 'equivalentTo' attributes in InterChainCrosslinkedResidues
			// This was because they were never populated. This block is the fix.
			if (instanceUtilities.isSchemaClassType(infModifiedResidueInst, ReactomeJavaConstants.InterChainCrosslinkedResidue)) {
				if (modifiedResidueInst.getAttribute(ReactomeJavaConstants.secondReferenceSequence) != null) {
					List<SimpleInstance> secondReferenceSequenceInstances = new ArrayList<>();
					for (SimpleInstance secondRefSequenceInst : (Collection<SimpleInstance>)
						modifiedResidueInst.getAttribute(ReactomeJavaConstants.secondReferenceSequence)) {

						secondReferenceSequenceInstances.add(secondRefSequenceInst);
					}
					infModifiedResidueInst.setAttribute(
						ReactomeJavaConstants.secondReferenceSequence, secondReferenceSequenceInstances);
				}
				if (modifiedResidueInst.getAttribute("equivalentTo") != null) {
					List<SimpleInstance> equivalentToInstances = new ArrayList<>();
					for (SimpleInstance equivalentToInst : (Collection<SimpleInstance>)
						modifiedResidueInst.getAttribute("equivalentTo")) {
						equivalentToInstances.add(equivalentToInst);
					}
					infModifiedResidueInst.setAttribute("equivalentTo", equivalentToInstances);

				}
			}
			// Caching based on an instance's defining attributes. This reduces the number of
			// 'checkForIdenticalInstance' calls, which slows things.
			String cacheKey = instanceUtilities.getCacheKey(infModifiedResidueInst);
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
		infEWASInst.setAttribute(ReactomeJavaConstants.hasModifiedResidue, infModifiedResidueInstances);
		// Caching based on an instance's defining attributes. This reduces the number of
		// 'checkForIdenticalInstance' calls, which slows things.
		String cacheKey = instanceUtilities.getCacheKey(infEWASInst);
		if (ewasIdenticals.get(cacheKey) != null) {
			infEWASInst = ewasIdenticals.get(cacheKey);
		} else {
			infEWASInst = instanceUtilities.checkForIdenticalInstances(infEWASInst, ewasInst);
			ewasIdenticals.put(cacheKey, infEWASInst);
		}

		instanceUtilities.addInferredFromValue(infEWASInst, ewasInst);
		instanceUtilities.addInferredToValue(ewasInst, infEWASInst);
		logger.info("Successfully inferred EWAS instance for " + homologue + " homologue");
		return infEWASInst;
	}

	private List<SimpleInstance> getModifiedResidues(SimpleInstance ewasInst) {
		return (List<SimpleInstance>) ewasInst.getAttribute(ReactomeJavaConstants.hasModifiedResidue);
	}

//	private static GKInstance inferModifiedResidue(GKInstance modifiedResidueInst, GKInstance ewasInst, GKInstance infEWASInst, GKInstance infReferenceGeneProductInst, boolean phosFlag) throws Exception {
//
//	}

//	private String getModifiedResidueDisplayName() {
//		return "";
//	}

	// Mirrors InstanceDisplayNameGenerator's PhysicalEntity display name -- "<first name> [<first compartment>]", or
	// the first name alone when there is no compartment -- which can't be used as it only takes a GKInstance.
	// A SimpleInstance holds a single value as itself rather than in a List, so both forms are handled.
	private String generateEWASDisplayName(SimpleInstance ewasInst) {
		Object nameValue = ewasInst.getAttribute(ReactomeJavaConstants.name);
		if (nameValue instanceof List) {
			List<?> names = (List<?>) nameValue;
			nameValue = names.isEmpty() ? null : names.get(0);
		}
		String displayName = nameValue != null ? nameValue.toString() : "";

		Object compartmentValue = ewasInst.getAttribute(ReactomeJavaConstants.compartment);
		if (compartmentValue instanceof List) {
			List<?> compartments = (List<?>) compartmentValue;
			compartmentValue = compartments.isEmpty() ? null : compartments.get(0);
		}
		if (compartmentValue != null) {
			displayName += " [" + ((SimpleInstance) compartmentValue).getDisplayName() + "]";
		}
		return displayName;
	}

	private String getReferenceEntityIdentifier(SimpleInstance ewasInst) {
		return ((SimpleInstance) ewasInst
			.getAttribute(ReactomeJavaConstants.referenceEntity))
			.getAttribute(ReactomeJavaConstants.identifier).toString();
	}

	private SimpleInstance inferReferenceGeneProduct(SimpleInstance ewasInst, String homologue) throws Exception {
		String homologueId = getHomologueId(homologue);

		SimpleInstance infReferenceGeneProductInst = instanceUtilities.createNewInferredGKInstance(
			(SimpleInstance) ewasInst.getAttribute(ReactomeJavaConstants.referenceEntity));
		infReferenceGeneProductInst.setAttribute(ReactomeJavaConstants.identifier, homologueId);
		// Reference DB can differ between homologue mappings, but can be differentiated by the
		// 'homologueSource' found in each mapping.
		// With PANTHER data, the Protein IDs are exclusively UniProt
		SimpleInstance referenceDatabaseInst = getReferenceDatabase().getReferenceDatabase(getHomologueSource(homologue));
		infReferenceGeneProductInst.setAttribute(ReactomeJavaConstants.referenceDatabase, referenceDatabaseInst);

		// Creates ReferenceDNASequence instance from ReferenceEntity
		List<SimpleInstance> inferredReferenceDNAInstances = createReferenceDNASequences(homologueId);
		infReferenceGeneProductInst.setAttribute(ReactomeJavaConstants.referenceGene, inferredReferenceDNAInstances);

		infReferenceGeneProductInst.setAttribute(ReactomeJavaConstants.species, getSpeciesInstance());
		String referenceGeneProductSource = getReferenceDatabase().getReferenceDatabaseSourceName(getHomologueSource(homologue));
		infReferenceGeneProductInst.setDisplayName(referenceGeneProductSource + ":" + homologueId);

		// GeneName value comes from UniProt's identifier mapping service.
		if (getMappings().homologueHasGeneName(homologueId)) {
			infReferenceGeneProductInst.setAttribute(ReactomeJavaConstants.geneName, getMappings().getGeneName(homologueId));
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

	private void inferStartAndEndCoordinates(SimpleInstance ewasInst, SimpleInstance infEWASInst) {
		// Method for adding start/end coordinates. It is convoluted due to a quirk with assigning the
		// name differently based on coordinate value (see infer_events.pl lines 1190-1192).
		// The name of the entity needs to be at the front of the 'name' array if the coordinate is over 1,
		// and rearranging arrays in Java for this was a bit tricky.
		List<Integer> inferredStartCoordinates = new ArrayList<>();
		for (int startCoord : (Collection<Integer>) ewasInst.getAttribute(ReactomeJavaConstants.startCoordinate)) {
			inferredStartCoordinates.add(startCoord);
		}
		infEWASInst.setAttribute(ReactomeJavaConstants.startCoordinate, inferredStartCoordinates);

		List<Integer> inferredEndCoordinates = new ArrayList<>();
		for (int endCoord : (Collection<Integer>) ewasInst.getAttribute(ReactomeJavaConstants.endCoordinate)) {
			inferredEndCoordinates.add(endCoord);
		}
		infEWASInst.setAttribute(ReactomeJavaConstants.endCoordinate, inferredEndCoordinates);

		if (infEWASInst.getAttribute(ReactomeJavaConstants.startCoordinate) != null &&
			(int) infEWASInst.getAttribute(ReactomeJavaConstants.startCoordinate) > 1 ||
			infEWASInst.getAttribute(ReactomeJavaConstants.endCoordinate) != null &&
				(int) infEWASInst.getAttribute(ReactomeJavaConstants.endCoordinate) > 1) {
			List<String> infEWASInstNames = (ArrayList<String>) (ewasInst).getAttribute(ReactomeJavaConstants.name);
			infEWASInst.setAttribute(ReactomeJavaConstants.name, infEWASInstNames.get(0));
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
	private List<SimpleInstance> createReferenceDNASequences(String homologueId) throws Exception {
		List<SimpleInstance> referenceDNAInstances = new ArrayList<>();
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

	private SimpleInstance createReferenceDNASequence(String ensgId) throws Exception {
		logger.info("Creating ReferenceDNASequence for " + ensgId);
		//SchemaClass referenceDNAClass = getCurrentDBA().getSchema().getClassByName(ReactomeJavaConstants.ReferenceDNASequence);
		SimpleInstance referenceDNAInst = new SimpleInstance();
		referenceDNAInst.setSchemaClassName(ReactomeJavaConstants.ReferenceDNASequence);
		referenceDNAInst.setAttribute(ReactomeJavaConstants.identifier, ensgId);
		referenceDNAInst.setAttribute(ReactomeJavaConstants.referenceDatabase, getReferenceDatabase().getEnsemblDBInst());
		referenceDNAInst.setAttribute(ReactomeJavaConstants.species, getSpeciesInstance());
		referenceDNAInst.setDisplayName("ENSEMBL:" + ensgId);
		referenceDNAInst = instanceUtilities.checkForIdenticalInstances(referenceDNAInst, null);
		return referenceDNAInst;
	}

	private SimpleInstance createReferenceDNASequenceWithAlternateReferenceDatabase(String ensgId) throws Exception {
		logger.info("Creating ReferenceDNASequence for " + ensgId + " using alternate reference database");
		SimpleInstance alternateRefDNAInst = new SimpleInstance();
		alternateRefDNAInst.setSchemaClassName(ReactomeJavaConstants.ReferenceDNASequence);

		String altDbIdentifier = ensgId;
		if (getReferenceDatabase().getAlternateReferenceDbId() != null) {
			altDbIdentifier = altDbIdentifier.replaceAll(getReferenceDatabase().getAlternateReferenceDbId(), "");
		}
		alternateRefDNAInst.setAttribute(ReactomeJavaConstants.identifier, altDbIdentifier);
		alternateRefDNAInst.setAttribute(ReactomeJavaConstants.referenceDatabase, getReferenceDatabase().getAlternateReferenceDatabase());
		alternateRefDNAInst.setAttribute(ReactomeJavaConstants.species, getSpeciesInstance());
		alternateRefDNAInst.setDisplayName(
			getReferenceDatabase().getAlternateReferenceDatabase().getAttribute(ReactomeJavaConstants.name) + ":" + ensgId);
		alternateRefDNAInst = instanceUtilities.checkForIdenticalInstances(alternateRefDNAInst, null);
		return alternateRefDNAInst;
	}

	private Mappings getMappings() {
		return this.mappings;
	}

	private OrthologyReferenceDatabase getReferenceDatabase() {
		return this.orthologyReferenceDatabase;
	}

	private SimpleInstance getSpeciesInstance() throws Exception {
		return this.instanceUtilities.getSpeciesInstance();
	}
}
