package org.reactome.orthoinference;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.sql.SQLException;
import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;

import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.InvalidAttributeException;

import org.json.simple.parser.ParseException;

/**
 *
 * @author jcook
 *
 * The Java version of infer_events.pl -- The gist of this module is that it looks at all existing Human
 * ReactionlikeEvent (RlE) instances (mostly Reactions and BlackBoxEvents) in the Test_Reactome database,
 * and attempts to computationally infer them in each of Reactome's model organisms. Each RlE is broken down into its
 * primary components (input, output, catalyst, and regulator), which are themselves broken into their PhysicalEntity
 * subunits. The homology data used for the inference process comes from PANTHER (www.pantherdb.org) and is generated
 * during the 'Orthopairs' step of the Reactome release process.  After all inference attempts for each RlE has been
 * completed in an organism, the pathways that contain the reactions are filled with these newly inferred ones.
 */

public class EventsInferrer {
	private static final Logger logger = LogManager.getLogger();

	private ConfigProperties configProperties;
	private String speciesCode;

	private ReactionInferrer reactionInferrer;
	private PathwaysInferrer pathwaysInferrer;

	private static GKInstance speciesInst;
	private static Map<GKInstance,GKInstance> manualEventToNonHumanSource = new HashMap<>();
	private static List<GKInstance> manualHumanEvents = new ArrayList<>();
	//private static StableIdentifierGenerator stableIdentifierGenerator;
	private static OrthologousPathwayDiagramGenerator orthologousPathwayDiagramGenerator;

	public EventsInferrer(ConfigProperties configProperties, String speciesCode) throws IOException {
		this.configProperties = configProperties;
		this.speciesCode = speciesCode;

		this.reactionInferrer = new ReactionInferrer(configProperties, speciesCode);
		this.pathwaysInferrer = new PathwaysInferrer();
	}

	@SuppressWarnings("unchecked")
	public void inferEvents() throws Exception {
		logger.info("Beginning orthoinference of " + getSpeciesName());

/**
 *  Start of ReactionlikeEvent inference. Retrieves all human ReactionlikeEvents, and attempts to infer each for the species.
 */
		// Gets DB instance of source species (human)
		GKInstance humanSpeciesInstance = getHumanSpeciesInstance();

		// Gets Reaction instances of source species (human)
		Collection<GKInstance> reactionInstances = (Collection<GKInstance>)
			getCurrentDBA().fetchInstanceByAttribute("ReactionlikeEvent", "species", "=", humanSpeciesInstance.getDBID());

		List<Long> dbids = new ArrayList<>();
		Map<Long, GKInstance> reactionMap = new HashMap<>();
		for (GKInstance reactionInst : reactionInstances) {
			dbids.add(reactionInst.getDBID());
			reactionMap.put(reactionInst.getDBID(), reactionInst);
		}
		Collections.sort(dbids);

		logger.info(humanSpeciesInstance.getDisplayName() + " ReactionlikeEvent instances: " + dbids.size());
		for (Long dbid : dbids) {
			GKInstance reactionInst = reactionMap.get(dbid);
			logger.info("Attempting RlE inference: " + reactionInst);
			// Check if the current Reaction already exists for this species, that it is a valid instance (passes
			// some filters), and that it doesn't have a Disease attribute.
			// Adds to manualHumanEvents array if it passes conditions. This code block allows you to re-run the code
			// without re-inferring instances.
			List<GKInstance> previouslyInferredInstances = new ArrayList<GKInstance>();
			previouslyInferredInstances.addAll(
				checkIfPreviouslyInferred(reactionInst, ReactomeJavaConstants.orthologousEvent, previouslyInferredInstances));
			previouslyInferredInstances.addAll(
				checkIfPreviouslyInferred(reactionInst, ReactomeJavaConstants.inferredFrom, previouslyInferredInstances));
			if (previouslyInferredInstances.size() > 0) {
				GKInstance prevInfInst = previouslyInferredInstances.get(0);
				if (prevInfInst.getAttributeValue(ReactomeJavaConstants.disease) == null) {
					GKInstance evidenceTypeInst = (GKInstance) prevInfInst.getAttributeValue(ReactomeJavaConstants.evidenceType);
					if (evidenceTypeInst != null &&
						evidenceTypeInst.getDisplayName().contains("electronic")) {
						getReactionInferrer().addAlreadyInferredEvents(reactionInst, prevInfInst);
					} else {
						logger.info("Inferred RlE already exists, skipping inference");
						manualEventToNonHumanSource.put(reactionInst, prevInfInst);
						manualHumanEvents.add(reactionInst);
					}
				} else {
					logger.info("Disease reaction, skipping inference");
				}
				continue;
			}

			// An inferred ReactionlikeEvent doesn't already exist for this species, and an orthologous inference will
			// be attempted.
			try {
				getReactionInferrer().inferReaction(reactionInst);
				logger.info("Successfully inferred " + reactionInst);
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		PathwaysInferrer.setInferredEvent(getReactionInferrer().getInferredEvent());
		PathwaysInferrer.inferPathways(getReactionInferrer().getInferrableHumanEvents());

		orthologousPathwayDiagramGenerator = new OrthologousPathwayDiagramGenerator(
			getCurrentDBA(), getPrevDBA(), speciesInst, humanSpeciesInstance, configProperties.getPersonId());
		orthologousPathwayDiagramGenerator.generateOrthologousPathwayDiagrams();
		outputReport(speciesCode, configProperties.getReleaseVersion());
		logger.info("Finished orthoinference of " + getSpeciesName());
	}

//	public StableIdentifierGenerator getStableIdentifierGenerator() {
//		return stableIdentifierGenerator;
//	}

	private GKInstance getHumanSpeciesInstance() throws Exception {
		Collection<GKInstance> sourceSpeciesInst = (Collection<GKInstance>)
			getCurrentDBA().fetchInstanceByAttribute("Species", "name", "=", "Homo sapiens");
		if (sourceSpeciesInst.isEmpty()) {
			logger.fatal("Could not find Species instance for Homo sapiens");
			System.exit(1);
		}
		return sourceSpeciesInst.iterator().next();
	}

//	/**
//	 * Create mapping of UniProt accessions to species-specific gene names, and then set this mapping for use in
//	 * EWASInferrer.
//	 * @param species String - 4-letter shortened version of species name (eg: Homo sapiens --> hsap).
//	 * @throws IOException - Thrown if file is not found.
//	 */
//	private void readAndSetGeneNameMappingFile(String species) throws IOException {
//		Map<String, String> geneNameMappings = readGeneNameMappingFile(species);
//
//		EWASInferrer.setGeneNameMappingFile(geneNameMappings);
//	}

	private void createNewFile(String filename) throws IOException {
		File file = new File(filename);
		if (file.exists()) {
			file.delete();
		}
		file.createNewFile();
	}

//	private void setReleaseDates(String dateOfRelease) {
//		ReactionInferrer.setReleaseDate(dateOfRelease);
//		PathwaysInferrer.setReleaseDate(dateOfRelease);
//	}

	@SuppressWarnings("unchecked")
	private List<GKInstance> checkIfPreviouslyInferred(
		GKInstance reactionInst, String attribute, List<GKInstance> previouslyInferredInstances)
		throws InvalidAttributeException, Exception {
		for (GKInstance attributeInst : (Collection<GKInstance>) reactionInst.getAttributeValuesList(attribute)) {
			GKInstance reactionSpeciesInst = (GKInstance) attributeInst.getAttributeValue(ReactomeJavaConstants.species);
			if (reactionSpeciesInst.getDBID() == speciesInst.getDBID() &&
				attributeInst.getAttributeValue(ReactomeJavaConstants.isChimeric) == null) {
				previouslyInferredInstances.add(attributeInst);
			}
		}
		return previouslyInferredInstances;
	}

	private void outputReport(String species, String releaseVersion) throws IOException {
		int eligibleCount = getReactionInferrer().getEligibleCount();
		int inferredCount = getReactionInferrer().getInferredCount();
		float percentInferred = (float) 100 * inferredCount / eligibleCount;
		// Create file if it doesn't exist
		String reportFilename = "report_ortho_inference_test_reactome_" + releaseVersion + ".txt";
		logger.info("Updating " + reportFilename);
		if (!Files.exists(Paths.get(reportFilename))) {
			createNewFile(reportFilename);
		}
		String results = "hsap to " + species + ":\t" + inferredCount + " out of " + eligibleCount +
			" eligible reactions (" + String.format("%.2f", percentInferred) + "%)\n";
		Files.write(Paths.get(reportFilename), results.getBytes(), StandardOpenOption.APPEND);
	}

//	// Statically store the adaptor variable in each class
//	private void setDbAdaptors(MySQLAdaptor dbAdaptor) {
//		ReactionInferrer.setAdaptor(dbAdaptor);
//		InstanceUtilities.setAdaptor(dbAdaptor);
//		OrthologousEntityGenerator.setAdaptor(dbAdaptor);
//		EWASInferrer.setAdaptor(dbAdaptor);
//		PathwaysInferrer.setAdaptor(dbAdaptor);
//
//	}

//	private void readAndSetHomologueMappingFile(String species, String fromSpecies)
//		throws IOException {
//		Map<String,String[]> homologueMappings = readHomologueMappingFile(species, fromSpecies);
//		ProteinCountUtility.setHomologueMappingFile(homologueMappings);
//		EWASInferrer.setHomologueMappingFile(homologueMappings);
//	}

	private ConfigProperties getConfigProperties() {
		return this.configProperties;
	}

	private String getSpeciesCode() {
		return this.speciesCode;
	}

	private ReactionInferrer getReactionInferrer() {
		return this.reactionInferrer;
	}

	private String getSpeciesName() throws IOException, ParseException {
		return getConfigProperties().getSpeciesConfig().getSpeciesName(getSpeciesCode());
	}

	private MySQLAdaptor getCurrentDBA() throws SQLException {
		return getConfigProperties().getCurrentDBA();
	}

	private MySQLAdaptor getPrevDBA() throws SQLException {
		return getConfigProperties().getPreviousDBA();
	}
}
