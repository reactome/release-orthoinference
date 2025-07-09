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
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

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
@Component
public class EventsInferrer {
	private static final Logger logger = LogManager.getLogger();

	private ConfigProperties configProperties;
	private String speciesCode;

	private ReactionInferrer reactionInferrer;
	private PathwaysInferrer pathwaysInferrer;

	private InstanceUtilities utils;

	private static Map<GKInstance,GKInstance> manualEventToNonHumanSource = new HashMap<>();
	private static List<GKInstance> manualHumanEvents = new ArrayList<>();
	private static OrthologousPathwayDiagramGenerator orthologousPathwayDiagramGenerator;

	public EventsInferrer(
		ConfigProperties configProperties,
		@Qualifier("targetSpeciesCode") String speciesCode,
		ReactionInferrer reactionInferrer,
		PathwaysInferrer pathwaysInferrer,
		InstanceUtilities utils
	) {
		this.configProperties = configProperties;
		this.speciesCode = speciesCode;

		this.reactionInferrer = reactionInferrer;
		this.pathwaysInferrer = pathwaysInferrer;

		this.utils = utils;
	}

	@SuppressWarnings("unchecked")
	public void inferEvents() throws Exception {
		logger.info("Beginning orthoinference of " + getSpeciesName());

		GKInstance humanSpeciesInstance = getHumanSpeciesInstance();
		Map<Long, GKInstance> reactionMap = getHumanReactionMap(humanSpeciesInstance);
		processReactions(reactionMap);
		generatePathwaysAndDiagrams(humanSpeciesInstance);

		outputReport(speciesCode, configProperties.getReleaseVersion());
		logger.info("Finished orthoinference of " + getSpeciesName());
	}

	private Map<Long, GKInstance> getHumanReactionMap(GKInstance humanSpeciesInstance) throws Exception {
		Collection<GKInstance> reactionInstances = (Collection<GKInstance>)
				getCurrentDBA().fetchInstanceByAttribute("ReactionlikeEvent", "species", "=", humanSpeciesInstance.getDBID());

		Map<Long, GKInstance> reactionMap = new HashMap<>();
		List<Long> dbids = new ArrayList<>();

		for (GKInstance reactionInst : reactionInstances) {
			dbids.add(reactionInst.getDBID());
			reactionMap.put(reactionInst.getDBID(), reactionInst);
		}
		Collections.sort(dbids);

		logger.info(humanSpeciesInstance.getDisplayName() + " ReactionlikeEvent instances: " + dbids.size());
		return reactionMap;
	}

	private void processReactions(Map<Long, GKInstance> reactionMap) throws Exception {
		for (Map.Entry<Long, GKInstance> entry : reactionMap.entrySet()) {
			GKInstance reactionInst = entry.getValue();
			logger.info("Attempting RlE inference: " + reactionInst);

			if (isAlreadyInferred(reactionInst)) {
				continue;
			}

			inferReactionWithErrorHandling(reactionInst);
		}
	}

	private boolean isAlreadyInferred(GKInstance reactionInst) throws Exception {
		List<GKInstance> previouslyInferredInstances = new ArrayList<>();
		previouslyInferredInstances.addAll(
				checkIfPreviouslyInferred(reactionInst, ReactomeJavaConstants.orthologousEvent, previouslyInferredInstances));
		previouslyInferredInstances.addAll(
				checkIfPreviouslyInferred(reactionInst, ReactomeJavaConstants.inferredFrom, previouslyInferredInstances));

		if (previouslyInferredInstances.isEmpty()) {
			return false;
		}

		return handlePreviouslyInferredInstance(reactionInst, previouslyInferredInstances.get(0));
	}

	private boolean handlePreviouslyInferredInstance(GKInstance reactionInst, GKInstance prevInfInst) throws Exception {
		if (prevInfInst.getAttributeValue(ReactomeJavaConstants.disease) != null) {
			logger.info("Disease reaction, skipping inference");
			return true;
		}

		GKInstance evidenceTypeInst = (GKInstance) prevInfInst.getAttributeValue(ReactomeJavaConstants.evidenceType);
		if (evidenceTypeInst != null && evidenceTypeInst.getDisplayName().contains("electronic")) {
			getReactionInferrer().addAlreadyInferredEvents(reactionInst, prevInfInst);
		} else {
			logger.info("Inferred RlE already exists, skipping inference");
			manualEventToNonHumanSource.put(reactionInst, prevInfInst);
			manualHumanEvents.add(reactionInst);
		}
		return true;
	}

	private void inferReactionWithErrorHandling(GKInstance reactionInst) {
		try {
			getReactionInferrer().inferReaction(reactionInst);
			logger.info("Successfully inferred " + reactionInst);
		} catch (Exception e) {
			logger.error("Failed to infer reaction: " + reactionInst, e);
			throw new RuntimeException("Failed to infer reaction", e);
		}
	}

	private void generatePathwaysAndDiagrams(GKInstance humanSpeciesInstance) throws Exception {
		getPathwaysInferrer().setInferredEvent(getReactionInferrer().getInferredEvent());
		getPathwaysInferrer().inferPathways(getReactionInferrer().getInferrableHumanEvents());

		orthologousPathwayDiagramGenerator = new OrthologousPathwayDiagramGenerator(
				getCurrentDBA(),
				getPrevDBA(),
				getSpeciesInstance(),
				humanSpeciesInstance,
				configProperties.getPersonId()
		);
		orthologousPathwayDiagramGenerator.generateOrthologousPathwayDiagrams();
	}

	private GKInstance getHumanSpeciesInstance() throws Exception {
		Collection<GKInstance> sourceSpeciesInst = (Collection<GKInstance>)
			getCurrentDBA().fetchInstanceByAttribute("Species", "name", "=", "Homo sapiens");
		if (sourceSpeciesInst.isEmpty()) {
			logger.fatal("Could not find Species instance for Homo sapiens");
			System.exit(1);
		}
		return sourceSpeciesInst.iterator().next();
	}

	private void createNewFile(String filename) throws IOException {
		File file = new File(filename);
		if (file.exists()) {
			file.delete();
		}
		file.createNewFile();
	}

	@SuppressWarnings("unchecked")
	private List<GKInstance> checkIfPreviouslyInferred(
		GKInstance reactionInst, String attribute, List<GKInstance> previouslyInferredInstances)
		throws InvalidAttributeException, Exception {
		for (GKInstance attributeInst : (Collection<GKInstance>) reactionInst.getAttributeValuesList(attribute)) {
			GKInstance reactionSpeciesInst = (GKInstance) attributeInst.getAttributeValue(ReactomeJavaConstants.species);
			if (reactionSpeciesInst.getDBID() == getSpeciesInstance().getDBID() &&
				attributeInst.getAttributeValue(ReactomeJavaConstants.isChimeric) == null) {
				previouslyInferredInstances.add(attributeInst);
			}
		}
		return previouslyInferredInstances;
	}

	private void outputReport(String species, int releaseVersion) throws IOException {
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

	private ConfigProperties getConfigProperties() {
		return this.configProperties;
	}

	private String getSpeciesCode() {
		return this.speciesCode;
	}

	private ReactionInferrer getReactionInferrer() {
		return this.reactionInferrer;
	}

	private PathwaysInferrer getPathwaysInferrer() {
		return this.pathwaysInferrer;
	}

	private InstanceUtilities getUtils() {
		return this.utils;
	}

	private String getSpeciesName() throws IOException, ParseException {
		return getConfigProperties().getSpeciesConfig().getSpeciesName(getSpeciesCode());
	}

	private GKInstance getSpeciesInstance() throws Exception {
		return getUtils().getSpeciesInstance();
	}

	private MySQLAdaptor getCurrentDBA() throws SQLException {
		return getConfigProperties().getCurrentDBA();
	}

	private MySQLAdaptor getPrevDBA() throws SQLException {
		return getConfigProperties().getPreviousDBA();
	}
}
