package org.reactome.orthoinference;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import org.gk.model.ReactomeJavaConstants;

import org.json.simple.parser.ParseException;
import org.reactome.curation.model.SimpleInstance;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

import static org.reactome.util.general.CollectionUtils.safeList;

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
	// How many RlEs are processed between each log of the curator-tool-ws call figures.
	private static final int CALL_STATS_INTERVAL = 100;

	private ConfigProperties configProperties;
	private String speciesCode;

	private ReactionInferrer reactionInferrer;
	private PathwaysInferrer pathwaysInferrer;

	private InstanceUtilities utils;
	private CuratorToolAPI curatorToolAPI;

	private static Map<Long,SimpleInstance> manualEventToNonHumanSource = new HashMap<>();
	private static List<SimpleInstance> manualHumanEvents = new ArrayList<>();
	private static OrthologousPathwayDiagramGenerator orthologousPathwayDiagramGenerator;

	public EventsInferrer(
		ConfigProperties configProperties,
		@Qualifier("targetSpeciesCode") String speciesCode,
		ReactionInferrer reactionInferrer,
		PathwaysInferrer pathwaysInferrer,
		InstanceUtilities utils,
		CuratorToolAPI curatorToolAPI
	) {
		this.configProperties = configProperties;
		this.speciesCode = speciesCode;

		this.reactionInferrer = reactionInferrer;
		this.pathwaysInferrer = pathwaysInferrer;

		this.utils = utils;
		this.curatorToolAPI = curatorToolAPI;
	}

	@SuppressWarnings("unchecked")
	public void inferEvents() throws Exception {
		logger.info("Beginning orthoinference of " + getSpeciesName());

		//SimpleInstance humanSpeciesInstance = getHumanSpeciesInstance();
		Map<Long, SimpleInstance> reactionMap = getHumanReactionMap();
		processReactions(reactionMap);
		generatePathways();
		//generateDiagrams();

		outputReport(speciesCode, configProperties.getReleaseVersion());
		logger.info("Finished orthoinference of " + getSpeciesName());
	}

	private Map<Long, SimpleInstance> getHumanReactionMap() {
		List<SimpleInstance> reactionInstances = getCuratorToolAPI().fetchHumanReactionlikeEvents();
		// curator-tool-ws reports a failed search as an empty result rather than an error, and there is always at
		// least one human reaction, so nothing found means the search failed.
		if (reactionInstances.isEmpty()) {
			logger.fatal("Could not find any Homo sapiens ReactionlikeEvent instances");
			System.exit(1);
		}

		// The search returns shell instances -- dbId, displayName and schemaClass only -- so each is inflated to have
		// the attributes the inference reads.
		Map<Long, SimpleInstance> reactionMap = new LinkedHashMap<>();
		for (SimpleInstance reactionShell : reactionInstances) {
			SimpleInstance reactionInst = getCuratorToolAPI().inflate(reactionShell);
			if (reactionInst == null) {
				throw new IllegalStateException("No ReactionlikeEvent found for db id " + reactionShell.getDbId() +
					" (" + reactionShell.getDisplayName() + ") returned by the human reaction search");
			}
			reactionMap.put(reactionInst.getDbId(), reactionInst);
		}

		logger.info("Homo sapiens ReactionlikeEvent instances: " + reactionMap.size());
		return reactionMap;
	}

	private void processReactions(Map<Long, SimpleInstance> reactionMap) throws Exception {
		int processedCount = 0;
		for (Map.Entry<Long, SimpleInstance> entry : reactionMap.entrySet()) {
			processedCount++;
			if (processedCount % CALL_STATS_INTERVAL == 0) {
				getCuratorToolAPI().logCallStats(processedCount + "/" + reactionMap.size() + " RlEs processed");
			}

			SimpleInstance reactionInst = entry.getValue();
			logger.info("Attempting RlE inference: " + reactionInst);

			if (isAlreadyInferred(reactionInst)) {
				continue;
			}

			inferReactionWithErrorHandling(reactionInst);
		}
		getCuratorToolAPI().logCallStats("All " + reactionMap.size() + " RlEs processed");
	}

	private boolean isAlreadyInferred(SimpleInstance reactionInst) throws Exception {
		List<SimpleInstance> previouslyInferredInstances = new ArrayList<>();
		previouslyInferredInstances.addAll(
			checkIfPreviouslyInferred(reactionInst, ReactomeJavaConstants.orthologousEvent, previouslyInferredInstances));
		previouslyInferredInstances.addAll(
			checkIfPreviouslyInferred(reactionInst, ReactomeJavaConstants.inferredFrom, previouslyInferredInstances));

		if (previouslyInferredInstances.isEmpty()) {
			return false;
		}

		return handlePreviouslyInferredInstance(reactionInst, previouslyInferredInstances.get(0));
	}

	private boolean handlePreviouslyInferredInstance(SimpleInstance reactionInst, SimpleInstance prevInfInst) {
		if (prevInfInst.getAttribute(ReactomeJavaConstants.disease) != null) {
			logger.info("Disease reaction, skipping inference");
			return true;
		}

		SimpleInstance evidenceTypeInst = (SimpleInstance) prevInfInst.getAttribute(ReactomeJavaConstants.evidenceType);
		if (evidenceTypeInst != null && evidenceTypeInst.getDisplayName().contains("electronic")) {
			getReactionInferrer().addAlreadyInferredEvents(reactionInst, prevInfInst);
		} else {
			logger.info("Inferred RlE already exists, skipping inference");
			manualEventToNonHumanSource.put(reactionInst.getDbId(), prevInfInst);
			manualHumanEvents.add(reactionInst);
		}
		return true;
	}

	private void inferReactionWithErrorHandling(SimpleInstance reactionInst) {
		try {
			getReactionInferrer().inferReaction(reactionInst);
			logger.info("Successfully inferred " + reactionInst);
		} catch (Exception e) {
			logger.error("Failed to infer reaction: " + reactionInst, e);
			throw new RuntimeException("Failed to infer reaction " + reactionInst, e);
		}
	}

	private void generatePathways() throws Exception {
		getPathwaysInferrer().setInferredEvent(getReactionInferrer().getInferredEvent());
		getPathwaysInferrer().inferPathways(getReactionInferrer().getInferrableHumanEvents());
	}

//	private void generateDiagrams() {
//		SimpleInstance humanSpeciesInstance = getHumanSpeciesInstance();
//		orthologousPathwayDiagramGenerator = new OrthologousPathwayDiagramGenerator(
//				getCurrentDBA(),
//				getPrevDBA(),
//				getSpeciesInstance(),
//				humanSpeciesInstance,
//				configProperties.getPersonId()
//		);
//		orthologousPathwayDiagramGenerator.generateOrthologousPathwayDiagrams();
//	}

	private SimpleInstance getHumanSpeciesInstance() {
		SimpleInstance humanSpeciesInstance = getCuratorToolAPI().fetchSpecies("Homo sapiens");

		if (humanSpeciesInstance == null) {
			logger.fatal("Could not find Species instance for Homo sapiens");
			System.exit(1);
		}
		return humanSpeciesInstance;
	}

	private void createNewFile(String filename) throws IOException {
		File file = new File(filename);
		if (file.exists()) {
			file.delete();
		}
		file.createNewFile();
	}

	@SuppressWarnings("unchecked")
	private List<SimpleInstance> checkIfPreviouslyInferred(
		SimpleInstance reactionInst, String attribute, List<SimpleInstance> previouslyInferredInstances)
		throws Exception {
		for (SimpleInstance attributeShell : safeList((List<SimpleInstance>) reactionInst.getAttribute(attribute))) {
			// The reaction's attribute values are shells -- dbId, displayName and schemaClass only -- so each is
			// inflated before its species can be read. The inflated instance is the one kept, since the caller goes
			// on to read its disease and evidenceType.
			SimpleInstance attributeInst = getCuratorToolAPI().inflate(attributeShell);
			if (attributeInst == null) {
				logger.warn("No instance found for " + attributeShell + " -- skipping it");
				continue;
			}
			if (hasSpecies(attributeInst, getSpeciesInstance()) &&
				attributeInst.getAttribute(ReactomeJavaConstants.isChimeric) == null) {
				previouslyInferredInstances.add(attributeInst);
			}
		}
		return previouslyInferredInstances;
	}

	// An Event's species is multi-valued, so its value is a list; a single instance is also accepted.
	private boolean hasSpecies(SimpleInstance instance, SimpleInstance speciesInst) {
		Object speciesValue = instance.getAttribute(ReactomeJavaConstants.species);
		Collection<?> speciesValues = speciesValue instanceof Collection ?
			(Collection<?>) speciesValue : Collections.singletonList(speciesValue);
		return speciesValues.stream().anyMatch(value -> value instanceof SimpleInstance &&
			Objects.equals(((SimpleInstance) value).getDbId(), speciesInst.getDbId()));
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

	private CuratorToolAPI getCuratorToolAPI() {
		return this.curatorToolAPI;
	}

	private String getSpeciesName() throws IOException, ParseException {
		return getConfigProperties().getSpeciesConfig().getSpeciesName(getSpeciesCode());
	}

	private SimpleInstance getSpeciesInstance() throws Exception {
		return getUtils().getSpeciesInstance();
	}

//	private MySQLAdaptor getCurrentDBA() throws SQLException {
//		return getConfigProperties().getCurrentDBA();
//	}
//
//	private MySQLAdaptor getPrevDBA() throws SQLException {
//		return getConfigProperties().getPreviousDBA();
//	}
}
