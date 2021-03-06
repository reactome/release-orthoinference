package org.reactome.orthoinference;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.*;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.FileUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;
import static org.gk.model.ReactomeJavaConstants.*;

import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.InvalidAttributeException;
import org.gk.schema.SchemaClass;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.reactome.release.common.database.InstanceEditUtils;

/**
 *
 * @author jcook
 *
 * The Java version of infer_events.pl -- The gist of this module is that it looks at all existing Human ReactionlikeEvent (RlE) instances (mostly Reactions and BlackBoxEvents) in the Test_Reactome database,
 * and attempts to computationally infer them in each of Reactome's model organisms. Each RlE is broken down into its primary components (input, output, catalyst, and regulator), which are themselves broken
 * into their PhysicalEntity subunits. The homology data used for the inference process comes from PANTHER (www.pantherdb.org) and is generated during the 'Orthopairs' step of the Reactome release process.
 * After all inference attempts for each RlE has been completed in an organism, the pathways that contain the reactions are filled with these newly inferred ones.
 *
 *
 */

public class EventsInferrer
{
	private static final Logger logger = LogManager.getLogger();
	private static MySQLAdaptor dbAdaptor;
	private static MySQLAdaptor dbAdaptorPrev;
	private static String releaseVersion;
	private static GKInstance instanceEditInst;
	private static GKInstance speciesInst;
	private static Map<GKInstance,GKInstance> manualEventToNonHumanSource = new HashMap<>();
	private static List<GKInstance> manualHumanEvents = new ArrayList<>();
	private static StableIdentifierGenerator stableIdentifierGenerator;
	private static OrthologousPathwayDiagramGenerator orthologousPathwayDiagramGenerator;
	private static final String INFERRED_EVIDENCE_TYPE_DISPLAY_NAME = "inferred by electronic annotation";

	@SuppressWarnings("unchecked")
	public static void inferEvents(Properties props, String species) throws Exception
	{
		logger.info("Preparing DB Adaptor and setting project variables");
		// Set up DB adaptor using config.properties file
		String username = props.getProperty("release.database.user");
		String password = props.getProperty("release.database.password");
		String database = props.getProperty("release_current.name");
		String prevDatabase = props.getProperty("release_previous.name");
		String host = props.getProperty("release.database.host");
		int port = Integer.valueOf(props.getProperty("release.database.port"));

		dbAdaptor = new MySQLAdaptor(host, database, username, password, port);
		dbAdaptorPrev = new MySQLAdaptor(host, prevDatabase, username, password, port);
		if (dbAdaptor == null || dbAdaptorPrev == null) {
			logger.fatal("Null MySQLAdaptor, terminating orthoinference");
			System.exit(1);
		}
		setDbAdaptors(dbAdaptor);

		releaseVersion = props.getProperty("releaseNumber");
		String pathToOrthopairs = props.getProperty("pathToOrthopairs", "orthopairs");
		String pathToSpeciesConfig = props.getProperty("pathToSpeciesConfig", "src/main/resources/Species.json");
		String dateOfRelease = props.getProperty("dateOfRelease");
		int personId = Integer.valueOf(props.getProperty("personId"));
		setReleaseDates(dateOfRelease);

		// Finds all skippable ReactionlikeEvents based on a static list of skippable Pathways.
		SkipInstanceChecker.buildStaticSkipList();

		JSONParser parser = new JSONParser();
		Object obj = parser.parse(new FileReader(pathToSpeciesConfig));
		JSONObject jsonObject = (JSONObject) obj;

		// Parse Species information (found in Species.json config file)
		JSONObject speciesObject = (JSONObject) jsonObject.get(species);
		JSONArray speciesNames = (JSONArray) speciesObject.get("name");
		String speciesName = (String) speciesNames.get(0);
		logger.info("Beginning orthoinference of " + speciesName);

		JSONObject refDb = (JSONObject) speciesObject.get("refdb");
		String refDbUrl = (String) refDb.get("url");
		String refDbProteinUrl = (String) refDb.get("access");
		String refDbGeneUrl = (String) refDb.get("ensg_access");

		// Creates two files that a) list reactions that are eligible for inference and b) those that are successfully inferred
		String eligibleFilename = "eligible_" + species	+ "_75.txt";
		String inferredFilename = "inferred_" + species + "_75.txt";
		createNewFile(eligibleFilename);
		createNewFile(inferredFilename);
		ReactionInferrer.setEligibleFilename(eligibleFilename);
		ReactionInferrer.setInferredFilename(inferredFilename);

		stableIdentifierGenerator = new StableIdentifierGenerator(dbAdaptor, (String) speciesObject.get("abbreviation"));
		// Set static variables (DB/Species Instances, mapping files) that will be repeatedly used
		setInstanceEdits(personId);
		try {
			readAndSetHomologueMappingFile(species, "hsap", pathToOrthopairs);
			readAndSetGeneNameMappingFile(species, pathToOrthopairs);
		} catch (Exception e) {
			logger.fatal("Unable to locate " + speciesName +" mapping file: hsap_" + species + "_mapping.tsv. Orthology prediction not possible.");
			e.printStackTrace();
			System.exit(1);
		}
		EWASInferrer.readENSGMappingFile(species, pathToOrthopairs);
		EWASInferrer.fetchAndSetUniprotDbInstance();
		EWASInferrer.createEnsemblProteinDbInstance(speciesName, refDbUrl, refDbProteinUrl);
		EWASInferrer.createEnsemblGeneDBInstance(speciesName, refDbUrl, refDbGeneUrl);

		JSONObject altRefDbJSON = (JSONObject) speciesObject.get("alt_refdb");
		if (altRefDbJSON != null)
		{
			logger.info("Alternate DB exists for " + speciesName);
			EWASInferrer.createAlternateReferenceDBInstance(altRefDbJSON);
		} else {
			EWASInferrer.setAltRefDbToFalse();
		}
		createAndSetSpeciesInstance(speciesName);
		setSummationInstance();
		setEvidenceTypeInstance();
		OrthologousEntityGenerator.setComplexSummationInstance();

/**
 *  Start of ReactionlikeEvent inference. Retrieves all human ReactionlikeEvents, and attempts to infer each for the species.
 */
		// Gets DB instance of source species (human)
		Collection<GKInstance> sourceSpeciesInst = (Collection<GKInstance>) dbAdaptor.fetchInstanceByAttribute("Species", "name", "=", "Homo sapiens");
		if (sourceSpeciesInst.isEmpty())
		{
			logger.fatal("Could not find Species instance for Homo sapiens");
			System.exit(1);
		}
		long humanInstanceDbId = sourceSpeciesInst.iterator().next().getDBID();
		orthologousPathwayDiagramGenerator = new OrthologousPathwayDiagramGenerator(dbAdaptor, dbAdaptorPrev, speciesInst, personId, humanInstanceDbId);
		// Gets Reaction instances of source species (human)
		Collection<GKInstance> reactionInstances = (Collection<GKInstance>) dbAdaptor.fetchInstanceByAttribute("ReactionlikeEvent", "species", "=", humanInstanceDbId);

		List<Long> dbids = new ArrayList<>();
		Map<Long, GKInstance> reactionMap = new HashMap<>();
		for (GKInstance reactionInst : reactionInstances) {
			dbids.add(reactionInst.getDBID());
			reactionMap.put(reactionInst.getDBID(), reactionInst);
		}
		Collections.sort(dbids);

		logger.info(sourceSpeciesInst.iterator().next().getDisplayName() + " ReactionlikeEvent instances: " + dbids.size());
		for (Long dbid : dbids)
		{
			GKInstance reactionInst = reactionMap.get(dbid);
			logger.info("Attempting RlE inference: " + reactionInst);
			// Check if the current Reaction already exists for this species, that it is a valid instance (passes some filters), and that it doesn't have a Disease attribute.
			// Adds to manualHumanEvents array if it passes conditions. This code block allows you to re-run the code without re-inferring instances.
			List<GKInstance> previouslyInferredInstances = new ArrayList<GKInstance>();
			previouslyInferredInstances.addAll(checkIfPreviouslyInferred(reactionInst, orthologousEvent, previouslyInferredInstances));
			previouslyInferredInstances.addAll(checkIfPreviouslyInferred(reactionInst, inferredFrom, previouslyInferredInstances));
			if (previouslyInferredInstances.size() > 0)
			{
				GKInstance prevInfInst = previouslyInferredInstances.get(0);
				if (prevInfInst.getAttributeValue(disease) == null)
				{
					GKInstance evidenceTypeInst = (GKInstance) prevInfInst.getAttributeValue(evidenceType);
					if (evidenceTypeInst != null && evidenceTypeInst.getDisplayName().contains(INFERRED_EVIDENCE_TYPE_DISPLAY_NAME)) {
						ReactionInferrer.addAlreadyInferredEvents(reactionInst, prevInfInst);
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

			// An inferred ReactionlikeEvent doesn't already exist for this species, and an orthologous inference will be attempted.
			try {
				ReactionInferrer.inferReaction(reactionInst);
				logger.info("Successfully inferred " + reactionInst);
			} catch (Exception e) {
				e.printStackTrace();
				System.exit(1);
			}
		}
		PathwaysInferrer.setInferredEvent(ReactionInferrer.getInferredEvent());
		PathwaysInferrer.inferPathways(ReactionInferrer.getInferrableHumanEvents());
		orthologousPathwayDiagramGenerator.generateOrthologousPathwayDiagrams();
		outputReport(species);
		logger.info("Finished orthoinference of " + speciesName);
	}

	/**
	 * Create mapping of UniProt accessions to species-specific gene names, and then set this mapping for use in EWASInferrer.
	 * @param species String - 4-letter shortened version of species name (eg: Homo sapiens --> hsap).
	 * @param pathToOrthopairs String - Path to directory containing orthopairs files.
	 * @throws IOException - Thrown if file is not found.
	 */
	private static void readAndSetGeneNameMappingFile(String species, String pathToOrthopairs) throws IOException {
		Map<String, String> geneNameMappings = readGeneNameMappingFile(species, pathToOrthopairs);

		EWASInferrer.setGeneNameMappingFile(geneNameMappings);
	}

	/**
	 * Read in the {species}_gene_name_mapping.tsv file and create a Map of UniProt identifiers to gene names.
	 * @param species String - 4-letter shortened version of species name (eg: Homo sapiens --> hsap).
	 * @param pathToOrthopairs String - Path to directory containing orthopairs files.
	 * @return pathToOrthopairs String - Path to directory containing orthopairs files.
	 * @throws IOException - Thrown if file is not found.
	 */
	private static Map<String, String> readGeneNameMappingFile(String species, String pathToOrthopairs) throws IOException {
		Path geneNameMappingFilePath = Paths.get(pathToOrthopairs, species + "_gene_name_mapping.tsv");
		Map<String, String> geneNameMappings = new HashMap<>();
		for (String line : Files.readAllLines(geneNameMappingFilePath)) {
			String[] tabSplit = line.split("\t");
			if (tabSplit.length == 2) {
				String uniprotId = tabSplit[0];
				String geneName = tabSplit[1];
				geneNameMappings.put(uniprotId, geneName);
			}
		}

		return geneNameMappings;
	}

	private static void createNewFile(String filename) throws IOException {
		File file = new File(filename);
		if (file.exists()) {
			file.delete();
		}
		file.createNewFile();
	}

	public static StableIdentifierGenerator getStableIdentifierGenerator() {
		return stableIdentifierGenerator;
	}

	private static void setReleaseDates(String dateOfRelease)
	{
		ReactionInferrer.setReleaseDate(dateOfRelease);
		PathwaysInferrer.setReleaseDate(dateOfRelease);

	}

	@SuppressWarnings("unchecked")
	private static List<GKInstance> checkIfPreviouslyInferred(GKInstance reactionInst, String attribute, List<GKInstance> previouslyInferredInstances) throws InvalidAttributeException, Exception
	{
		for (GKInstance attributeInst : (Collection<GKInstance>) reactionInst.getAttributeValuesList(attribute))
		{
			GKInstance reactionSpeciesInst = (GKInstance) attributeInst.getAttributeValue(species);
			if (reactionSpeciesInst.getDBID() == speciesInst.getDBID() && attributeInst.getAttributeValue(isChimeric) == null)
			{
				previouslyInferredInstances.add(attributeInst);
			}
		}
		return previouslyInferredInstances;
	}

	private static void outputReport(String species) throws IOException
	{
		int eligibleCount = ReactionInferrer.getEligibleCount();
		int inferredCount = ReactionInferrer.getInferredCount();
		float percentInferred = (float) 100*inferredCount/eligibleCount;
		// Create file if it doesn't exist
		String reportFilename = "report_ortho_inference_test_reactome_" + releaseVersion + ".txt";
		logger.info("Updating " + reportFilename);
		if (!Files.exists(Paths.get(reportFilename))) {
			createNewFile(reportFilename);
		}
		String results = "hsap to " + species + ":\t" + inferredCount + " out of " + eligibleCount + " eligible reactions (" + String.format("%.2f", percentInferred) + "%)\n";
		Files.write(Paths.get(reportFilename), results.getBytes(), StandardOpenOption.APPEND);
	}

	// Statically store the adaptor variable in each class
	private static void setDbAdaptors(MySQLAdaptor dbAdaptor)
	{
		ReactionInferrer.setAdaptor(dbAdaptor);
		SkipInstanceChecker.setAdaptor(dbAdaptor);
		InstanceUtilities.setAdaptor(dbAdaptor);
		OrthologousEntityGenerator.setAdaptor(dbAdaptor);
		EWASInferrer.setAdaptor(dbAdaptor);
		PathwaysInferrer.setAdaptor(dbAdaptor);

	}

	private static void readAndSetHomologueMappingFile(String species, String fromSpecies, String pathToOrthopairs) throws IOException {
		Map<String,String[]> homologueMappings = readHomologueMappingFile(species, fromSpecies, pathToOrthopairs);
		ProteinCountUtility.setHomologueMappingFile(homologueMappings);
		EWASInferrer.setHomologueMappingFile(homologueMappings);
	}

	// Read the species-specific orthopair 'mapping' file, and create a HashMap with the contents
	private static Map<String, String[]> readHomologueMappingFile(String toSpecies, String fromSpecies, String pathToOrthopairs) throws IOException
	{
		String orthopairsFileName = fromSpecies + "_" + toSpecies + "_mapping.tsv";
		String orthopairsFilePath = Paths.get(pathToOrthopairs, orthopairsFileName).toString();
		logger.info("Reading in " + orthopairsFilePath);
		FileReader fr = new FileReader(orthopairsFilePath);
		BufferedReader br = new BufferedReader(fr);

		Map<String, String[]> homologueMappings = new HashMap<>();
		String currentLine;
		while ((currentLine = br.readLine()) != null)
		{
			String[] tabSplit = currentLine.split("\t");
			String mapKey = tabSplit[0];
			String[] spaceSplit = tabSplit[1].split(" ");
			homologueMappings.put(mapKey, spaceSplit);
		}
		br.close();
		fr.close();
		return homologueMappings;
	}

	// Find the instance specific to this species
	private static void createAndSetSpeciesInstance(String toSpeciesLong) throws Exception
	{
		SchemaClass referenceDb = dbAdaptor.getSchema().getClassByName(Species);
		speciesInst = new GKInstance(referenceDb);
		speciesInst.setDbAdaptor(dbAdaptor);
		speciesInst.addAttributeValue(created, instanceEditInst);
		speciesInst.addAttributeValue(name, toSpeciesLong);
		speciesInst.addAttributeValue(_displayName, toSpeciesLong);
		speciesInst = InstanceUtilities.checkForIdenticalInstances(speciesInst, null);
		logger.info("Using species instance: " + speciesInst);
		OrthologousEntityGenerator.setSpeciesInstance(speciesInst);
		EWASInferrer.setSpeciesInstance(speciesInst);
		InstanceUtilities.setSpeciesInstance(speciesInst);
	}
	// Create and set static Summation instance
	private static void setSummationInstance() throws Exception
	{
		GKInstance summationInst = new GKInstance(dbAdaptor.getSchema().getClassByName(Summation));
		summationInst.setDbAdaptor(dbAdaptor);
		summationInst.addAttributeValue(created, instanceEditInst);
		String summationText = "This event has been computationally inferred from an event that has been demonstrated in another species.<p>The inference is based on the homology mapping from PANTHER. Briefly, reactions for which all involved PhysicalEntities (in input, output and catalyst) have a mapped orthologue/paralogue (for complexes at least 75% of components must have a mapping) are inferred to the other species. High level events are also inferred for these events to allow for easier navigation.<p><a href='/electronic_inference_compara.html' target = 'NEW'>More details and caveats of the event inference in Reactome.</a> For details on PANTHER see also: <a href='http://www.pantherdb.org/about.jsp' target='NEW'>http://www.pantherdb.org/about.jsp</a>";
		summationInst.addAttributeValue(text, summationText);
		summationInst.addAttributeValue(_displayName, summationText);
		summationInst = InstanceUtilities.checkForIdenticalInstances(summationInst, null);
		ReactionInferrer.setSummationInstance(summationInst);
		PathwaysInferrer.setSummationInstance(summationInst);
	}
	// Create and set static EvidenceType instance
	private static void setEvidenceTypeInstance() throws Exception
	{
		GKInstance evidenceTypeInst = new GKInstance(dbAdaptor.getSchema().getClassByName(EvidenceType));
		evidenceTypeInst.setDbAdaptor(dbAdaptor);
		evidenceTypeInst.addAttributeValue(created, instanceEditInst);
		String evidenceTypeText = INFERRED_EVIDENCE_TYPE_DISPLAY_NAME;
		evidenceTypeInst.addAttributeValue(name, evidenceTypeText);
		evidenceTypeInst.addAttributeValue(name, "IEA");
		evidenceTypeInst.addAttributeValue(_displayName, evidenceTypeText);
		evidenceTypeInst = InstanceUtilities.checkForIdenticalInstances(evidenceTypeInst, null);
		ReactionInferrer.setEvidenceTypeInstance(evidenceTypeInst);
		PathwaysInferrer.setEvidenceTypeInstance(evidenceTypeInst);
	}

	private static void setInstanceEdits(int personId) throws Exception
	{
		instanceEditInst = InstanceEditUtils.createInstanceEdit(dbAdaptor, personId, "org.reactome.orthoinference");
		logger.info("Instance edit: " + instanceEditInst);
		InstanceUtilities.setInstanceEdit(instanceEditInst);
		OrthologousEntityGenerator.setInstanceEdit(instanceEditInst);
		EWASInferrer.setInstanceEdit(instanceEditInst);
		PathwaysInferrer.setInstanceEdit(instanceEditInst);
	}

	/**
	 * Download the Wormbase file that contains the WBGene IDs mapped to gene names, and then create a mapping of IDs to names.
	 * @param wormbaseUrl -- URL to Wormbase file to be downloaded and processed.
	 * @return -- Map<String, List,<String>> of WBGene IDs to gene names.
	 * @throws IOException
	 */
	private static Map<String, List<String>> downloadAndProcessWormbaseFile(URL wormbaseUrl) throws IOException {
		File wormbaseFile = Paths.get(wormbaseUrl.toString()).getFileName().toFile();
		Map<String, List<String>> wormbaseMappings = new HashMap<>();
		FileUtils.copyURLToFile(wormbaseUrl, wormbaseFile);
		BufferedReader br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(wormbaseFile))));
		String line = null;
		while ((line = br.readLine()) != null) {
			String wormbaseGeneId = line.split(",")[1];
			String wormbaseGeneName = line.split(",")[2];
			if (!wormbaseGeneId.isEmpty() && !wormbaseGeneName.isEmpty()) {
				if (wormbaseMappings.containsKey(wormbaseGeneId)) {
					wormbaseMappings.get(wormbaseGeneId).add(wormbaseGeneName);
				} else {
					wormbaseMappings.computeIfAbsent(wormbaseGeneId, k -> new ArrayList<>()).add(wormbaseGeneName);
				}
			}
		}
		br.close();
		wormbaseFile.delete();
		return wormbaseMappings;
	}
}
