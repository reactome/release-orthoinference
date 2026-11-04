package org.reactome.orthoinference;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.ReactomeJavaConstants;
import org.reactome.curation.model.SimpleInstance;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

import static org.reactome.util.general.CollectionUtils.safeList;

@Component
public class ReactionInferrer {

	private static final Logger logger = LogManager.getLogger();
	private final CuratorToolAPI curatorToolAPI;

	private ConfigProperties configProperties;
	private String speciesCode;
	private OrthologousEntityGenerator orthologousEntityGenerator;
	private StableIdentifierGenerator stableIdentifierGenerator;
	private InstanceUtilities instanceUtilities;
	private ProteinCountUtility proteinCountUtility;

	private static Map<Long, SimpleInstance> inferredCatalystMap = new HashMap<>();
	private static Map<Long, SimpleInstance> inferredEventMap = new HashMap<>();
	private static Integer eligibleCount = 0;
	private static Integer inferredCount = 0;
	private static List<SimpleInstance> inferrableHumanEvents = new ArrayList<>();
	private SkipInstanceChecker skipInstanceChecker;

	public ReactionInferrer(
		ConfigProperties configProperties,
		@Qualifier("targetSpeciesCode") String speciesCode,
		OrthologousEntityGenerator orthologousEntityGenerator,
		InstanceUtilities instanceUtilities,
		ProteinCountUtility proteinCountUtility,
		StableIdentifierGenerator stableIdentifierGenerator,
		SkipInstanceChecker skipInstanceChecker,
		CuratorToolAPI curatorToolAPI) throws Exception {

		this.configProperties = configProperties;
		this.speciesCode = speciesCode;

		this.orthologousEntityGenerator = orthologousEntityGenerator;
		this.instanceUtilities = instanceUtilities;
		this.proteinCountUtility = proteinCountUtility;

		this.stableIdentifierGenerator = stableIdentifierGenerator;
		this.skipInstanceChecker = skipInstanceChecker;

		initializeEligibleAndInferredFiles();
		this.curatorToolAPI = curatorToolAPI;
	}

	// Infers PhysicalEntity instances of input, output, catalyst activity, and regulations that are associated with
	// incoming reactionInst.
	public void inferReaction(SimpleInstance reactionInst) throws Exception {
		if (skipInstanceChecker.instanceShouldBeSkipped(reactionInst) ||
			alreadyInferred(reactionInst) ||
			!hasProteins(reactionInst)) {
			return;
		}

		logger.info("Passed skip tests, RlE eligible for inference");

		SimpleInstance infReactionInst = createInitialInferredInstance(reactionInst);

		recordEligibleEvent(reactionInst);

		if (!inferInputsOutputsAndCatalysts(reactionInst, infReactionInst)) {
			return;
		}

		finalizeInferredReaction(reactionInst, infReactionInst);
	}

	public Map<Long, SimpleInstance> getInferredEvent() {
		return inferredEventMap;
	}

	public List<SimpleInstance> getInferrableHumanEvents() {
		return inferrableHumanEvents;
	}

	public int getEligibleCount() {
		return eligibleCount;
	}

	public int getInferredCount() {
		return inferredCount;
	}

	public void addAlreadyInferredEvents(SimpleInstance reactionInst, SimpleInstance previouslyInferredReactionInst) {
		inferredEventMap.put(reactionInst.getDbId(), previouslyInferredReactionInst);
		inferrableHumanEvents.add(reactionInst);
	}

	private boolean alreadyInferred(SimpleInstance reactionInst) {
		SimpleInstance inferredEvent = inferredEventMap.get(reactionInst.getDbId());
		if (inferredEvent != null) {
			logger.info("Already inferred RlE instance: " + inferredEvent);
			return true;
		}
		return false;
	}

	private SimpleInstance createInitialInferredInstance(SimpleInstance reactionInst) throws Exception {
		SimpleInstance infReactionInst = instanceUtilities.createNewInferredGKInstance(reactionInst);

		infReactionInst.setAttribute(ReactomeJavaConstants.name,
				reactionInst.getAttribute(ReactomeJavaConstants.name));
		infReactionInst.setAttribute(ReactomeJavaConstants.goBiologicalProcess,
				reactionInst.getAttribute(ReactomeJavaConstants.goBiologicalProcess));
		infReactionInst.setAttribute(ReactomeJavaConstants.summation,
				instanceUtilities.getSummationInstance());
		infReactionInst.setAttribute(ReactomeJavaConstants.evidenceType,
				instanceUtilities.getEvidenceType());
		infReactionInst.setDisplayName(reactionInst.getDisplayName());

		if (instanceUtilities.isValidAttribute(infReactionInst, ReactomeJavaConstants.releaseDate)) {
			infReactionInst.setAttribute(ReactomeJavaConstants.releaseDate, getDateOfRelease());
		}

		SimpleInstance orthoStableIdentifierInst = stableIdentifierGenerator
				.generateOrthologousStableId(infReactionInst, reactionInst);
		infReactionInst.setAttribute(ReactomeJavaConstants.stableIdentifier, orthoStableIdentifierInst);

		return infReactionInst;
	}

	private boolean hasProteins(SimpleInstance reactionInst) throws Exception {
		List<Integer> reactionProteinCounts = getProteinCountUtility().getDistinctProteinCounts(reactionInst);
		int reactionTotalProteinCounts = reactionProteinCounts.get(0);

		if (reactionTotalProteinCounts == 0) {
			logger.info("No distinct proteins found in instance -- terminating inference for " + reactionInst);
			return false;
		}

		logger.info("Total protein count for RlE: " + reactionTotalProteinCounts);
		return true;
	}

	private void recordEligibleEvent(SimpleInstance reactionInst) throws Exception {
		String eligibleEventName = reactionInst.getDbId() + "\t" + reactionInst.getDisplayName() + "\n";
		eligibleCount++;
		Files.write(getEligibleFile(), eligibleEventName.getBytes(),
				StandardOpenOption.APPEND, StandardOpenOption.CREATE);
	}

	private boolean inferInputsOutputsAndCatalysts(SimpleInstance reactionInst, SimpleInstance infReactionInst) throws Exception {
		logger.info("Inferring inputs...");
		if (!inferReactionInputsOrOutputs(reactionInst, infReactionInst, ReactomeJavaConstants.input)) {
			logger.info("Input inference unsuccessful -- terminating inference for " + reactionInst);
			return false;
		}

		logger.info("Inferring outputs...");
		if (!inferReactionInputsOrOutputs(reactionInst, infReactionInst, ReactomeJavaConstants.output)) {
			logger.info("Output inference unsuccessful -- terminating inference for " + reactionInst);
			return false;
		}

		logger.info("Inferring catalysts...");
		if (!inferReactionCatalysts(reactionInst, infReactionInst)) {
			logger.info("Catalyst inference unsuccessful -- terminating inference for " + reactionInst);
			return false;
		}

		return true;
	}

	private void finalizeInferredReaction(SimpleInstance reactionInst, SimpleInstance infReactionInst) throws Exception {
		logger.info("Inferring regulations...");
		List<SimpleInstance> inferredRegulations;
		try {
			inferredRegulations = inferReactionRegulations(reactionInst);
		} catch (MissingRegulatorForRequiredRegulationException e) {
			logger.info(e.getMessage());
			return;
		}

		// The inferred reaction is new, so everything it links to is set on it before it is committed, and it is
		// committed once. Each link added after a commit would cost a refresh and another commit.
		addRegulations(infReactionInst, inferredRegulations);
		addLinksToSourceReaction(infReactionInst, reactionInst);
		infReactionInst = storeInferredReaction(infReactionInst);

		// The source reaction's link back needs the inferred reaction's dbId, so it can only be added once the
		// inferred reaction is stored. Added rather than set, so the orthologous events already inferred for other
		// species are kept.
		getUtils().addOrthologousEventValue(reactionInst, infReactionInst);

		inferredEventMap.put(reactionInst.getDbId(), infReactionInst);
		recordInferredEvent(reactionInst, infReactionInst);
	}

	private SimpleInstance storeInferredReaction(SimpleInstance infReactionInst) throws Exception {
		return curatorToolAPI.commit(infReactionInst);
	}

	// Links the inferred reaction, which is not stored yet, to its source reaction through inferredFrom and
	// orthologousEvent.
	private void addLinksToSourceReaction(SimpleInstance infReactionInst, SimpleInstance reactionInst) {
		if (getUtils().isValidAttribute(infReactionInst, ReactomeJavaConstants.inferredFrom)) {
			getUtils().addAttributeValueIfNecessary(infReactionInst, reactionInst, ReactomeJavaConstants.inferredFrom);
		}
		getUtils().addAttributeValueIfNecessary(infReactionInst, reactionInst, ReactomeJavaConstants.orthologousEvent);
	}

	// Stores the inferred regulations (or finds identical stored ones) and adds them to the regulatedBy attribute of
	// the inferred reaction, which is not stored yet. A regulation refers only to its regulator, not to the reaction,
	// so it can be stored before the reaction is.
	private void addRegulations(SimpleInstance infReactionInst, List<SimpleInstance> inferredRegulations) throws Exception {
		if (inferredRegulations.isEmpty()) {
			return;
		}

		logger.info("Number of regulator(s) inferred: " + inferredRegulations.size());
		for (SimpleInstance infRegulation : inferredRegulations) {
			infRegulation = instanceUtilities.checkForIdenticalInstances(infRegulation, null);
			// Added rather than set, so each regulation is kept instead of replacing the one before it.
			getUtils().addAttributeValueIfNecessary(infReactionInst, infRegulation, ReactomeJavaConstants.regulatedBy);
		}
	}

	private void recordInferredEvent(SimpleInstance reactionInst, SimpleInstance infReactionInst) throws Exception {
		inferredCount++;
		inferrableHumanEvents.add(reactionInst);
		String inferredEventStr = infReactionInst.getDbId() + "\t" + infReactionInst.getDisplayName() + "\n";
		Files.write(getInferredFile(), inferredEventStr.getBytes(),
				StandardOpenOption.APPEND, StandardOpenOption.CREATE);
	}
	
	// Function used to create inferred PhysicalEntities contained in the 'input' or 'output' attributes of the current
	// reaction instance.
	@SuppressWarnings("unchecked")
	private boolean inferReactionInputsOrOutputs(
		SimpleInstance reactionInst, SimpleInstance infReactionInst, String attribute) throws Exception {
		Collection<SimpleInstance> attributeInstances = getUtils().getInflatedInstances(reactionInst, attribute);
		logger.info("Total " + attribute + " instances: " + attributeInstances.size());
		logger.info(capitalizeFirstLetter(attribute) + " instances: " + attributeInstances);

		List<SimpleInstance> infAttributeInstances = new ArrayList<>();
		for (SimpleInstance attributeInst : attributeInstances) {
			SimpleInstance infAttributeInst = getOrthologousEntityGenerator().createOrthoEntity(attributeInst, false);
			if (infAttributeInst == null) {
				return false; // Exit with failure if input/output not inferrable
			}
			infAttributeInstances.add(infAttributeInst);
		}
		infReactionInst.setAttribute(attribute, infAttributeInstances);
		logger.info("Completed " + attribute + " inference");
		return true;
	}
	
	// Function used to create inferred catalysts associated with the current reaction instance.
	// Infers all PhysicalEntity's associated with the reaction's 'catalystActivity' and 'activeUnit' attributes

	@SuppressWarnings("unchecked")
	private boolean inferReactionCatalysts(SimpleInstance reactionInst, SimpleInstance infReactionInst) throws Exception {
		Collection<SimpleInstance> catalystInstances = getUtils().getInflatedInstances(reactionInst, ReactomeJavaConstants.catalystActivity);
		logCatalystInfo(catalystInstances);

		for (SimpleInstance catalystInst : catalystInstances) {
			logger.info("Attempting catalyst inference: " + catalystInst);
			SimpleInstance infCatalystInst = getOrCreateInferredCatalyst(catalystInst);
			if (infCatalystInst == null) {
				return false;
			}
			infReactionInst.setAttribute(ReactomeJavaConstants.catalystActivity, infCatalystInst);
		}

		logger.info("Completed catalyst(s) inference for " + reactionInst);
		return true;
	}

	private void logCatalystInfo(Collection<SimpleInstance> catalystInstances) {
		logger.info("Total CatalystActivity instances: " + catalystInstances.size());
		if (!catalystInstances.isEmpty()) {
			logger.info("Catalyst instance(s): " + catalystInstances);
		}
	}
	private SimpleInstance getOrCreateInferredCatalyst(SimpleInstance sourceCatalyst) throws Exception {
		if (inferredCatalystMap.get(sourceCatalyst.getDbId()) != null) {
			logger.info("Inferred catalyst already exists");
			return inferredCatalystMap.get(sourceCatalyst.getDbId());
		}

		SimpleInstance inferredCatalyst = createInferredCatalyst(sourceCatalyst);
		inferredCatalystMap.put(sourceCatalyst.getDbId(), inferredCatalyst);
		return inferredCatalyst;
	}

	private SimpleInstance createInferredCatalyst(SimpleInstance sourceCatalyst) throws Exception {
		SimpleInstance inferredCatalyst = instanceUtilities.createNewInferredGKInstance(sourceCatalyst);

		inferActivityAndDisplayName(sourceCatalyst, inferredCatalyst);

		if (!inferPhysicalEntity(sourceCatalyst, inferredCatalyst)) {
			return null;
		}

		// Handle active units
		inferActiveUnits(sourceCatalyst, inferredCatalyst);

		// Finalize catalyst instance
		inferredCatalyst = finalizeCatalystInstance(inferredCatalyst);

		return inferredCatalyst;
	}

	private void inferActivityAndDisplayName(SimpleInstance sourceCatalyst, SimpleInstance inferredCatalyst) {
		inferredCatalyst.setAttribute(ReactomeJavaConstants.activity,
				sourceCatalyst.getAttribute(ReactomeJavaConstants.activity));
		inferredCatalyst.setDisplayName(sourceCatalyst.getDisplayName());
	}

	private boolean inferPhysicalEntity(SimpleInstance sourceCatalyst, SimpleInstance inferredCatalyst) throws Exception {
		SimpleInstance catalystPEInst = getUtils().getInflatedInstance(sourceCatalyst, ReactomeJavaConstants.physicalEntity);
		if (catalystPEInst != null) {
			logger.info("Catalyst PE instance: " + catalystPEInst);
			SimpleInstance infCatalystPEInst = getOrthologousEntityGenerator().createOrthoEntity(catalystPEInst, false);
			if (infCatalystPEInst != null) {
				inferredCatalyst.setAttribute(ReactomeJavaConstants.physicalEntity, infCatalystPEInst);
				return true;
			}
			return false;
		}
		return true;
	}

	@SuppressWarnings("unchecked")
	private void inferActiveUnits(SimpleInstance sourceCatalyst, SimpleInstance inferredCatalyst) throws Exception {
		List<SimpleInstance> activeUnits = new ArrayList<>();
		Collection<SimpleInstance> activeUnitInstances =
			getUtils().getInflatedInstances(sourceCatalyst, ReactomeJavaConstants.activeUnit);

		logger.info("Total active unit instances: " + activeUnitInstances);
		if (!activeUnitInstances.isEmpty()) {
			logger.info("Active unit instance(s): " + activeUnitInstances);
			for (SimpleInstance activeUnitInst : activeUnitInstances) {
				logger.info("Active Unit instance: " + activeUnitInst);
				SimpleInstance infActiveUnitInst = getOrthologousEntityGenerator().createOrthoEntity(activeUnitInst, false);
				if (infActiveUnitInst != null) {
					activeUnits.add(infActiveUnitInst);
				}
			}
		}
		inferredCatalyst.setAttribute(ReactomeJavaConstants.activeUnit, activeUnits);
	}

	private SimpleInstance finalizeCatalystInstance(SimpleInstance inferredCatalyst) throws Exception {
		inferredCatalyst = getUtils().checkForIdenticalInstances(inferredCatalyst, null);
		return inferredCatalyst;
	}

	// Function used to infer regulation instances. Logic existed for regulators that had CatalystActivity and
	// Event instances, but they have never come up in the many times this has been run.
	@SuppressWarnings("unchecked")
	private List<SimpleInstance> inferReactionRegulations(SimpleInstance reactionInst) throws Exception {
		List<SimpleInstance> inferredRegulations = new ArrayList<>();

		Collection<SimpleInstance> regulations = getRegulationsForReaction(reactionInst);
		if (regulations.isEmpty()) {
			logger.info("No regulations found for " + reactionInst);
			return inferredRegulations;
		}

		for (SimpleInstance regulation : regulations) {
			SimpleInstance inferredRegulation = inferSingleRegulation(regulation);
			if (inferredRegulation != null) {
				inferredRegulations.add(inferredRegulation);
			}
		}

		return inferredRegulations;
	}

	@SuppressWarnings("unchecked")
	private Collection<SimpleInstance> getRegulationsForReaction(SimpleInstance reactionInst) {
		return getUtils().getInflatedInstances(reactionInst, ReactomeJavaConstants.regulatedBy);
	}

	private SimpleInstance inferSingleRegulation(SimpleInstance regulation) throws Exception {
		SimpleInstance regulator = getUtils().getInflatedInstance(regulation, ReactomeJavaConstants.regulator);
		if (regulator == null) {
			logger.warn("Regulation instance has no regulator, skipping");
			return null;
		}

		SimpleInstance inferredRegulator = inferRegulator(regulator);
		if (inferredRegulator == null) {
			logger.info("Could not infer regulator for " + regulation);

			if (isRequirement(regulation)) {
				throw new MissingRegulatorForRequiredRegulationException(
					"No regulator inferrable for required regulation - " + regulation
				);
			}
			return null;
		}

		return createInferredRegulation(regulation, inferredRegulator);
	}

	private SimpleInstance inferRegulator(SimpleInstance regulator) throws Exception {
		if (getUtils().isSchemaClassType(regulator, ReactomeJavaConstants.PhysicalEntity)) {
			return getOrthologousEntityGenerator().createOrthoEntity(regulator, false);
		} else {
			String errorMessage = regulator.getSchemaClassName() + " regulators are not currently supported";
			logger.error(errorMessage);
			throw new RuntimeException(errorMessage);
		}
	}

	private boolean isRequirement(SimpleInstance regulation) {
		return getUtils().isSchemaClassType(regulation, ReactomeJavaConstants.Requirement);
	}

	private SimpleInstance createInferredRegulation(SimpleInstance sourceRegulation, SimpleInstance inferredRegulator)
			throws Exception {
		SimpleInstance inferredRegulation = getUtils().createNewInferredGKInstance(sourceRegulation);
		inferredRegulation.setAttribute(ReactomeJavaConstants.regulator, inferredRegulator);
		inferredRegulation.setDisplayName(sourceRegulation.getDisplayName());

		return inferredRegulation;
	}

	private void initializeEligibleAndInferredFiles() throws IOException {
		Files.deleteIfExists(getEligibleFile());
		Files.deleteIfExists(getInferredFile());
	}

	private Path getEligibleFile() {
		return Paths.get("eligible_" + getSpeciesCode() + "_75.txt");
	}

	private Path getInferredFile() {
		return Paths.get("inferred_" + getSpeciesCode() + "_75.txt");
	}

	private String getSpeciesCode() {
		return this.speciesCode;
	}

	private OrthologousEntityGenerator getOrthologousEntityGenerator() {
		return this.orthologousEntityGenerator;
	}

	private InstanceUtilities getUtils() {
		return this.instanceUtilities;
	}

	private ConfigProperties getConfigProperties() {
		return this.configProperties;
	}

	private String getDateOfRelease() {
		return getConfigProperties().getDateOfRelease();
	}

	private ProteinCountUtility getProteinCountUtility() {
		return this.proteinCountUtility;
	}

	private String capitalizeFirstLetter(String str) {
		return str.substring(0, 1).toUpperCase() + str.substring(1);
	}

	private static class MissingRegulatorForRequiredRegulationException extends Exception {
		public MissingRegulatorForRequiredRegulationException(String message) {
			super(message);
		}
	}
}
