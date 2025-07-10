package org.reactome.orthoinference;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.sql.SQLException;
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
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

@Component
public class ReactionInferrer {

	private static final Logger logger = LogManager.getLogger();

	private ConfigProperties configProperties;
	private String speciesCode;
	private OrthologousEntityGenerator orthologousEntityGenerator;
	private StableIdentifierGenerator stableIdentifierGenerator;
	private InstanceUtilities instanceUtilities;

	private static Map<GKInstance, GKInstance> inferredCatalystMap = new HashMap<>();
	private static Map<GKInstance, GKInstance> inferredEventMap = new HashMap<>();
	private static Integer eligibleCount = 0;
	private static Integer inferredCount = 0;
	private static List<GKInstance> inferrableHumanEvents = new ArrayList<>();
	private SkipInstanceChecker skipInstanceChecker;

	public ReactionInferrer(
		ConfigProperties configProperties,
		@Qualifier("targetSpeciesCode") String speciesCode,
		OrthologousEntityGenerator orthologousEntityGenerator,
		InstanceUtilities instanceUtilities,
		StableIdentifierGenerator stableIdentifierGenerator,
		SkipInstanceChecker skipInstanceChecker
	) throws Exception {

		this.configProperties = configProperties;
		this.speciesCode = speciesCode;

		this.orthologousEntityGenerator = orthologousEntityGenerator;
		this.instanceUtilities = instanceUtilities;

		this.stableIdentifierGenerator = stableIdentifierGenerator;
		this.skipInstanceChecker = skipInstanceChecker;

		initializeEligibleAndInferredFiles();
	}

	// Infers PhysicalEntity instances of input, output, catalyst activity, and regulations that are associated with
	// incoming reactionInst.
	public void inferReaction(GKInstance reactionInst) throws Exception {
		if (skipInstanceChecker.instanceShouldBeSkipped(reactionInst) ||
			alreadyInferred(reactionInst) ||
			!hasProteins(reactionInst)) {
			return;
		}

		logger.info("Passed skip tests, RlE eligible for inference");

		GKInstance infReactionInst = createInitialInferredInstance(reactionInst);

		recordEligibleEvent(reactionInst);

		if (!inferInputsOutputsAndCatalysts(reactionInst, infReactionInst)) {
			return;
		}

		finalizeInferredReaction(reactionInst, infReactionInst);
	}

	public Map<GKInstance, GKInstance> getInferredEvent() {
		return inferredEventMap;
	}

	public List<GKInstance> getInferrableHumanEvents() {
		return inferrableHumanEvents;
	}

	public int getEligibleCount() {
		return eligibleCount;
	}

	public int getInferredCount() {
		return inferredCount;
	}

	public void addAlreadyInferredEvents(GKInstance reactionInst, GKInstance previouslyInferredReactionInst) {
		inferredEventMap.put(reactionInst, previouslyInferredReactionInst);
		inferrableHumanEvents.add(reactionInst);
	}

	private boolean alreadyInferred(GKInstance reactionInst) throws Exception {
		GKInstance inferredEvent = inferredEventMap.get(reactionInst);
		if (inferredEvent != null) {
			logger.info("Already inferred RlE instance: " + inferredEvent);
			return true;
		}
		return false;
	}

	private GKInstance createInitialInferredInstance(GKInstance reactionInst) throws Exception {
		GKInstance infReactionInst = instanceUtilities.createNewInferredGKInstance(reactionInst);

		infReactionInst.addAttributeValue(ReactomeJavaConstants.name,
				reactionInst.getAttributeValuesList(ReactomeJavaConstants.name));
		infReactionInst.addAttributeValue(ReactomeJavaConstants.goBiologicalProcess,
				reactionInst.getAttributeValue(ReactomeJavaConstants.goBiologicalProcess));
		infReactionInst.addAttributeValue(ReactomeJavaConstants.summation,
				instanceUtilities.getSummationInstance());
		infReactionInst.addAttributeValue(ReactomeJavaConstants.evidenceType,
				instanceUtilities.getEvidenceType());
		infReactionInst.addAttributeValue(ReactomeJavaConstants._displayName,
				reactionInst.getAttributeValue(ReactomeJavaConstants._displayName));

		if (infReactionInst.getSchemClass().isValidAttribute(ReactomeJavaConstants.releaseDate)) {
			infReactionInst.addAttributeValue(ReactomeJavaConstants.releaseDate, getDateOfRelease());
		}

		GKInstance orthoStableIdentifierInst = stableIdentifierGenerator
				.generateOrthologousStableId(infReactionInst, reactionInst);
		infReactionInst.addAttributeValue(ReactomeJavaConstants.stableIdentifier, orthoStableIdentifierInst);

		return infReactionInst;
	}

	private boolean hasProteins(GKInstance reactionInst) throws Exception {
		List<Integer> reactionProteinCounts = getProteinCountUtility().getDistinctProteinCounts(reactionInst);
		int reactionTotalProteinCounts = reactionProteinCounts.get(0);

		if (reactionTotalProteinCounts == 0) {
			logger.info("No distinct proteins found in instance -- terminating inference for " + reactionInst);
			return false;
		}

		logger.info("Total protein count for RlE: " + reactionTotalProteinCounts);
		return true;
	}

	private void recordEligibleEvent(GKInstance reactionInst) throws Exception {
		String eligibleEventName = reactionInst.getAttributeValue(ReactomeJavaConstants.DB_ID).toString() + "\t" +
				reactionInst.getDisplayName() + "\n";
		eligibleCount++;
		Files.write(getEligibleFile(), eligibleEventName.getBytes(),
				StandardOpenOption.APPEND, StandardOpenOption.CREATE);
	}

	private boolean inferInputsOutputsAndCatalysts(GKInstance reactionInst, GKInstance infReactionInst) throws Exception {
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

	private void finalizeInferredReaction(GKInstance reactionInst, GKInstance infReactionInst) throws Exception {
		logger.info("Inferring regulations...");
		List<GKInstance> inferredRegulations;
		try {
			inferredRegulations = inferReactionRegulations(reactionInst);
		} catch (MissingRegulatorForRequiredRegulationException e) {
			logger.info(e.getMessage());
			return;
		}

		storeInferredReaction(infReactionInst);
		updateInferredFromAndOrthologousEventAttributes(reactionInst, infReactionInst);

		inferredEventMap.put(reactionInst, infReactionInst);

		processInferredRegulations(infReactionInst, inferredRegulations);
		recordInferredEvent(reactionInst, infReactionInst);
	}

	private void storeInferredReaction(GKInstance infReactionInst) throws Exception {
		getCurrentDBA().storeInstance(infReactionInst);
	}

	private void updateInferredFromAndOrthologousEventAttributes(GKInstance reactionInst, GKInstance infReactionInst) throws Exception {
		if (infReactionInst.getSchemClass().isValidAttribute(ReactomeJavaConstants.inferredFrom)) {
			infReactionInst = instanceUtilities.addAttributeValueIfNecessary(
					infReactionInst, reactionInst, ReactomeJavaConstants.inferredFrom);
			getCurrentDBA().updateInstanceAttribute(infReactionInst, ReactomeJavaConstants.inferredFrom);
		}

		infReactionInst = instanceUtilities.addAttributeValueIfNecessary(
				infReactionInst, reactionInst, ReactomeJavaConstants.orthologousEvent);
		getCurrentDBA().updateInstanceAttribute(infReactionInst, ReactomeJavaConstants.orthologousEvent);

		reactionInst.addAttributeValue(ReactomeJavaConstants.orthologousEvent, infReactionInst);
		getCurrentDBA().updateInstanceAttribute(reactionInst, ReactomeJavaConstants.orthologousEvent);

	}

	private void processInferredRegulations(GKInstance infReactionInst, List<GKInstance> inferredRegulations) throws Exception {
		if (!inferredRegulations.isEmpty()) {
			logger.info("Number of regulator(s) inferred: " + inferredRegulations.size());
			for (GKInstance infRegulation : inferredRegulations) {
				infRegulation = instanceUtilities.checkForIdenticalInstances(infRegulation, null);
				infReactionInst.addAttributeValue("regulatedBy", infRegulation);
				getCurrentDBA().updateInstanceAttribute(infReactionInst, "regulatedBy");
			}
		}
	}

	private void recordInferredEvent(GKInstance reactionInst, GKInstance infReactionInst) throws Exception {
		inferredCount++;
		inferrableHumanEvents.add(reactionInst);
		String inferredEventStr = infReactionInst.getAttributeValue(ReactomeJavaConstants.DB_ID).toString() + "\t" +
				infReactionInst.getDisplayName() + "\n";
		Files.write(getInferredFile(), inferredEventStr.getBytes(),
				StandardOpenOption.APPEND, StandardOpenOption.CREATE);
	}
	
	// Function used to create inferred PhysicalEntities contained in the 'input' or 'output' attributes of the current
	// reaction instance.
	@SuppressWarnings("unchecked")
	private boolean inferReactionInputsOrOutputs(
		GKInstance reactionInst, GKInstance infReactionInst, String attribute) throws Exception {
		Collection<GKInstance> attributeInstances =
			(Collection<GKInstance>) reactionInst.getAttributeValuesList(attribute);
		logger.info("Total " + attribute + " instances: " + attributeInstances.size());
		logger.info(capitalizeFirstLetter(attribute) + " instances: " + attributeInstances);

		List<GKInstance> infAttributeInstances = new ArrayList<>();
		for (GKInstance attributeInst : attributeInstances) {
			GKInstance infAttributeInst = getOrthologousEntityGenerator().createOrthoEntity(attributeInst, false);
			if (infAttributeInst == null) {
				return false; // Exit with failure if input/output not inferrable
			}
			infAttributeInstances.add(infAttributeInst);
		}
		infReactionInst.addAttributeValue(attribute, infAttributeInstances);
		logger.info("Completed " + attribute + " inference");
		return true;
	}
	
	// Function used to create inferred catalysts associated with the current reaction instance.
	// Infers all PhysicalEntity's associated with the reaction's 'catalystActivity' and 'activeUnit' attributes
	@SuppressWarnings("unchecked")
	private boolean inferReactionCatalysts(GKInstance reactionInst, GKInstance infReactionInst)
		throws Exception {
		Collection<GKInstance> catalystInstances =
			(Collection<GKInstance>) reactionInst.getAttributeValuesList(ReactomeJavaConstants.catalystActivity);
		logger.info("Total CatalystActivity instances: " + catalystInstances.size());
		if (!catalystInstances.isEmpty()) {
			logger.info("Catalyst instance(s): " + catalystInstances);
		}
		for (GKInstance catalystInst : catalystInstances) {
			logger.info("Attempting catalyst inference: " + catalystInst);
			if (inferredCatalystMap.get(catalystInst) == null) {
				GKInstance infCatalystInst = instanceUtilities.createNewInferredGKInstance(catalystInst);
				infCatalystInst.setDbAdaptor(getCurrentDBA());
				infCatalystInst.addAttributeValue(ReactomeJavaConstants.activity, catalystInst.getAttributeValue(ReactomeJavaConstants.activity));
				GKInstance catalystPEInst = (GKInstance) catalystInst.getAttributeValue(ReactomeJavaConstants.physicalEntity);
				if (catalystPEInst != null) {
					logger.info("Catalyst PE instance: " + catalystPEInst);
					GKInstance infCatalystPEInst = getOrthologousEntityGenerator().createOrthoEntity(
						catalystPEInst, false);
					if (infCatalystPEInst != null) {
						infCatalystInst.addAttributeValue(ReactomeJavaConstants.physicalEntity, infCatalystPEInst);
					} else {
						return false;
					}
				}

				List<GKInstance> activeUnits = new ArrayList<>();
				Collection<GKInstance> activeUnitInstances =
					(Collection<GKInstance>) catalystInst.getAttributeValuesList(ReactomeJavaConstants.activeUnit);
				logger.info("Total active unit instances: " + activeUnitInstances);
				if (!activeUnitInstances.isEmpty()) {
					logger.info("Active unit instance(s): " + activeUnitInstances);
					for (GKInstance activeUnitInst : activeUnitInstances) {
						logger.info("Active Unit instance: " + activeUnitInst);
						GKInstance infActiveUnitInst = getOrthologousEntityGenerator().createOrthoEntity(
							activeUnitInst, false);
						if (infActiveUnitInst != null) {
							activeUnits.add(infActiveUnitInst);
						}
					}
				}
				infCatalystInst.addAttributeValue(ReactomeJavaConstants.activeUnit, activeUnits);
				infCatalystInst.addAttributeValue(ReactomeJavaConstants._displayName, catalystInst.getAttributeValue(ReactomeJavaConstants._displayName));
				infCatalystInst = instanceUtilities.checkForIdenticalInstances(infCatalystInst, null);
				inferredCatalystMap.put(catalystInst, infCatalystInst);
			} else {
				logger.info("Inferred catalyst already exists");
			}
			infReactionInst.addAttributeValue(ReactomeJavaConstants.catalystActivity, inferredCatalystMap.get(catalystInst));
		}
		logger.info("Completed catalyst inference");
		return true;
	}
	

	// Function used to infer regulation instances. Logic existed for regulators that had CatalystActivity and
	// Event instances, but they have never come up in the many times this has been run.
	@SuppressWarnings("unchecked")
	private List<GKInstance> inferReactionRegulations(GKInstance reactionInst) throws Exception {
		List<GKInstance> inferredRegulations = new ArrayList<>();

		Collection<GKInstance> regulations = getRegulationsForReaction(reactionInst);
		if (regulations.isEmpty()) {
			logger.info("No regulations found for " + reactionInst);
			return inferredRegulations;
		}

		for (GKInstance regulation : regulations) {
			GKInstance inferredRegulation = inferSingleRegulation(regulation);
			if (inferredRegulation != null) {
				inferredRegulations.add(inferredRegulation);
			}
		}

		return inferredRegulations;
	}

	@SuppressWarnings("unchecked")
	private Collection<GKInstance> getRegulationsForReaction(GKInstance reactionInst) throws Exception {
		return reactionInst.getAttributeValuesList(ReactomeJavaConstants.regulatedBy);
	}

	private GKInstance inferSingleRegulation(GKInstance regulation) throws Exception {
		GKInstance regulator = (GKInstance) regulation.getAttributeValue(ReactomeJavaConstants.regulator);
		if (regulator == null) {
			logger.warn("Regulation instance has no regulator, skipping");
			return null;
		}

		GKInstance inferredRegulator = inferRegulator(regulator);
		if (inferredRegulator == null) {
			if (isRequirement(regulation)) {
				throw new MissingRegulatorForRequiredRegulationException(
					"No regulator inferrable for required regulation - " + regulation
				);
			}
			logger.info("Could not infer regulator for " + regulation);
			return null;
		}

		return createInferredRegulation(regulation, inferredRegulator);
	}

	private GKInstance inferRegulator(GKInstance regulator) throws Exception {
		if (regulator.getSchemClass().isa(ReactomeJavaConstants.PhysicalEntity)) {
			return orthologousEntityGenerator.createOrthoEntity(regulator, false);
		} else {
			String errorMessage = regulator.getSchemClass().getName() + " regulators are not currently supported";
			logger.error(errorMessage);
			throw new RuntimeException(errorMessage);
		}
	}

	private boolean isRequirement(GKInstance regulation) {
		return regulation.getSchemClass().isa(ReactomeJavaConstants.Requirement);
	}

	private GKInstance createInferredRegulation(GKInstance sourceRegulation, GKInstance inferredRegulator)
			throws Exception {
		GKInstance inferredRegulation = instanceUtilities.createNewInferredGKInstance(sourceRegulation);
		inferredRegulation.setDbAdaptor(getCurrentDBA());
		inferredRegulation.addAttributeValue(ReactomeJavaConstants.regulator, inferredRegulator);
		inferredRegulation.addAttributeValue(ReactomeJavaConstants._displayName, sourceRegulation.getAttributeValue(ReactomeJavaConstants._displayName));

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

	private MySQLAdaptor getCurrentDBA() throws SQLException {
		return getConfigProperties().getCurrentDBA();
	}

	private String getDateOfRelease() {
		return getConfigProperties().getDateOfRelease();
	}

	private ProteinCountUtility getProteinCountUtility() {
		return getUtils().getProteinCountUtility();
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
