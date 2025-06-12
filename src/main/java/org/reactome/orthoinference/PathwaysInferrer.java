package org.reactome.orthoinference;

import java.sql.SQLException;
import java.util.*;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;

import static org.reactome.util.general.CollectionUtils.safeList;

import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.springframework.stereotype.Component;

@Component
public class PathwaysInferrer {

	private static final Logger logger = LogManager.getLogger();

	private ConfigProperties configProperties;
	private InstanceUtilities instanceUtilities;
	private Utils utils;

	private static List<GKInstance> updatedInferrableHumanEvents = new ArrayList<>();
	private static Map<GKInstance, GKInstance> sourceInstanceToInferredInstance = new HashMap<>();
	//private static GKInstance diseasePathwayInst;

	public PathwaysInferrer(ConfigProperties configProperties, InstanceUtilities instanceUtilities, Utils utils) {
		this.configProperties = configProperties;
		this.instanceUtilities = instanceUtilities;
		this.utils = utils;
	}

	// This class populates species pathways with the instances that have been inferred.
	// This was copied heavily from the Perl, so my explanations are a little sparse here.
	public void inferPathways(List<GKInstance> inferrableHumanReactionLikeEvents) throws Exception {
		//diseasePathwayInst = getDiseaseInstance();

		logger.info("Beginning Pathway inference");
		updatedInferrableHumanEvents.addAll(inferrableHumanReactionLikeEvents);

		// First, go through each of the inferred RlE instances and generate the entire pathway hierarchy it is
		// associated with. Inferred Reactions are not added to the Pathway at this point. This includes the immediate
		// Pathway, but also all parent pathways up to its TopLevelPathway.
		logger.info("Building inferred Pathway hierarchies");
		Set<Long> seenPathwayHierarchy = new HashSet<>();
		for (GKInstance inferrableHumanReactionLikeEvent : inferrableHumanReactionLikeEvents) {
			logger.info("Building inferred pathways for RlE: " + inferrableHumanReactionLikeEvent);
			if (!seenPathwayHierarchy.contains(inferrableHumanReactionLikeEvent.getDBID())) {
				createInferredPathwayHierarchy(inferrableHumanReactionLikeEvent);
				seenPathwayHierarchy.add(inferrableHumanReactionLikeEvent.getDBID());
			} else {
				logger.info("Inferred pathways already exist for RlE");
			}
		}
		logger.info("Finished building inferred Pathway hierarchies");

		// After generating the inferred Pathways hierarchys, the associated inferred Events (RlEs & Pathways) need
		// to be added to them.
		logger.info("Populating inferred Pathways with inferred Events");
		addInferredEventsToInferredPathways();
		logger.info("Finished populating inferred Pathways with inferred Events");

		//TODO: LOG starting HERE

		// Connect preceding events to RlEs, if they have any in the source species.
		logger.info("Adding preceding events to inferred Events");
		inferPrecedingEvents();
		logger.info("Finished adding preceding events to inferred Events");

		// Any source species Events (Pathway or RlEs) that were modified during Pathway inference are updated with a
		// 'modified' instance edit.
		updateModifiedAttributeIfNecessary();
	}

	@SuppressWarnings("unchecked")
	// This generates the inferred Pathway of an inferred RlE. It iterates, inferring parent Pathways until reaching
	// the TopLevelPathway.  Inferred Reactions are not added to the Pathway at this step.
	private void createInferredPathwayHierarchy(GKInstance humanEvent) throws Exception {
		List<GKInstance> humanPathwayReferralInstances =
			safeList(humanEvent.getReferers(ReactomeJavaConstants.hasEvent));

		if (humanPathwayReferralInstances.isEmpty()) {
			logger.info("Top Level Pathway inferred: " + humanEvent);
			return;
		}

		for (GKInstance humanPathwayReferralInstance : humanPathwayReferralInstances) {
			logger.info("Generating inferred Pathway: " + humanPathwayReferralInstance);
			// Pathways that have been inferred already are skipped, as are Pathways that are only children of the
			// Disease TopLevelPathway.
			if (hasNotBeenInferred(humanPathwayReferralInstance) &&
				!instanceUtilities.onlyInDiseasePathway(humanPathwayReferralInstance)) {

				inferPathway(humanPathwayReferralInstance);
			}
			createInferredPathwayHierarchy(humanPathwayReferralInstance);
		}
	}

	private boolean hasNotBeenInferred(GKInstance sourcePathwayReferralInst) {
		return sourceInstanceToInferredInstance.get(sourcePathwayReferralInst) == null;
	}

	private void inferPathway(GKInstance humanPathway) throws Exception {
		GKInstance inferredPathway = instanceUtilities.createNewInferredGKInstance(humanPathway);
		inferredPathway.addAttributeValue(
			ReactomeJavaConstants.name, humanPathway.getAttributeValuesList(ReactomeJavaConstants.name));
		inferredPathway.addAttributeValue(ReactomeJavaConstants.summation, getUtils().getSummationInstance());
		if (inferredPathway.getSchemClass().isValidAttribute(ReactomeJavaConstants.releaseDate)) {
			inferredPathway.addAttributeValue(ReactomeJavaConstants.releaseDate, getConfigProperties().getDateOfRelease());
		}
		inferredPathway.addAttributeValue(ReactomeJavaConstants.inferredFrom, humanPathway);
		inferredPathway.addAttributeValue(ReactomeJavaConstants.evidenceType, getUtils().getEvidenceType());
		for (GKInstance goBioProcessInst : getGoBiologicalProcess(humanPathway)) {
			inferredPathway.addAttributeValue(ReactomeJavaConstants.goBiologicalProcess, goBioProcessInst);
		}
		inferredPathway.addAttributeValue(ReactomeJavaConstants.orthologousEvent, humanPathway);

		if (humanPathway.getSchemClass().isa(ReactomeJavaConstants.ReactionlikeEvent)) {
			logger.warn(humanPathway + " is a ReactionLikeEvent, which is unexpected -- refer to infer_events.pl");
		}
		inferredPathway.setDisplayName(humanPathway.getDisplayName());
		sourceInstanceToInferredInstance.put(humanPathway, inferredPathway);
		GKInstance orthoStableIdentifierInst = getUtils().getStableIdentifierGenerator()
			.generateOrthologousStableId(inferredPathway, humanPathway);
		inferredPathway.addAttributeValue(ReactomeJavaConstants.stableIdentifier, orthoStableIdentifierInst);
		getCurrentDBA().storeInstance(inferredPathway);

		// This was replaced with addAttributeValueIfNecessary due to a bug where a Pathway instance's
		// 'OrthologousEvent' attribute was being replaced, instead of being added to the existing array when the
		// script was executed from a jar (rather than from Eclipse) (Justin Cook 2018)
		humanPathway = instanceUtilities.addAttributeValueIfNecessary(
			humanPathway, inferredPathway, ReactomeJavaConstants.orthologousEvent);
		getCurrentDBA().updateInstanceAttribute(humanPathway, ReactomeJavaConstants.orthologousEvent);

		//TODO: At this point, sourcePathwayReferralInst is always a Pathway.
		// Perhaps move to its own data structure? Holdout from Perl...
		updatedInferrableHumanEvents.add(humanPathway);
	}

	// This populates the hasEvent slot of all inferred Pathways that were just generated with corresponding inferred
	// reactions
	private void addInferredEventsToInferredPathways() throws Exception {
		Set<Long> seenInferredPathway = new HashSet<>();
		for (GKInstance humanPathwayInst : getPathways(updatedInferrableHumanEvents)) {
			if (seenInferredPathway.contains(humanPathwayInst.getDBID())) {
				logger.info("Inferred Pathway has already been populated with inferred Events");
				continue;
			}

			GKInstance inferredPathwayInst = sourceInstanceToInferredInstance.get(humanPathwayInst);
			List<GKInstance> inferredEventInstances = getInferredEventInstances(humanPathwayInst);
			if (!isPathway(inferredPathwayInst)) {
				logger.info(humanPathwayInst + " and " +
					sourceInstanceToInferredInstance.get(humanPathwayInst) + " have different classes " +
					"(likely connected via manual inference");
				continue;
			}

			// Add inferred Events to inferred Pathway
			logger.info("Adding " + inferredEventInstances.size() + " inferred Event(s) to " +
				"inferred Pathway: " + sourceInstanceToInferredInstance.get(humanPathwayInst));
			for (GKInstance infEventInst : inferredEventInstances) {
				inferredPathwayInst = instanceUtilities.addAttributeValueIfNecessary(
					inferredPathwayInst, infEventInst, ReactomeJavaConstants.hasEvent);
				sourceInstanceToInferredInstance.remove(humanPathwayInst);
				sourceInstanceToInferredInstance.put(humanPathwayInst, inferredPathwayInst);
			}
			getCurrentDBA().updateInstanceAttribute(
				sourceInstanceToInferredInstance.get(humanPathwayInst), ReactomeJavaConstants.hasEvent);

			seenInferredPathway.add(humanPathwayInst.getDBID());
		}
	}

	// Collect inferred Events associated with source Event
	private List<GKInstance> getInferredEventInstances(GKInstance humanPathwayInst) throws Exception {
		List<GKInstance> inferredEventInstances = new ArrayList<>();
		for (GKInstance eventInst :
			(Collection<GKInstance>) humanPathwayInst.getAttributeValuesList(ReactomeJavaConstants.hasEvent)) {
			if (sourceInstanceToInferredInstance.get(eventInst) != null) {
				inferredEventInstances.add(sourceInstanceToInferredInstance.get(eventInst));
			}
		}
		return inferredEventInstances;
	}


	@SuppressWarnings("unchecked")
	private void inferPrecedingEvents() throws Exception {
		Set<GKInstance> seenPrecedingEvent = new HashSet<>();
		for (GKInstance inferrableEventInst : updatedInferrableHumanEvents) {
			if (!seenPrecedingEvent.contains(inferrableEventInst)) {
				if (inferrableEventInst.getAttributeValue(ReactomeJavaConstants.precedingEvent)!= null) {
					logger.info("Adding preceding event to " + inferrableEventInst);
					List<GKInstance> precedingEventInstances = new ArrayList<>();
					// Find all preceding events for source instance that have an inferred counterpart
					for (GKInstance precedingEventInst :
						(Collection<GKInstance>) inferrableEventInst.getAttributeValuesList(ReactomeJavaConstants.precedingEvent)) {
						if (sourceInstanceToInferredInstance.get(precedingEventInst) != null) {
							precedingEventInstances.add(sourceInstanceToInferredInstance.get(precedingEventInst));
						}
					}
					Set<String> inferredPrecedingEvents = new HashSet<>();
					// Find any inferred preceding events that already exist for the inferred instance
					// (don't want to add any redundant preceding events)
					for (GKInstance precedingEventInst : (Collection<GKInstance>)
						sourceInstanceToInferredInstance.get(inferrableEventInst)
							.getAttributeValuesList(ReactomeJavaConstants.precedingEvent)) {
						inferredPrecedingEvents.add(precedingEventInst.getDBID().toString());
					}
					List<GKInstance> updatedPrecedingEventInstances = new ArrayList<>();
					// Find existing preceding events that haven't already been attached to the inferred instance
					for (GKInstance precedingEventInst : precedingEventInstances) {
						if (!inferredPrecedingEvents.contains(precedingEventInst.getDBID().toString())) {
							updatedPrecedingEventInstances.add(precedingEventInst);
						}
					}
					// Add preceding event to inferred instance
					if (!updatedPrecedingEventInstances.isEmpty()) {
						sourceInstanceToInferredInstance.get(inferrableEventInst).addAttributeValue(
							ReactomeJavaConstants.precedingEvent, updatedPrecedingEventInstances);
						getCurrentDBA().updateInstanceAttribute(sourceInstanceToInferredInstance.get(inferrableEventInst),
							ReactomeJavaConstants.precedingEvent);
					}
				}
				seenPrecedingEvent.add(inferrableEventInst);
			}
		}
	}

	private void updateModifiedAttributeIfNecessary() throws Exception {

		Set<Long> seenInstanceEditInst = new HashSet<>();
		for (GKInstance humanPathwayInst : updatedInferrableHumanEvents) {
			if (!seenInstanceEditInst.contains(humanPathwayInst.getDBID())) {
				GKInstance createdInst = (GKInstance) humanPathwayInst.getAttributeValue(ReactomeJavaConstants.created);
				if (createdInst == null ||
					!createdInst.getDBID().toString().matches(getInstanceEdit().getDBID().toString())) {

					boolean modifiedExists = false;
					for (GKInstance modifiedInst :
						(Collection<GKInstance>) humanPathwayInst.getAttributeValuesList(ReactomeJavaConstants.modified)) {
						if (modifiedInst.getDBID().toString().matches(getInstanceEdit().getDBID().toString())) {
							modifiedExists = true;
						}
					}
					if (!modifiedExists) {
						humanPathwayInst.addAttributeValue(ReactomeJavaConstants.modified, getInstanceEdit());
						getCurrentDBA().updateInstanceAttribute(humanPathwayInst, ReactomeJavaConstants.modified);
					}
					seenInstanceEditInst.add(humanPathwayInst.getDBID());
				}
			}
		}
	}

	private Collection<GKInstance> getGoBiologicalProcess(GKInstance humanPathway) throws Exception {
		return (Collection<GKInstance>) humanPathway.getAttributeValuesList(ReactomeJavaConstants.goBiologicalProcess);
	}

	private List<GKInstance> getPathways(List<GKInstance> events) {
		return events.stream().filter(this::isPathway).collect(Collectors.toList());
	}

	private boolean isPathway(GKInstance event) {
		return event.getSchemClass().isValidAttribute(ReactomeJavaConstants.hasEvent);
	}

	public void setInferredEvent(Map<GKInstance,GKInstance> inferredEventCopy) {
		sourceInstanceToInferredInstance = inferredEventCopy;
	}

	private Utils getUtils() {
		return this.utils;
	}

	private ConfigProperties getConfigProperties() {
		return this.configProperties;
	}

	private MySQLAdaptor getCurrentDBA() throws SQLException {
		return getConfigProperties().getCurrentDBA();
	}

	private GKInstance getInstanceEdit() throws Exception {
		return this.instanceUtilities.getInstanceEdit();
	}

//	private static GKInstance getDiseaseInstance() throws Exception {
//		if (diseasePathwayInst == null) {
//			diseasePathwayInst = dba.fetchInstance(InstanceUtilities.getDiseasePathwayDbId());
//		}
//
//		return diseasePathwayInst;
//	}
}
