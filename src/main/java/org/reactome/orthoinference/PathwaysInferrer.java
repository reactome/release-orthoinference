package org.reactome.orthoinference;

import java.util.*;
import java.util.stream.Collectors;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import static org.reactome.util.general.CollectionUtils.safeList;

import org.gk.model.ReactomeJavaConstants;
import org.reactome.curation.model.SimpleInstance;
import org.springframework.stereotype.Component;

@Component
public class PathwaysInferrer {

	private static final Logger logger = LogManager.getLogger();

	private ConfigProperties configProperties;
	private StableIdentifierGenerator stableIdentifierGenerator;
	private InstanceUtilities instanceUtilities;
	private CuratorToolAPI curatorToolAPI;

	private static List<SimpleInstance> updatedInferrableHumanEvents = new ArrayList<>();
	private static Map<SimpleInstance, SimpleInstance> sourceInstanceToInferredInstance = new HashMap<>();

	public PathwaysInferrer(
		ConfigProperties configProperties,
		StableIdentifierGenerator stableIdentifierGenerator,
		InstanceUtilities instanceUtilities,
		CuratorToolAPI curatorToolAPI
	) {
		this.configProperties = configProperties;
		this.stableIdentifierGenerator = stableIdentifierGenerator;
		this.instanceUtilities = instanceUtilities;
		this.curatorToolAPI = curatorToolAPI;
	}

	// This class populates species pathways with the instances that have been inferred.
	// This was copied heavily from the Perl, so my explanations are a little sparse here.
	public void inferPathways(List<SimpleInstance> inferrableHumanReactionLikeEvents) throws Exception {
		//diseasePathwayInst = getDiseaseInstance();

		logger.info("Beginning Pathway inference");
		updatedInferrableHumanEvents.addAll(inferrableHumanReactionLikeEvents);

		// First, go through each of the inferred RlE instances and generate the entire pathway hierarchy it is
		// associated with. Inferred Reactions are not added to the Pathway at this point. This includes the immediate
		// Pathway, but also all parent pathways up to its TopLevelPathway.
		logger.info("Building inferred Pathway hierarchies");
		Set<Long> seenPathwayHierarchy = new HashSet<>();
		for (SimpleInstance inferrableHumanReactionLikeEvent : inferrableHumanReactionLikeEvents) {
			logger.info("Building inferred pathways for RlE: " + inferrableHumanReactionLikeEvent);
			if (!seenPathwayHierarchy.contains(inferrableHumanReactionLikeEvent.getDbId())) {
				createInferredPathwayHierarchy(inferrableHumanReactionLikeEvent);
				seenPathwayHierarchy.add(inferrableHumanReactionLikeEvent.getDbId());
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
		//updateModifiedAttributeIfNecessary();
	}

	@SuppressWarnings("unchecked")
	// This generates the inferred Pathway of an inferred RlE. It iterates, inferring parent Pathways until reaching
	// the TopLevelPathway.  Inferred Reactions are not added to the Pathway at this step.
	private void createInferredPathwayHierarchy(SimpleInstance humanEvent) throws Exception {
		List<SimpleInstance> humanPathwayReferralInstances =
			safeList(curatorToolAPI.getReferrers(humanEvent, ReactomeJavaConstants.hasEvent));

		if (humanPathwayReferralInstances.isEmpty()) {
			logger.info("Top Level Pathway inferred: " + humanEvent);
			return;
		}

		for (SimpleInstance humanPathwayReferralInstance : humanPathwayReferralInstances) {
			logger.info("Generating inferred Pathway: " + humanPathwayReferralInstance);
			// Pathways that have been inferred already are skipped, as are Pathways that are only children of the
			// Disease TopLevelPathway.
			if (hasNotBeenInferred(humanPathwayReferralInstance) &&
				!getUtils().onlyInDiseasePathway(humanPathwayReferralInstance)) {

				inferPathway(humanPathwayReferralInstance);
			}
			createInferredPathwayHierarchy(humanPathwayReferralInstance);
		}
	}

	private boolean hasNotBeenInferred(SimpleInstance sourcePathwayReferralInst) {
		return sourceInstanceToInferredInstance.get(sourcePathwayReferralInst) == null;
	}

	private void inferPathway(SimpleInstance humanPathway) throws Exception {
//		if (humanPathway.(ReactomeJavaConstants.ReactionlikeEvent)) {
//			logger.warn(humanPathway + " is a ReactionLikeEvent, which is unexpected -- refer to infer_events.pl");
//		}

		SimpleInstance inferredPathway = getUtils().createNewInferredGKInstance(humanPathway);
		inferredPathway.setAttribute(
			ReactomeJavaConstants.name, humanPathway.getAttribute(ReactomeJavaConstants.name));
		inferredPathway.setAttribute(ReactomeJavaConstants.summation, getUtils().getSummationInstance());
		//if (inferredPathway.getSchemClass().isValidAttribute(ReactomeJavaConstants.releaseDate)) {
			inferredPathway.setAttribute(ReactomeJavaConstants.releaseDate, getConfigProperties().getDateOfRelease());
		//}
		inferredPathway.setAttribute(ReactomeJavaConstants.inferredFrom, humanPathway);
		inferredPathway.setAttribute(ReactomeJavaConstants.evidenceType, getUtils().getEvidenceType());
		inferredPathway.setAttribute(ReactomeJavaConstants.goBiologicalProcess,
			humanPathway.getAttribute(ReactomeJavaConstants.goBiologicalProcess));

		inferredPathway.setAttribute(ReactomeJavaConstants.orthologousEvent, humanPathway);

		inferredPathway.setDisplayName(humanPathway.getDisplayName());
		sourceInstanceToInferredInstance.put(humanPathway, inferredPathway);
		SimpleInstance orthoStableIdentifierInst = stableIdentifierGenerator
			.generateOrthologousStableId(inferredPathway, humanPathway);
		inferredPathway.setAttribute(ReactomeJavaConstants.stableIdentifier, orthoStableIdentifierInst);
		inferredPathway = curatorToolAPI.commit(inferredPathway);

		// This was replaced with addAttributeValueIfNecessary due to a bug where a Pathway instance's
		// 'OrthologousEvent' attribute was being replaced, instead of being added to the existing array when the
		// script was executed from a jar (rather than from Eclipse) (Justin Cook 2018)
		getUtils().addOrthologousEventValue(humanPathway, inferredPathway);

		//TODO: At this point, sourcePathwayReferralInst is always a Pathway.
		// Perhaps move to its own data structure? Holdout from Perl...
		updatedInferrableHumanEvents.add(humanPathway);
	}

	// This populates the hasEvent slot of all inferred Pathways that were just generated with corresponding inferred
	// reactions
	private void addInferredEventsToInferredPathways() throws Exception {
		Set<Long> seenInferredPathway = new HashSet<>();
		for (SimpleInstance humanPathwayInst : getPathways(updatedInferrableHumanEvents)) {
			if (seenInferredPathway.contains(humanPathwayInst.getDbId())) {
				logger.info("Inferred Pathway has already been populated with inferred Events");
				continue;
			}

			SimpleInstance inferredPathwayInst = sourceInstanceToInferredInstance.get(humanPathwayInst);
			List<SimpleInstance> inferredEventInstances = getInferredEventInstances(humanPathwayInst);
			if (!isPathway(inferredPathwayInst)) {
				logger.info(humanPathwayInst + " and " +
					sourceInstanceToInferredInstance.get(humanPathwayInst) + " have different classes " +
					"(likely connected via manual inference");
				continue;
			}

			// Add inferred Events to an up-to-date copy of the inferred Pathway and commit them together. The copy
			// held in the map is stale after the commit, so it is replaced with a fresh one.
			logger.info("Adding " + inferredEventInstances.size() + " inferred Event(s) to " +
				"inferred Pathway: " + inferredPathwayInst);
			SimpleInstance currentInferredPathwayInst = curatorToolAPI.refresh(inferredPathwayInst);
			for (SimpleInstance infEventInst : inferredEventInstances) {
				currentInferredPathwayInst = getUtils().addAttributeValueIfNecessary(
					currentInferredPathwayInst, infEventInst, ReactomeJavaConstants.hasEvent);
			}
			curatorToolAPI.commit(currentInferredPathwayInst);
			sourceInstanceToInferredInstance.put(humanPathwayInst, curatorToolAPI.refresh(currentInferredPathwayInst));

			seenInferredPathway.add(humanPathwayInst.getDbId());
		}
	}

	// Collect inferred Events associated with source Event
	private List<SimpleInstance> getInferredEventInstances(SimpleInstance humanPathwayInst) {
		List<SimpleInstance> inferredEventInstances = new ArrayList<>();
		for (SimpleInstance eventInst :
			getUtils().getAttributeValues(humanPathwayInst, ReactomeJavaConstants.hasEvent)) {
			if (sourceInstanceToInferredInstance.get(eventInst) != null) {
				inferredEventInstances.add(sourceInstanceToInferredInstance.get(eventInst));
			}
		}
		return inferredEventInstances;
	}


	private void inferPrecedingEvents() throws Exception {
		Set<SimpleInstance> seenPrecedingEvent = new HashSet<>();
		for (SimpleInstance inferrableEventInst : updatedInferrableHumanEvents) {
			if (!seenPrecedingEvent.contains(inferrableEventInst)) {
				if (inferrableEventInst.getAttribute(ReactomeJavaConstants.precedingEvent)!= null) {
					logger.info("Adding preceding event to " + inferrableEventInst);
					List<SimpleInstance> precedingEventInstances = new ArrayList<>();
					// Find all preceding events for source instance that have an inferred counterpart
					for (SimpleInstance precedingEventInst :
						instanceUtilities.getAttributeValues(inferrableEventInst, ReactomeJavaConstants.precedingEvent)) {
						if (sourceInstanceToInferredInstance.get(precedingEventInst) != null) {
							precedingEventInstances.add(sourceInstanceToInferredInstance.get(precedingEventInst));
						}
					}
					// The inferred instance may have been committed since it was put in the map, so an up-to-date
					// copy is checked and added to.
					SimpleInstance inferredEventInst =
						curatorToolAPI.refresh(sourceInstanceToInferredInstance.get(inferrableEventInst));
					Set<String> inferredPrecedingEvents = new HashSet<>();
					// Find any inferred preceding events that already exist for the inferred instance
					// (don't want to add any redundant preceding events)
					for (SimpleInstance precedingEventInst :
						instanceUtilities.getAttributeValues(inferredEventInst, ReactomeJavaConstants.precedingEvent)) {
						inferredPrecedingEvents.add(precedingEventInst.getDbId().toString());
					}
					List<SimpleInstance> updatedPrecedingEventInstances = new ArrayList<>();
					// Find existing preceding events that haven't already been attached to the inferred instance
					for (SimpleInstance precedingEventInst : precedingEventInstances) {
						if (!inferredPrecedingEvents.contains(precedingEventInst.getDbId().toString())) {
							updatedPrecedingEventInstances.add(precedingEventInst);
						}
					}
					// Add preceding events to inferred instance and commit them together. The copy held in the map is
					// stale after the commit, so it is replaced with a fresh one.
					if (!updatedPrecedingEventInstances.isEmpty()) {
						for (SimpleInstance precedingEventInst : updatedPrecedingEventInstances) {
							inferredEventInst = instanceUtilities.addAttributeValueIfNecessary(
								inferredEventInst, precedingEventInst, ReactomeJavaConstants.precedingEvent);
						}
						curatorToolAPI.commit(inferredEventInst);
						sourceInstanceToInferredInstance.put(inferrableEventInst, curatorToolAPI.refresh(inferredEventInst));
					}
				}
				seenPrecedingEvent.add(inferrableEventInst);
			}
		}
	}

	private List<SimpleInstance> getPathways(List<SimpleInstance> events) {
		return events.stream().filter(this::isPathway).collect(Collectors.toList());
	}

	private boolean isPathway(SimpleInstance event) {
		return Arrays.asList(ReactomeJavaConstants.Pathway, ReactomeJavaConstants.CellLineagePath, "TopLevelPathway")
			.contains(event.getSchemaClassName());
	}

	public void setInferredEvent(Map<SimpleInstance, SimpleInstance> inferredEventCopy) {
		sourceInstanceToInferredInstance = inferredEventCopy;
	}

	private InstanceUtilities getUtils() {
		return this.instanceUtilities;
	}

	private ConfigProperties getConfigProperties() {
		return this.configProperties;
	}
}
