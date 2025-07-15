package org.reactome.orthoinference.eventsInferrer;

import java.util.*;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;

import org.reactome.orthoinference.*;
import org.reactome.orthoinference.eventsInferrer.services.DatabaseService;
import org.reactome.orthoinference.eventsInferrer.services.InferenceServices;
import org.reactome.orthoinference.eventsInferrer.services.ReportingService;
import org.reactome.orthoinference.eventsInferrer.services.SpeciesService;
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

	private String speciesCode;

	private DatabaseService databaseService;
	private ReportingService reportingService;
	private SpeciesService speciesService;
	private InferenceServices inferenceServices;

	private static OrthologousPathwayDiagramGenerator orthologousPathwayDiagramGenerator;

	public EventsInferrer(
		@Qualifier("targetSpeciesCode") String speciesCode,
		DatabaseService databaseService,
		ReportingService reportingService,
		SpeciesService speciesService,
		InferenceServices inferenceServices
	) {
		this.speciesCode = speciesCode;

		this.databaseService = databaseService;
		this.reportingService = reportingService;
		this.speciesService = speciesService;
		this.inferenceServices = inferenceServices;
	}

	@SuppressWarnings("unchecked")
	public void inferEvents() throws Exception {
		String speciesName = speciesService.getSpeciesName(speciesCode);
		logger.info("Beginning orthoinference of " + speciesName);

		GKInstance humanSpeciesInstance = databaseService.getHumanSpeciesInstance();
		Map<Long, GKInstance> reactionMap = getHumanReactionMap(humanSpeciesInstance);
		inferenceServices.inferReactions(reactionMap);

		inferenceServices.inferPathways();
		inferenceServices.generatePathwayDiagrams(humanSpeciesInstance, databaseService.getCurrentDBA(), databaseService.getPrevDBA());

		reportingService.outputReport(
			speciesCode,
			inferenceServices.getEligibleCount(),
			inferenceServices.getInferredCount()
		);

		logger.info("Finished orthoinference of " + speciesName);
	}

	private Map<Long, GKInstance> getHumanReactionMap(GKInstance humanSpeciesInstance) throws Exception {
		Collection<GKInstance> reactionInstances = (Collection<GKInstance>)
			databaseService.getCurrentDBA().fetchInstanceByAttribute("ReactionlikeEvent", "species", "=", humanSpeciesInstance.getDBID());

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
}
