package org.reactome.orthoinference.eventsInferrer.services;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.reactome.orthoinference.PathwaysInferrer;
import org.reactome.orthoinference.ReactionInferrer;
import org.springframework.stereotype.Component;

import java.util.*;

@Component
public class InferenceServices {
    private static final Logger logger = LogManager.getLogger();

    private final ReactionInferrer reactionInferrer;
    private final PathwaysInferrer pathwaysInferrer;
    private final PathwayDiagramService pathwayDiagramService;
    private final SpeciesService speciesService;

    private static Map<GKInstance,GKInstance> manualEventToNonHumanSource = new HashMap<>();
    private static List<GKInstance> manualHumanEvents = new ArrayList<>();

    public InferenceServices(
        ReactionInferrer reactionInferrer,
        PathwaysInferrer pathwaysInferrer,
        PathwayDiagramService pathwayDiagramService,
        SpeciesService speciesService
    ) {
        this.reactionInferrer = reactionInferrer;
        this.pathwaysInferrer = pathwaysInferrer;
        this.pathwayDiagramService = pathwayDiagramService;
        this.speciesService = speciesService;
    }

    public void inferReactions(Map<Long, GKInstance> reactionMap) throws Exception {
        for (Map.Entry<Long, GKInstance> entry : reactionMap.entrySet()) {
            GKInstance reactionInst = entry.getValue();
            logger.info("Attempting RlE inference: " + reactionInst);

            if (isAlreadyInferred(reactionInst)) {
                continue;
            }

            inferReactionWithErrorHandling(reactionInst);
        }
    }

    public void generatePathwayDiagrams(
        GKInstance humanSpeciesInstance,
        MySQLAdaptor currentDBA,
        MySQLAdaptor prevDBA
    ) throws Exception {
        pathwayDiagramService.generateOrthologousPathwayDiagrams(currentDBA, prevDBA, humanSpeciesInstance);
    }

    public void inferPathways() throws Exception {
        pathwaysInferrer.setInferredEvent(reactionInferrer.getInferredEvent());
        pathwaysInferrer.inferPathways(reactionInferrer.getInferrableHumanEvents());
    }

    public int getEligibleCount() {
        return reactionInferrer.getEligibleCount();
    }

    public int getInferredCount() {
        return reactionInferrer.getInferredCount();
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
            reactionInferrer.addAlreadyInferredEvents(reactionInst, prevInfInst);
        } else {
            logger.info("Inferred RlE already exists, skipping inference");
            manualEventToNonHumanSource.put(reactionInst, prevInfInst);
            manualHumanEvents.add(reactionInst);
        }
        return true;
    }


    @SuppressWarnings("unchecked")
    private List<GKInstance> checkIfPreviouslyInferred(
        GKInstance reactionInst,
        String attribute,
        List<GKInstance> previouslyInferredInstances
    ) throws Exception {
        GKInstance speciesInstance = speciesService.getSpeciesInstance();

        for (GKInstance attributeInst : (Collection<GKInstance>) reactionInst.getAttributeValuesList(attribute)) {
            GKInstance reactionSpeciesInst = (GKInstance) attributeInst.getAttributeValue(ReactomeJavaConstants.species);
            if (reactionSpeciesInst.getDBID() == speciesInstance.getDBID() &&
                    attributeInst.getAttributeValue(ReactomeJavaConstants.isChimeric) == null) {
                previouslyInferredInstances.add(attributeInst);
            }
        }
        return previouslyInferredInstances;
    }

    private void inferReactionWithErrorHandling(GKInstance reactionInst) {
        try {
            reactionInferrer.inferReaction(reactionInst);
            logger.info("Successfully inferred " + reactionInst);
        } catch (Exception e) {
            logger.error("Failed to infer reaction: " + reactionInst, e);
            throw new RuntimeException("Failed to infer reaction", e);
        }
    }
}
