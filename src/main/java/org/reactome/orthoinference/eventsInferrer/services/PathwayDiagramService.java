package org.reactome.orthoinference.eventsInferrer.services;

import org.gk.model.GKInstance;
import org.gk.persistence.MySQLAdaptor;
import org.reactome.orthoinference.ConfigProperties;
import org.reactome.orthoinference.OrthologousPathwayDiagramGenerator;
import org.springframework.stereotype.Component;

@Component
public class PathwayDiagramService {
    private final ConfigProperties configProperties;
    private final SpeciesService speciesService;

    public PathwayDiagramService(ConfigProperties configProperties, SpeciesService speciesService) {
        this.configProperties = configProperties;
        this.speciesService = speciesService;
    }

    public void generateOrthologousPathwayDiagrams(
        MySQLAdaptor currentDBA,
        MySQLAdaptor prevDBA,
        GKInstance humanSpeciesInstance
    ) throws Exception {
        OrthologousPathwayDiagramGenerator generator = new OrthologousPathwayDiagramGenerator(
                currentDBA,
                prevDBA,
                speciesService.getSpeciesInstance(),
                humanSpeciesInstance,
                configProperties.getPersonId()
        );
        generator.generateOrthologousPathwayDiagrams();
    }
}
