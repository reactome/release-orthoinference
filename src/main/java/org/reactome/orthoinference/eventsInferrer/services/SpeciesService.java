package org.reactome.orthoinference.eventsInferrer.services;

import org.gk.model.GKInstance;
import org.reactome.orthoinference.ConfigProperties;
import org.reactome.orthoinference.InstanceUtilities;

import java.io.IOException;
import org.json.simple.parser.ParseException;
import org.springframework.stereotype.Component;

@Component
public class SpeciesService {
    private final ConfigProperties configProperties;
    private final InstanceUtilities utils;

    public SpeciesService(ConfigProperties configProperties, InstanceUtilities utils) {
        this.configProperties = configProperties;
        this.utils = utils;
    }

    public String getSpeciesName(String speciesCode) throws IOException, ParseException {
        return configProperties.getSpeciesConfig().getSpeciesName(speciesCode);
    }

    public GKInstance getSpeciesInstance() throws Exception {
        return utils.getSpeciesInstance();
    }
}
