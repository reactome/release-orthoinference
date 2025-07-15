package org.reactome.orthoinference.eventsInferrer.services;

import org.gk.model.GKInstance;
import org.gk.persistence.MySQLAdaptor;
import org.reactome.orthoinference.ConfigProperties;
import org.reactome.orthoinference.InstanceUtilities;
import org.springframework.stereotype.Component;

import java.sql.SQLException;
import java.util.Collection;

@Component
public class DatabaseService {
    private final ConfigProperties configProperties;
    private final InstanceUtilities utils;

    public DatabaseService(ConfigProperties configProperties, InstanceUtilities utils) {
        this.configProperties = configProperties;
        this.utils = utils;
    }

    public GKInstance getHumanSpeciesInstance() throws Exception {
        Collection<GKInstance> sourceSpeciesInst = configProperties.getCurrentDBA()
            .fetchInstanceByAttribute("Species", "name", "=", "Homo sapiens");
        if (sourceSpeciesInst.isEmpty()) {
            throw new SpeciesNotFoundException("Could not find Species instance for Homo sapiens");
        }
        return sourceSpeciesInst.iterator().next();
    }

    public MySQLAdaptor getCurrentDBA() throws SQLException {
        return configProperties.getCurrentDBA();
    }

    public MySQLAdaptor getPrevDBA() throws SQLException {
        return configProperties.getPreviousDBA();
    }

    public GKInstance getSpeciesInstance() throws Exception {
        return utils.getSpeciesInstance();
    }

    public static class SpeciesNotFoundException extends Exception {
        public SpeciesNotFoundException(String message) {
            super(message);
        }
    }

}
