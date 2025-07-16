package org.reactome.orthoinference;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.InstanceDisplayNameGenerator;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.SchemaClass;
import org.json.simple.parser.ParseException;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

import java.io.IOException;
import java.sql.SQLException;
import java.util.*;

/*
 *  All PhysicalEntitys, ReactionlikeEvents and Pathways are routed to this class to generate their stable identifiers
 */
@Component
public class StableIdentifierGenerator {
    private static final Logger logger = LogManager.getLogger();
    private static Map<String,Integer> orthoIdsCount = new HashMap<>();

    private ConfigProperties configProperties;
    private String speciesCode;

    public StableIdentifierGenerator(ConfigProperties configProperties, @Qualifier("targetSpeciesCode") String speciesCode) {
        this.configProperties = configProperties;
        this.speciesCode = speciesCode;
    }

    public GKInstance generateOrthologousStableId(GKInstance inferredInstance, GKInstance sourceInstance)
        throws Exception {

        // Sometimes there already exists a StableIdentifier value for an instance, if there are multiple instances
        // that can create one instance.
        if (inferredInstanceAlreadyHasStableIdentifier(inferredInstance)) {
            return null;
        }

        // All Human PhysicalEntities and Events will have a StableIdentifier instance in the stableIdentifier
        // attribute
        return getOrthoStableIdentiferInstance(getSourceStableIdentifier(sourceInstance));
    }

    private boolean inferredInstanceAlreadyHasStableIdentifier(GKInstance inferredInstance) throws Exception {
        return inferredInstance.getAttributeValue(ReactomeJavaConstants.stableIdentifier) != null;
    }

    private GKInstance getSourceStableIdentifier(GKInstance sourceInstance) throws Exception {
        GKInstance sourceStableIdentifier =
            (GKInstance) sourceInstance.getAttributeValue(ReactomeJavaConstants.stableIdentifier);

        if (sourceStableIdentifier == null) {
            logAndThrow("No stable identifier instance found for " + sourceInstance);
        }

        return sourceStableIdentifier;
    }

    private GKInstance getOrthoStableIdentiferInstance(GKInstance sourceStableIdentifier) throws Exception {
        logger.info("Generating orthologous stable identifier for " + sourceStableIdentifier.getDisplayName());

        String orthoStableIdIdentifier = getOrthoStableIdIdentifier(sourceStableIdentifier);

        // Check that the stable identifier instance does not already exist in DB
        List<GKInstance> existingStableIdentifiers = getExistingStableIdentifiers(orthoStableIdIdentifier);

        GKInstance orthoStableIdentifierInstance;
        if (existingStableIdentifiers.isEmpty()) {
            orthoStableIdentifierInstance =
                createAndStoreOrthoStableIdentifier(sourceStableIdentifier, orthoStableIdIdentifier);
        } else {
            orthoStableIdentifierInstance = existingStableIdentifiers.get(0);
        }
        return orthoStableIdentifierInstance;
    }

    private String getOrthoStableIdIdentifier(GKInstance sourceStableIdentifier) throws Exception {
        // For now, Human is hard-coded as the source species, so we replace the stableIdentifier source species
        // based on that assumption
        String sourceIdentifier = (String) sourceStableIdentifier.getAttributeValue(ReactomeJavaConstants.identifier);
        String targetIdentifier = sourceIdentifier.replace("HSA", getSpeciesAbbreviation());
        // Paralogs will have the same base stable identifier, but we want to denote when that happens.
        // We pull the value from `orthoIdsCount`, increment it and then add it to the stable identifier name
        // (eg: R-MMU-123456-2)
        targetIdentifier = addParalogCountIfNecessary(targetIdentifier);

        return targetIdentifier;
    }

    private String addParalogCountIfNecessary(String targetIdentifier) {
        int paralogCount = getParalogCount(targetIdentifier);
        orthoIdsCount.put(targetIdentifier, paralogCount);
        if (paralogCount > 1) {
            targetIdentifier += "-" + paralogCount;
        }
        return targetIdentifier;
    }

    private int getParalogCount(String targetIdentifier) {
        return Optional.ofNullable(orthoIdsCount.get(targetIdentifier)).orElse(0) + 1;
    }

    private List<GKInstance> getExistingStableIdentifiers(String targetIdentifier) throws Exception {
        Collection<GKInstance> existingStableIdentifiers = (Collection<GKInstance>)
            getDBA().fetchInstanceByAttribute("StableIdentifier", "identifier", "=", targetIdentifier);
        if (existingStableIdentifiers == null) {
            return new ArrayList<>();
        }
        return new ArrayList<>(existingStableIdentifiers);
    }

    private GKInstance createAndStoreOrthoStableIdentifier(GKInstance stableIdentifierInst, String targetIdentifier)
        throws Exception {

        GKInstance orthoStableIdentifierInst =
            createOrthologousStableIdentifierInstance(targetIdentifier);
        getDBA().storeInstance(orthoStableIdentifierInst);
        logger.info("Stable identifier generated: " + orthoStableIdentifierInst.getDisplayName());

        return orthoStableIdentifierInst;
    }

    // Generates a new stable identifier instance
    private GKInstance createOrthologousStableIdentifierInstance(String targetIdentifier) throws Exception {
        SchemaClass instanceClass = getDBA().getSchema().getClassByName(ReactomeJavaConstants.StableIdentifier);
        GKInstance orthoStableIdentifierInst = new GKInstance(instanceClass);
        orthoStableIdentifierInst.setDbAdaptor(getDBA());

        orthoStableIdentifierInst.addAttributeValue(ReactomeJavaConstants.identifier, targetIdentifier);
        orthoStableIdentifierInst.addAttributeValue(ReactomeJavaConstants.identifierVersion, "1");

        InstanceDisplayNameGenerator.setDisplayName(orthoStableIdentifierInst);

        return orthoStableIdentifierInst;
    }


    private MySQLAdaptor getDBA() throws SQLException {
        return this.configProperties.getCurrentDBA();
    }

    private String getSpeciesAbbreviation() throws IOException, ParseException {
        return this.configProperties.getSpeciesConfig().getAbbreviation(getSpeciesCode());
    }

    private String getSpeciesCode() {
        return this.speciesCode;
    }

    private void logAndThrow(String errorMessage) {
        logger.fatal(errorMessage);
        throw new RuntimeException(errorMessage);
    }
}
