package org.reactome.orthoinference;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.ReactomeJavaConstants;
import org.json.simple.parser.ParseException;
import org.reactome.curation.model.SimpleInstance;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

import java.io.IOException;
import java.util.*;

/*
 *  All PhysicalEntitys, ReactionlikeEvents and Pathways are routed to this class to generate their stable identifiers
 */
@Component
public class StableIdentifierGenerator {
    private static final Logger logger = LogManager.getLogger();
    private static Map<String,Integer> seenOrthoIds = new HashMap<>();

    private static final String IDENTIFIER_VERSION = "1";
    private static final String STABLE_ID_ATTRIBUTE = "stId";

    private ConfigProperties configProperties;
    private String speciesCode;
    private CuratorToolAPI curatorToolAPI;

    public StableIdentifierGenerator(
        ConfigProperties configProperties,
        @Qualifier("targetSpeciesCode") String speciesCode,
        CuratorToolAPI curatorToolAPI
    ) {
        this.configProperties = configProperties;
        this.speciesCode = speciesCode;
        this.curatorToolAPI = curatorToolAPI;
    }

    public SimpleInstance generateOrthologousStableId(SimpleInstance inferredInst, SimpleInstance originalInst)
        throws Exception {

        // Sometimes there already exists a StableIdentifier value for an instance, if there are multiple instances
        // that can create one instance.
        if (inferredInst.getAttribute(ReactomeJavaConstants.stableIdentifier) != null) {
            return null;
        }

        String sourceIdentifier = getSourceIdentifier(originalInst);

        logger.info("Generating orthologous stable identifier for " + sourceIdentifier);

        SimpleInstance orthoStableIdentifierInst = getOrthoStableIdentiferInstance(sourceIdentifier);

        // Populate inferred instance with new StableIdentifier instance
        logger.info("Stable identifier generated: " + orthoStableIdentifierInst.getDisplayName());

        return orthoStableIdentifierInst;
    }

    private SimpleInstance getOrthoStableIdentiferInstance(String sourceIdentifier) throws Exception {
        String targetIdentifier = getTargetIdentifier(sourceIdentifier);

        // Check that the stable identifier instance does not already exist in DB
        SimpleInstance existingStableIdentifier = getCuratorToolAPI().fetchStableIdentifier(targetIdentifier);
        if (existingStableIdentifier == null) {
            return createAndStoreOrthoStableIdentifier(targetIdentifier);
        }
        return existingStableIdentifier;
    }

    private String getTargetIdentifier(String sourceIdentifier) throws Exception {
        // For now, Human is hard-coded as the source species, so we replace the stableIdentifier source species
        // based on that assumption
        String targetIdentifier = sourceIdentifier.replace("HSA", getSpeciesAbbreviation());
        // Paralogs will have the same base stable identifier, but we want to denote when that happens.
        // We pull the value from `seenOrthoIds`, increment it and then add it to the stable identifier name
        // (eg: R-MMU-123456-2)
        int paralogCount = Optional.ofNullable(seenOrthoIds.get(targetIdentifier)).orElse(0) + 1;
        seenOrthoIds.put(targetIdentifier, paralogCount);
        if (paralogCount > 1) {
            targetIdentifier += "-" + paralogCount;
        }

        return targetIdentifier;
    }

    private String getSourceIdentifier(SimpleInstance originalInst) {
        // Human PhysicalEntitys and Events have a StableIdentifier instance in the stableIdentifier attribute, and
        // its identifier is also kept on the instance itself as stId. Instances in the graph database do not all
        // have the stableIdentifier relationship, so stId is used when there is no StableIdentifier instance.
        SimpleInstance stableIdentifierInst =
            (SimpleInstance) originalInst.getAttribute(ReactomeJavaConstants.stableIdentifier);
        if (stableIdentifierInst != null) {
            return getIdentifier(stableIdentifierInst);
        }

        Object stId = originalInst.getAttribute(STABLE_ID_ATTRIBUTE);
        if (stId == null) {
            logAndThrow("No stable identifier found for " + originalInst);
        }
        return stId.toString();
    }

    private String getIdentifier(SimpleInstance stableIdentifierInst) {
        // An instance holds the instances it refers to as shells -- dbId, displayName and schemaClass only -- so the
        // StableIdentifier has to be inflated before its identifier can be read.
        String sourceIdentifier = (String) stableIdentifierInst.getAttribute(ReactomeJavaConstants.identifier);
        if (sourceIdentifier == null) {
            SimpleInstance inflatedStableIdentifierInst = getCuratorToolAPI().inflate(stableIdentifierInst);
            if (inflatedStableIdentifierInst != null) {
                sourceIdentifier =
                    (String) inflatedStableIdentifierInst.getAttribute(ReactomeJavaConstants.identifier);
            }
        }
        if (sourceIdentifier == null) {
            logAndThrow("No identifier found for stable identifier instance " + stableIdentifierInst);
        }
        return sourceIdentifier;
    }

    private SimpleInstance createAndStoreOrthoStableIdentifier(String targetIdentifier) {
        SimpleInstance orthoStableIdentifierInst = createOrthologousStableIdentifierInstance(targetIdentifier);
        // Commit copies the dbId the database assigns back onto orthoStableIdentifierInst.
        getCuratorToolAPI().commit(orthoStableIdentifierInst);

        return orthoStableIdentifierInst;
    }

    // Generates a new stable identifier instance
    private SimpleInstance createOrthologousStableIdentifierInstance(String targetIdentifier) {
        SimpleInstance orthoStableIdentifierInst = new SimpleInstance();
        orthoStableIdentifierInst.setSchemaClassName(ReactomeJavaConstants.StableIdentifier);

        orthoStableIdentifierInst.setAttribute(ReactomeJavaConstants.identifier, targetIdentifier);
        orthoStableIdentifierInst.setAttribute(ReactomeJavaConstants.identifierVersion, IDENTIFIER_VERSION);

        // Same display name InstanceDisplayNameGenerator gives a StableIdentifier: identifier.identifierVersion
        orthoStableIdentifierInst.setDisplayName(targetIdentifier + "." + IDENTIFIER_VERSION);

        return orthoStableIdentifierInst;
    }

    private CuratorToolAPI getCuratorToolAPI() {
        return this.curatorToolAPI;
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
