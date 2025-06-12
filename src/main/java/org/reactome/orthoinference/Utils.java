package org.reactome.orthoinference;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.SchemaClass;
import org.json.simple.parser.ParseException;
import org.reactome.release.common.database.InstanceEditUtils;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.context.annotation.Bean;
import org.springframework.stereotype.Component;

import java.io.IOException;
import java.sql.SQLException;

@Component
public class Utils {

    private ConfigProperties configProperties;
    private String speciesCode;

    private GKInstance speciesInst;
    private GKInstance instanceEdit;
    private GKInstance evidenceType;
    private GKInstance summation;

    private StableIdentifierGenerator stableIdentifierGenerator;
    private InstanceUtilities instanceUtilities;

    public Utils(
        ConfigProperties configProperties,
        @Qualifier("targetSpeciesCode") String speciesCode,
        StableIdentifierGenerator stableIdentifierGenerator,
        InstanceUtilities instanceUtilities
    ) {
        this.configProperties = configProperties;
        this.speciesCode = speciesCode;

        this.stableIdentifierGenerator = stableIdentifierGenerator;
        this.instanceUtilities = instanceUtilities;
    }



    public GKInstance getEvidenceType() throws Exception {
        if (evidenceType == null) {
            String evidenceTypeText = "inferred by electronic annotation";

            evidenceType = new GKInstance(
                    getCurrentDBA().getSchema().getClassByName(ReactomeJavaConstants.EvidenceType));
            evidenceType.setDbAdaptor(getCurrentDBA());
            evidenceType.addAttributeValue(ReactomeJavaConstants.created, instanceUtilities.getInstanceEdit());
            evidenceType.addAttributeValue(ReactomeJavaConstants.name, evidenceTypeText);
            evidenceType.addAttributeValue(ReactomeJavaConstants.name, "IEA");
            evidenceType.setDisplayName(evidenceTypeText);
            evidenceType = instanceUtilities.checkForIdenticalInstances(evidenceType, null);
        }
        return evidenceType;
    }

    public GKInstance getSummationInstance() throws Exception {
        if (summation == null) {
            summation = new GKInstance(getCurrentDBA().getSchema().getClassByName(ReactomeJavaConstants.Summation));
            summation.setDbAdaptor(getCurrentDBA());
            summation.addAttributeValue(ReactomeJavaConstants.created, instanceUtilities.getInstanceEdit());
            String summationText = "This event has been computationally inferred from an event that has been" +
                " demonstrated in another species.<p>The inference is based on the homology mapping from PANTHER." +
                " Briefly, reactions for which all involved PhysicalEntities (in input, output and catalyst) have a" +
                " mapped orthologue/paralogue (for complexes at least 75% of components must have a mapping) are" +
                " inferred to the other species. High level events are also inferred for these events to allow for" +
                " easier navigation.<p><a href='/electronic_inference_compara.html' target = 'NEW'>More details and" +
                " caveats of the event inference in Reactome.</a> For details on PANTHER see also:" +
                " <a href='http://www.pantherdb.org/about.jsp' target='NEW'>http://www.pantherdb.org/about.jsp</a>";
            summation.addAttributeValue(ReactomeJavaConstants.text, summationText);
            summation.setDisplayName(summationText);
            summation = instanceUtilities.checkForIdenticalInstances(summation, null);
        }
        return summation;
    }

    public StableIdentifierGenerator getStableIdentifierGenerator() {
        return stableIdentifierGenerator;
    }

    private ConfigProperties getConfigProperties() {
        return this.configProperties;
    }

    private String getSpeciesCode() {
        return this.speciesCode;
    }

    private MySQLAdaptor getCurrentDBA() throws SQLException {
        return getConfigProperties().getCurrentDBA();
    }
}
