package org.reactome.orthoinference;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.GKSchemaClass;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

import java.util.*;

@Component
public class ComplexPolymerInferrer {
    private static final Logger logger = LogManager.getLogger();

    private static Map<GKInstance, GKInstance> complexPolymerIdenticals = new HashMap<>();
    private static Map<String,GKInstance> complexIdenticals = new HashMap<>();

    private MySQLAdaptor dba;
    private OrthoEntityInferrer orthoEntityInferrer;
    private InstanceUtilities instanceUtilities;

    private GKInstance complexSummationInst;

    public ComplexPolymerInferrer(
        @Qualifier("currentDBA") MySQLAdaptor dba,
        OrthoEntityInferrer orthoEntityInferrer,
        InstanceUtilities instanceUtilities
    ) {
        this.dba = dba;
        this.orthoEntityInferrer = orthoEntityInferrer;
        this.instanceUtilities = instanceUtilities;
    }

    // Infers Complex or Polymer instances. These instances are generally comprised of more than 1 PhysicalEntity,
    // and calls 'createOrthoEntity' for each one. Complex/Polymer instances are also subject to the
    // 'countDistinctProteins' function. The result from this needs to have at least 75% of total proteins to be
    // inferrable for inference to continue.
    public GKInstance createInfComplexPolymer(GKInstance complexInst, boolean override) throws Exception {
        if (complexPolymerIdenticals.containsKey(complexInst)) {
            logger.info("Inferred Complex/Polymer already exists");
            return complexPolymerIdenticals.get(complexInst);
        }

        if (!hasValidProteinCounts(complexInst, override)) {
            return null;
        }

        GKInstance infComplexInst = createInitialInferredInstance(complexInst);
        if (!processComponents(complexInst, infComplexInst)) {
            return null;
        }

        infComplexInst = handleCachingAndUpdates(infComplexInst, complexInst);

        if (override) {
            return infComplexInst;
        }

        complexPolymerIdenticals.put(complexInst, infComplexInst);
        return infComplexInst;
    }


    private boolean hasValidProteinCounts(GKInstance complexInst, boolean override) throws Exception {
        List<Integer> complexProteinCounts = getProteinCountUtility().getDistinctProteinCounts(complexInst);
        int complexTotalProteinCounts = complexProteinCounts.get(0);
        int complexInferrableProteinCounts = complexProteinCounts.get(1);

        int percent = calculateProteinPercentage(complexTotalProteinCounts, complexInferrableProteinCounts);
        int percentThreshold = 75;

        if (!override && !meetsProteinThreshold(complexTotalProteinCounts, complexInferrableProteinCounts, percent, percentThreshold)) {
            logger.info(
                    "Complex/Polymer protein count is below " + percentThreshold + "% threshold (" + percent + "%) -- " +
                            "terminating inference"
            );
            return false;
        }

        logger.info("Complex protein counts. Total: " + complexTotalProteinCounts +
                "  Inferrable: " + complexInferrableProteinCounts);
        return true;
    }

    private int calculateProteinPercentage(int total, int inferrable) {
        if (total > 0) {
            return (inferrable * 100) / total;
        }
        return 0;
    }

    private boolean meetsProteinThreshold(int total, int inferrable, int percent, int percentThreshold) {
        return !(total > 0 && inferrable == 0) && percent >= percentThreshold;
    }

    private GKInstance createInitialInferredInstance(GKInstance complexInst) throws Exception {
        GKInstance infComplexInst = instanceUtilities.createNewInferredGKInstance(complexInst);
        infComplexInst.addAttributeValue(ReactomeJavaConstants.summation, getComplexSummationInstance());
        infComplexInst.addAttributeValue(ReactomeJavaConstants.name,
                complexInst.getAttributeValue(ReactomeJavaConstants.name));
        return infComplexInst;
    }

    private boolean processComponents(GKInstance complexInst, GKInstance infComplexInst) throws Exception {
        List<GKInstance> infComponentInstances = new ArrayList<>();

        if (complexInst.getSchemClass().isa(ReactomeJavaConstants.Complex)) {
            processComplexComponents(complexInst, infComplexInst, infComponentInstances);
        } else if (complexInst.getSchemClass().isa(ReactomeJavaConstants.Polymer)) {
            processPolymerComponents(complexInst, infComplexInst, infComponentInstances);
        } else {
            logger.warn(complexInst + " is not a Complex or a Polymer");
            return false;
        }

        infComplexInst.setAttributeValue(ReactomeJavaConstants._displayName,
                complexInst.getAttributeValue(ReactomeJavaConstants._displayName));
        return true;
    }

    private void processComplexComponents(GKInstance complexInst, GKInstance infComplexInst,
                                          List<GKInstance> infComponentInstances) throws Exception {
        Collection<GKInstance> componentInstances =
                complexInst.getAttributeValuesList(ReactomeJavaConstants.hasComponent);
        logger.info("Complex components: " + componentInstances);

        for (GKInstance componentInst : componentInstances) {
            infComponentInstances.add(this.orthoEntityInferrer.createOrthoEntity(componentInst, true));
        }
        infComplexInst.addAttributeValue(ReactomeJavaConstants.hasComponent, infComponentInstances);
    }

    private void processPolymerComponents(GKInstance complexInst, GKInstance infComplexInst,
                                          List<GKInstance> infComponentInstances) throws Exception {
        Collection<GKInstance> repeatedUnitInstances =
            complexInst.getAttributeValuesList(ReactomeJavaConstants.repeatedUnit);
        logger.info("Polymer repeated units: " + repeatedUnitInstances);

        for (GKInstance repeatedUnitInst : repeatedUnitInstances) {
            infComponentInstances.add(this.orthoEntityInferrer.createOrthoEntity(repeatedUnitInst, true));
        }
        infComplexInst.addAttributeValue(ReactomeJavaConstants.repeatedUnit, infComponentInstances);
    }

    private GKInstance handleCachingAndUpdates(GKInstance infComplexInst, GKInstance complexInst) throws Exception {
        String cacheKey = instanceUtilities.getCacheKey(
                (GKSchemaClass) infComplexInst.getSchemClass(), infComplexInst);

        if (complexIdenticals.containsKey(cacheKey)) {
            return complexIdenticals.get(cacheKey);
        }

        infComplexInst = instanceUtilities.checkForIdenticalInstances(infComplexInst, complexInst);
        complexIdenticals.put(cacheKey, infComplexInst);

        updateInferredAttributes(infComplexInst, complexInst);

        return infComplexInst;
    }

    private void updateInferredAttributes(GKInstance infComplexInst, GKInstance complexInst) throws Exception {
        infComplexInst = instanceUtilities.addAttributeValueIfNecessary(
                infComplexInst, complexInst, ReactomeJavaConstants.inferredFrom);
        dba.updateInstanceAttribute(infComplexInst, ReactomeJavaConstants.inferredFrom);

        complexInst = instanceUtilities.addAttributeValueIfNecessary(
                complexInst, infComplexInst, ReactomeJavaConstants.inferredTo);
        dba.updateInstanceAttribute(complexInst, ReactomeJavaConstants.inferredTo);
    }

    private GKInstance getComplexSummationInstance() throws Exception {
        if (complexSummationInst == null) {
            complexSummationInst = new GKInstance(dba.getSchema().getClassByName(ReactomeJavaConstants.Summation));
            complexSummationInst.setDbAdaptor(dba);
            complexSummationInst.addAttributeValue(ReactomeJavaConstants.created, getInstanceEdit());
            String complexSummationText = "This complex/polymer has been computationally inferred (based on PANTHER) " +
                    "from a complex/polymer involved in an event that has been demonstrated in another species.";
            complexSummationInst.addAttributeValue(ReactomeJavaConstants.text, complexSummationText);
            complexSummationInst.setAttributeValue(ReactomeJavaConstants._displayName, complexSummationText);
            complexSummationInst = instanceUtilities.checkForIdenticalInstances(complexSummationInst, null);
        }
        return complexSummationInst;
    }

    private GKInstance getInstanceEdit() throws Exception {
        return this.instanceUtilities.getInstanceEdit();
    }

    private ProteinCountUtility getProteinCountUtility() {
        return this.instanceUtilities.getProteinCountUtility();
    }
}
