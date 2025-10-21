package org.reactome.orthoinference;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.GKSchemaClass;
import org.gk.schema.SchemaClass;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

import java.util.*;

@Component
public class EntitySetInferrer {
    private static final Logger logger = LogManager.getLogger();

    private static Map<GKInstance, GKInstance> inferredEntitySetIdenticals = new HashMap<>();
    private static Map<String,GKInstance> entitySetIdenticals = new HashMap<>();

    private OrthoEntityInferrer orthoEntityInferrer;
    private InstanceUtilities instanceUtilities;
    private MySQLAdaptor dba;

    public EntitySetInferrer(
        OrthoEntityInferrer orthoEntityInferrer,
        InstanceUtilities instanceUtilities,
        @Qualifier("currentDBA") MySQLAdaptor dba
    ) {
        this.orthoEntityInferrer = orthoEntityInferrer;
        this.instanceUtilities = instanceUtilities;
        this.dba = dba;
    }

    // EntitySet inference function. This function will initially call createOrthoEntity on all 'members' before
    // filtering by the type of EntitySet (Open, Candidate, or Defined Sets) and completing a specific inference.
    // Important to note is that while there are multiple cases where createOrthoEntity is called (for members and
    // candidates) in createInfEntitySet, the override functionality is not used here.  Presumably, this is because
    // the instances aren't a constituent part of a single instance (as in Complexes), but rather are stand-alone ones
    // that also happen to be included in a Set.  This means they should be subject to the stringency of a typical
    // instance, rather than using override to create mock instances that allow an instance to be inferred more easily.
    @SuppressWarnings("unchecked")
    public GKInstance createInfEntitySet(GKInstance entitySetInst, boolean override) throws Exception {
        if (inferredEntitySetIdenticals.containsKey(entitySetInst)) {
            logger.info("Inferred EntitySet already exists");
            return inferredEntitySetIdenticals.get(entitySetInst);
        }

        List<GKInstance> infMembersList = inferMembers(entitySetInst);
        logMembersInference(entitySetInst, infMembersList);

        // Begin inference of EntitySet
        GKInstance infEntitySetInst = initializeInferredEntitySet(entitySetInst, infMembersList);

        // Validate protein counts
        if (!validateProteinCounts(entitySetInst, override)) {
            return null;
        }

        // Process based on EntitySet type
        infEntitySetInst = processEntitySetByType(entitySetInst, infEntitySetInst, infMembersList, override);
        if (infEntitySetInst == null) {
            return null;
        }

        // Set display name and handle caching
        infEntitySetInst = finalizeEntitySetInference(entitySetInst, infEntitySetInst, override);

        return infEntitySetInst;
    }

    private List<GKInstance> inferMembers(GKInstance entitySetInst) throws Exception {
        Set<String> existingMemberInstances = new HashSet<>();
        List<GKInstance> infMembersList = new ArrayList<>();
        Collection<GKInstance> memberInstances = entitySetInst.getAttributeValuesList(ReactomeJavaConstants.hasMember);

        for (GKInstance memberInst : memberInstances) {
            GKInstance infMemberInst = this.orthoEntityInferrer.createOrthoEntity(memberInst, false);
            if (isValidInferredMember(infMemberInst, existingMemberInstances)) {
                existingMemberInstances.add(infMemberInst.getAttributeValue(ReactomeJavaConstants.name).toString());
                infMembersList.add(infMemberInst);
            }
        }
        return infMembersList;
    }

    private void logMembersInference(GKInstance entitySetInst, List<GKInstance> infMembersList) throws Exception {
        if (!entitySetInst.getSchemClass().isa(ReactomeJavaConstants.CandidateSet)) {
            Collection<GKInstance> memberInstances = entitySetInst.getAttributeValuesList(ReactomeJavaConstants.hasMember);
            logger.info("Total member instances: " + memberInstances.size());
            logger.info("Member instances: " + memberInstances);
            logger.info("Total number of inferred members: " + infMembersList.size() + "/" + memberInstances.size());
        }
    }

    private GKInstance initializeInferredEntitySet(GKInstance entitySetInst, List<GKInstance> infMembersList) throws Exception {
        GKInstance infEntitySetInst = instanceUtilities.createNewInferredGKInstance(entitySetInst);
        infEntitySetInst.addAttributeValue(ReactomeJavaConstants.name,
                entitySetInst.getAttributeValuesList(ReactomeJavaConstants.name));
        infEntitySetInst.addAttributeValue(ReactomeJavaConstants.hasMember, infMembersList);
        return infEntitySetInst;
    }

    private boolean validateProteinCounts(GKInstance entitySetInst, boolean override) throws Exception {
        List<Integer> entitySetProteinCounts = getProteinCountUtility().getDistinctProteinCounts(entitySetInst);
        int entitySetTotalCount = entitySetProteinCounts.get(0);
        int entitySetInferrableCount = entitySetProteinCounts.get(1);

        if (!override && entitySetTotalCount > 0 && entitySetInferrableCount == 0) {
            logger.info("No distinct proteins found in EntitySet -- terminating inference");
            return false;
        }
        return true;
    }

    private GKInstance processEntitySetByType(GKInstance entitySetInst, GKInstance infEntitySetInst,
                                              List<GKInstance> infMembersList, boolean override) throws Exception {

        if (entitySetInst.getSchemClass().isa(ReactomeJavaConstants.CandidateSet)) {
            return processCandidateSet(entitySetInst, infEntitySetInst, infMembersList, override);
        } else if (entitySetInst.getSchemClass().isa(ReactomeJavaConstants.DefinedSet)) {
            return processDefinedSet(entitySetInst, infEntitySetInst, infMembersList, override);
        }
        return infEntitySetInst;
    }

    private GKInstance processCandidateSet(GKInstance entitySetInst, GKInstance infEntitySetInst,
                                           List<GKInstance> infMembersList, boolean override) throws Exception {

        List<GKInstance> infCandidatesList = inferCandidates(entitySetInst, infMembersList);

        if (!infCandidatesList.isEmpty()) {
            infEntitySetInst.addAttributeValue(ReactomeJavaConstants.hasCandidate, infCandidatesList);
            return infEntitySetInst;
        }

        return handleEmptyCandidates(entitySetInst, infEntitySetInst, infMembersList, override);
    }

    private GKInstance processDefinedSet(GKInstance entitySetInst, GKInstance infEntitySetInst,
                                         List<GKInstance> infMembersList, boolean override) throws Exception {

        if (infMembersList.isEmpty()) {
            if (override) {
                logger.info("Mock DefinedSet instance needed");
                return instanceUtilities.createMockGKInstance(entitySetInst);
            }
            logger.info("No member instances found -- terminating inference");
            return null;
        }

        if (infMembersList.size() == 1) {
            logger.info("Only 1 member from EntitySet was inferred, converting to PE: " + infMembersList.get(0));
            return infMembersList.get(0);
        }

        return infEntitySetInst;
    }

    private List<GKInstance> inferCandidates(GKInstance entitySetInst, List<GKInstance> infMembersList) throws Exception {
        Set<String> existingCandidateInstances = new HashSet<>();
        List<GKInstance> infCandidatesList = new ArrayList<>();

        Collection<GKInstance> candidateInstances =
                entitySetInst.getAttributeValuesList(ReactomeJavaConstants.hasCandidate);

        logger.info("Total candidate instances: " + candidateInstances.size());
        logger.info("Candidate instances: " + candidateInstances);

        for (GKInstance candidateInst : candidateInstances) {
            GKInstance infCandidateInst = this.orthoEntityInferrer.createOrthoEntity(candidateInst, false);
            if (isValidCandidate(infCandidateInst, infMembersList, existingCandidateInstances)) {
                existingCandidateInstances.add(infCandidateInst.getAttributeValue(ReactomeJavaConstants.name).toString());
                infCandidatesList.add(infCandidateInst);
            }
        }

        logger.info("Total number of inferred candidates: " +
                infCandidatesList.size() + "/" + candidateInstances.size());

        return infCandidatesList;
    }

    private boolean isValidCandidate(GKInstance infCandidateInst,
                                     List<GKInstance> infMembersList,
                                     Set<String> existingCandidateInstances) throws Exception {
        if (infCandidateInst == null) {
            return false;
        }

        String candidateName = infCandidateInst.getAttributeValue(ReactomeJavaConstants.name).toString();

        // Check if the candidate name doesn't exist in members or other candidates
        boolean notInMembers = infMembersList.stream()
                .noneMatch(member -> {
                    try {
                        return member.getAttributeValue(ReactomeJavaConstants.name).toString().equals(candidateName);
                    } catch (Exception e) {
                        logger.error("Error comparing member names", e);
                        return true; // Conservatively assume it matches to avoid duplicates
                    }
                });

        return notInMembers && !existingCandidateInstances.contains(candidateName);
    }

    private GKInstance handleEmptyCandidates(GKInstance entitySetInst,
                                             GKInstance infEntitySetInst,
                                             List<GKInstance> infMembersList,
                                             boolean override) throws Exception {
        if (infMembersList.isEmpty()) {
            if (override) {
                logger.info("Mock CandidateSet instance needed");
                return instanceUtilities.createMockGKInstance(entitySetInst);
            }
            return null;
        }

        if (infMembersList.size() == 1) {
            return infMembersList.get(0);
        }

        logger.info("No candidates inferred, but there are inferred members. Converting to DefinedSet");
        return convertToDefinedSet(entitySetInst, infEntitySetInst, infMembersList);
    }

    private GKInstance convertToDefinedSet(GKInstance entitySetInst,
                                           GKInstance infEntitySetInst,
                                           List<GKInstance> infMembersList) throws Exception {
        SchemaClass definedSetClass = dba.getSchema().getClassByName(ReactomeJavaConstants.DefinedSet);
        GKInstance infDefinedSetInst = new GKInstance(definedSetClass);
        infDefinedSetInst.setDbAdaptor(dba);

        // Set basic attributes
        infDefinedSetInst.addAttributeValue(ReactomeJavaConstants.created, getInstanceEdit());
        infDefinedSetInst.setAttributeValue(ReactomeJavaConstants.name,
                infEntitySetInst.getAttributeValuesList(ReactomeJavaConstants.name));
        infDefinedSetInst.setAttributeValue(ReactomeJavaConstants.hasMember, infMembersList);

        // Handle compartment
        if (hasValidCompartment(entitySetInst)) {
            addCompartmentToDefinedSet(entitySetInst, infDefinedSetInst);
        }

        infDefinedSetInst.addAttributeValue(ReactomeJavaConstants.species, getSpeciesInstance());

        logger.info("Successfully converted to DefinedSet");
        return infDefinedSetInst;
    }

    private boolean hasValidCompartment(GKInstance entitySetInst) throws Exception {
        return entitySetInst.getSchemClass().isValidAttribute(ReactomeJavaConstants.compartment) &&
                entitySetInst.getAttributeValue(ReactomeJavaConstants.compartment) != null;
    }

    private void addCompartmentToDefinedSet(GKInstance entitySetInst,
                                            GKInstance infDefinedSetInst) throws Exception {
        for (Object compartmentInst : entitySetInst.getAttributeValuesList(ReactomeJavaConstants.compartment)) {
            GKInstance compartmentInstGk = (GKInstance) compartmentInst;
            if (compartmentInstGk.getSchemClass().isa(ReactomeJavaConstants.Compartment)) {
                infDefinedSetInst.addAttributeValue(ReactomeJavaConstants.compartment, compartmentInstGk);
            } else {
                GKInstance newCompartmentInst =
                    instanceUtilities.createCompartmentInstance(compartmentInstGk);
                infDefinedSetInst.addAttributeValue(ReactomeJavaConstants.compartment, newCompartmentInst);
            }
        }
    }

    private GKInstance finalizeEntitySetInference(GKInstance entitySetInst, GKInstance infEntitySetInst,
                                                  boolean override) throws Exception {

        infEntitySetInst.setAttributeValue(ReactomeJavaConstants._displayName,
                entitySetInst.getAttributeValue(ReactomeJavaConstants._displayName));

        infEntitySetInst = handleCaching(entitySetInst, infEntitySetInst);
        updateEntitySetInferredAttributes(entitySetInst, infEntitySetInst);

        if (!override) {
            inferredEntitySetIdenticals.put(entitySetInst, infEntitySetInst);
        }

        return infEntitySetInst;
    }

    private boolean isValidInferredMember(GKInstance infMemberInst, Set<String> existingMemberInstances) throws Exception {
        return infMemberInst != null &&
                !existingMemberInstances.contains(infMemberInst.getAttributeValue(ReactomeJavaConstants.name).toString());
    }

    private GKInstance handleCaching(GKInstance entitySetInst, GKInstance infEntitySetInst) throws Exception {
        String cacheKey = instanceUtilities.getCacheKey((GKSchemaClass) infEntitySetInst.getSchemClass(), infEntitySetInst);
        if (entitySetIdenticals.containsKey(cacheKey)) {
            return entitySetIdenticals.get(cacheKey);
        }

        infEntitySetInst = instanceUtilities.checkForIdenticalInstances(infEntitySetInst, entitySetInst);
        entitySetIdenticals.put(cacheKey, infEntitySetInst);
        return infEntitySetInst;
    }

    private void updateEntitySetInferredAttributes(GKInstance entitySetInst, GKInstance infEntitySetInst) throws Exception {
        if (infEntitySetInst.getSchemClass().isValidAttribute(ReactomeJavaConstants.species) &&
                entitySetInst.getAttributeValue(ReactomeJavaConstants.species) != null) {

            infEntitySetInst = instanceUtilities.addAttributeValueIfNecessary(
                    infEntitySetInst, entitySetInst, ReactomeJavaConstants.inferredFrom);
            dba.updateInstanceAttribute(infEntitySetInst, ReactomeJavaConstants.inferredFrom);

            entitySetInst = instanceUtilities.addAttributeValueIfNecessary(
                    entitySetInst, infEntitySetInst, ReactomeJavaConstants.inferredTo);
            dba.updateInstanceAttribute(entitySetInst, ReactomeJavaConstants.inferredTo);
        }
    }

    private GKInstance getInstanceEdit() throws Exception {
        return this.instanceUtilities.getInstanceEdit();
    }

    private ProteinCountUtility getProteinCountUtility() {
        return this.instanceUtilities.getProteinCountUtility();
    }

    private GKInstance getSpeciesInstance() throws Exception {
        return this.instanceUtilities.getSpeciesInstance();
    }
}
