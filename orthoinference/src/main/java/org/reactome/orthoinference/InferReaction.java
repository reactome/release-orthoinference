package org.reactome.orthoinference;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.gk.model.GKInstance;
import org.gk.model.Instance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.InvalidAttributeException;

public class InferReaction {

	private static MySQLAdaptor dba;
	private static String eligibleFilehandle;
	private static GKInstance summationInst;
	private static GKInstance evidenceTypeInst;
	
	// This function mimics the Perl version of InferEvent, inferring attribute instances of input, output, catalyst activity, and regulations
	public static void inferEvent(GKInstance reactionInst) throws InvalidAttributeException, Exception
	{
		String dbId = reactionInst.getAttributeValue("DB_ID").toString();
		String stableId = reactionInst.getAttributeValue("name").toString();
		System.out.println("Reaction: [" + dbId + "] " + stableId);	
		
		if (SkipTests.skipInstance(reactionInst))
		{
			return;
		}

		// TODO: Release date/instance edit; %inferred_event (not till the end); %being_inferred; Global variables for summation/evidence type (rather than recreating it each time); 75% threshold as variable or hard-coded?
		// Creates inferred instance of reaction
		GKInstance infReactionInst = GenerateInstance.newInferredGKInstance(reactionInst);
		infReactionInst.addAttributeValue(ReactomeJavaConstants.name, reactionInst.getAttributeValuesList(ReactomeJavaConstants.name));
		infReactionInst.addAttributeValue(ReactomeJavaConstants.goBiologicalProcess, reactionInst.getAttributeValue(ReactomeJavaConstants.goBiologicalProcess));
		infReactionInst.addAttributeValue(ReactomeJavaConstants.summation, summationInst);
		infReactionInst.addAttributeValue(ReactomeJavaConstants.evidenceType, evidenceTypeInst);
		
		List<Integer> reactionProteinCounts = ProteinCount.countDistinctProteins(reactionInst);
//		System.out.println("Overall Counts: " + reactionProteinCounts);
		
		//TODO: Success measure; Verify no return needed; count_leaves (not till the end)
		if (reactionProteinCounts.get(0) > 0) {
			String eligibleEvent = reactionInst.getAttributeValue(ReactomeJavaConstants.DB_ID).toString() + "\t" + reactionInst.getDisplayName() + "\n";	
			Files.write(Paths.get(eligibleFilehandle), eligibleEvent.getBytes(), StandardOpenOption.APPEND);
			if (InferReaction.inferAttributes(reactionInst, infReactionInst, "input"))
			{
				if (InferReaction.inferAttributes(reactionInst, infReactionInst, "output"))
				{
					
				}
			}
	//		InferReaction.inferAttributes(reactionInst, inferredReaction, "output");
	//		InferReaction.inferCatalyst(reactionInst, inferredReaction);
		} 
	}
	
	// Function used to create inferred instances related to either 'input' or 'output'
	@SuppressWarnings("unchecked")
	public static boolean inferAttributes(GKInstance reactionInst, Instance infReactionInst, String attribute) throws InvalidAttributeException, Exception
	{
		ArrayList<GKInstance> infAttributeInstances = new ArrayList<GKInstance>();
		//TODO: Return true if there isn't any null inferred instances
		for (GKInstance attributeInst : (Collection<GKInstance>) reactionInst.getAttributeValuesList(attribute))
		{
			GKInstance infAttributeInst = OrthologousEntity.createOrthoEntity(attributeInst, false);
			if (infAttributeInst == null)
			{
				return false;
			} 
			infAttributeInstances.add(infAttributeInst);
		}
		infReactionInst.addAttributeValue(attribute, infAttributeInstances);
		return true;
	}
	
	// Function used to created inferred catalysts
	@SuppressWarnings("unchecked")
	public static void inferCatalyst(GKInstance reactionInst, Instance inferredReaction) throws InvalidAttributeException, Exception
	{
		for (@SuppressWarnings("unused") GKInstance attributeInst : (Collection<GKInstance>) reactionInst.getAttributeValuesList("catalystActivity"))
		{
//			System.out.println("  CatalystActivity");
		}
	}
	
	public static void setAdaptor(MySQLAdaptor dbAdaptor)
	{
		dba = dbAdaptor;
	}
	
	public static void setEligibleFilename(String eligibleFilename)
	{
		eligibleFilehandle = eligibleFilename;
	}
	
	public static void setEvidenceTypeInst() throws Exception
	{
		evidenceTypeInst = new GKInstance(dba.getSchema().getClassByName(ReactomeJavaConstants.EvidenceType));
		evidenceTypeInst.setDbAdaptor(dba);
		evidenceTypeInst.addAttributeValue(ReactomeJavaConstants.name, "inferred by electronic annotation");
		evidenceTypeInst.addAttributeValue(ReactomeJavaConstants.name, "IEA");
		evidenceTypeInst = GenerateInstance.checkForIdenticalInstances(evidenceTypeInst);
	}
	
	public static void setSummationInst() throws Exception
	{
		summationInst = new GKInstance(dba.getSchema().getClassByName(ReactomeJavaConstants.Summation));
		summationInst.setDbAdaptor(dba);
		summationInst.addAttributeValue(ReactomeJavaConstants.text, "This event has been computationally inferred from an event that has been demonstrated in another species.<p>The inference is based on the homology mapping in Ensembl Compara. Briefly, reactions for which all involved PhysicalEntities (in input, output and catalyst) have a mapped orthologue/paralogue (for complexes at least 75% of components must have a mapping) are inferred to the other species. High level events are also inferred for these events to allow for easier navigation.<p><a href='/electronic_inference_compara.html' target = 'NEW'>More details and caveats of the event inference in Reactome.</a> For details on the Ensembl Compara system see also: <a href='http://www.ensembl.org/info/docs/compara/homology_method.html' target='NEW'>Gene orthology/paralogy prediction method.</a>");
		summationInst = GenerateInstance.checkForIdenticalInstances(summationInst);
	}
	
}
