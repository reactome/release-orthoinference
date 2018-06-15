package org.reactome.orthoinference;

import java.util.Collection;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.gk.schema.InvalidAttributeException;
import org.gk.schema.InvalidAttributeValueException;
import org.gk.schema.SchemaClass;

public class GenerateInstance {
	
		private static MySQLAdaptor dba; 
		private static GKInstance speciesInst = null;
	
		public static void setAdaptor(MySQLAdaptor dbAdaptor)
		{
			dba = dbAdaptor;
		}
		
		public void setSpeciesInst(GKInstance speciesInstCopy)
		{
			speciesInst = speciesInstCopy;
		}
		
		// Creates new instance that will be inferred based on the incoming instances class		
		public static GKInstance newInferredGKInstance(GKInstance instanceToBeInferred)
//		TODO: Instance Edits; Valid Attribute comparment/species; check_intracellular 
		{
			GKInstance inferredInst = null;
			String reactionClass = instanceToBeInferred.getSchemClass().getName();
			SchemaClass referenceDNAClass = dba.getSchema().getClassByName(reactionClass);
			inferredInst = new GKInstance(referenceDNAClass);
			return inferredInst;
		}
		
		// create_ghost equivalent; Returns a mock homologue that is needed in cases of unsuccessful inference
		public static GKInstance newMockGKInstance(GKInstance instanceToBeMocked) throws InvalidAttributeException, InvalidAttributeValueException, Exception
		{
			SchemaClass geeClass = dba.getSchema().getClassByName(ReactomeJavaConstants.GenomeEncodedEntity);
			GKInstance mockedInst = new GKInstance(geeClass);
			String mockedName = (String) instanceToBeMocked.getAttributeValue("name");
			mockedInst.addAttributeValue(ReactomeJavaConstants.name, "Ghost homologue of " + mockedName);
			mockedInst.addAttributeValue(ReactomeJavaConstants.species, speciesInst);
			//TODO: Instance edit; check intracellular; inferred to/from; update;
			return mockedInst;
		}
		
		// Checks that equivalent instances don't already exist in the DB, substituting if they do
		//TODO: Go over the Perl version and make sure they match perfectly
		public static GKInstance checkForIdenticalInstances(GKInstance inferredInst) throws Exception
		{
			@SuppressWarnings("unchecked")
			Collection<GKInstance> identicalInstances = dba.fetchIdenticalInstances(inferredInst);
			if (identicalInstances.size() == 1) {
				return identicalInstances.iterator().next();
			} else if (identicalInstances.size() > 1) {
				return identicalInstances.iterator().next();
			} else {
				return inferredInst;
			}
		}
}
