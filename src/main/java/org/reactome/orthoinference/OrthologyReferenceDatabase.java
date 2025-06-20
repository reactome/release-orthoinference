package org.reactome.orthoinference;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.persistence.MySQLAdaptor;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.ParseException;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

import java.io.IOException;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 * Created 5/5/2025
 */
@Component
public class OrthologyReferenceDatabase {

	private ConfigProperties configProperties;
	private String speciesCode;
	private InstanceUtilities instanceUtilities;

	private GKInstance ensemblDBInst;
	private GKInstance uniprotDbInst;

	private GKInstance alternateDbInst;
	private boolean altRefDbExists;
	private String altRefDbId;

	public OrthologyReferenceDatabase(
		ConfigProperties configProperties,
		@Qualifier("targetSpeciesCode") String speciesCode,
		InstanceUtilities instanceUtilities
	) throws Exception {
		this.configProperties = configProperties;
		this.speciesCode = speciesCode;
		this.instanceUtilities = instanceUtilities;

		setEnsEMBLDatabase();
		fetchAndSetUniprotDbInstance();
		setAlternateReferenceDatabase();
	}

	public GKInstance getReferenceDatabase(String homologSource) {
		return homologSource.equals("ENSP") ? ensemblDBInst : uniprotDbInst;
	}

	public String getReferenceDatabaseSourceName(String homologSource) {
		return homologSource.equals("ENSP") ? "ENSEMBL" : "UniProt";
	}

	public GKInstance getEnsemblDBInst() {
		return this.ensemblDBInst;
	}

	public GKInstance getAlternateReferenceDatabase() {
		return this.alternateDbInst;
	}

	public boolean alternateReferenceDatabaseExists() {
		return this.altRefDbExists;
	}

	public String getAlternateReferenceDbId() {
		return this.altRefDbId;
	}

	private void setEnsEMBLDatabase() throws Exception {
		if (isFungiSpecies()) {
			fetchOrCreateFungiEnsemblDatabase();
		} else if (isProtistSpecies()) {
			fetchOrCreateProtistsEnsemblDatabase();
		} else {
			fetchOrCreateMainEnsemblDatabase();
		}
	}

	// Fetches Uniprot DB instance
	@SuppressWarnings("unchecked")
	private void fetchAndSetUniprotDbInstance() throws Exception {
		Collection<GKInstance> uniprotDbInstances = (Collection<GKInstance>)
			getCurrentDBA().fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceDatabase, ReactomeJavaConstants.name, "=", "UniProt");
		uniprotDbInst = uniprotDbInstances.iterator().next();
	}

	private void setAlternateReferenceDatabase() throws Exception {
		JSONObject altRefDbJSON = getAltRefDBJSON();
		if (altRefDbJSON != null) {
			createAlternateReferenceDBInstance(altRefDbJSON);
		} else {
			setAltRefDbToFalse();
		}
	}

	private void fetchOrCreateProtistsEnsemblDatabase() throws Exception {
		final String dbName = "EnsEMBL Protists";
		final String dbBaseUrl = "https://protists.ensembl.org/";

		fetchOrCreateEnsemblDatabase(dbName, dbBaseUrl);
	}

	private void fetchOrCreateFungiEnsemblDatabase() throws Exception {
		final String dbName = "EnsEMBL Fungi";
		final String dbBaseUrl = "https://fungi.ensembl.org/";

		fetchOrCreateEnsemblDatabase(dbName, dbBaseUrl);
	}

	public void fetchOrCreateMainEnsemblDatabase() throws Exception {
		final String dbName = "ENSEMBL";
		final String dbBaseUrl = "https://ensembl.org/";

		fetchOrCreateEnsemblDatabase(dbName, dbBaseUrl);
	}

	private boolean isProtistSpecies() {
		final List<String> protistsSpeciesCodes = Arrays.asList("ddis", "pfal");
		return protistsSpeciesCodes.contains(getSpeciesCode());
	}

	private boolean isFungiSpecies() {
		final List<String> fungiSpeciesCodes = Arrays.asList("scer", "spom");
		return fungiSpeciesCodes.contains(getSpeciesCode());
	}

	private void fetchOrCreateEnsemblDatabase(String dbName, String dbBaseUrl)
		throws Exception {

		Collection<GKInstance> ensemblDBInstances =
			(Collection<GKInstance>) getCurrentDBA().fetchInstanceByAttribute(ReactomeJavaConstants.ReferenceDatabase, ReactomeJavaConstants.name, "=", dbName);
		if (ensemblDBInstances == null || ensemblDBInstances.isEmpty()) {
			ensemblDBInst = createAndStoreEnsEMBLDatabase(dbName, dbBaseUrl);
		} else {
			ensemblDBInst = ensemblDBInstances.iterator().next();
		}
	}

	private GKInstance createAndStoreEnsEMBLDatabase(String dbName, String dbBaseUrl)
		throws Exception {

		GKInstance ensEMBLDbInst = new GKInstance(getCurrentDBA().getSchema().getClassByName(ReactomeJavaConstants.ReferenceDatabase));

		ensEMBLDbInst.setDbAdaptor(getCurrentDBA());
		ensEMBLDbInst.addAttributeValue(ReactomeJavaConstants.created, getInstanceEdit());
		ensEMBLDbInst.addAttributeValue(ReactomeJavaConstants.name, dbName);
		ensEMBLDbInst.addAttributeValue(ReactomeJavaConstants.url, dbBaseUrl);
		ensEMBLDbInst.addAttributeValue(ReactomeJavaConstants.accessUrl, dbBaseUrl + "id/###ID###");
		ensEMBLDbInst.setDisplayName(dbName);
		getCurrentDBA().storeInstance(ensEMBLDbInst);
		return ensEMBLDbInst;
	}

	private JSONObject getAltRefDBJSON() throws IOException, ParseException {
		return getConfigProperties().getSpeciesConfig().getAltRefDbJSON(getSpeciesCode());
	}

	// Create instance pertaining to any alternative reference DB for the species
	private void createAlternateReferenceDBInstance(JSONObject altRefDbJSON) throws Exception {
		alternateDbInst = new GKInstance(getCurrentDBA().getSchema().getClassByName(ReactomeJavaConstants.ReferenceDatabase));
		alternateDbInst.setDbAdaptor(getCurrentDBA());
		alternateDbInst.addAttributeValue(ReactomeJavaConstants.created, getInstanceEdit());
		alternateDbInst.addAttributeValue(ReactomeJavaConstants.name, ((JSONArray) altRefDbJSON.get("dbname")).get(0));
		alternateDbInst.addAttributeValue(ReactomeJavaConstants.url, altRefDbJSON.get("url"));
		alternateDbInst.addAttributeValue(ReactomeJavaConstants.accessUrl, altRefDbJSON.get("access"));
		alternateDbInst.setAttributeValue(ReactomeJavaConstants._displayName, ((JSONArray) altRefDbJSON.get("dbname")).get(0));
		alternateDbInst = instanceUtilities.checkForIdenticalInstances(alternateDbInst, null);
		if (altRefDbJSON.get("alt_id") != null)
		{
			altRefDbId = (String) altRefDbJSON.get("alt_id");
		}
		altRefDbExists = true;
	}

	private void setAltRefDbToFalse() {
		altRefDbExists = false;
	}

	private GKInstance getInstanceEdit() throws Exception {
		return this.instanceUtilities.getInstanceEdit();
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
