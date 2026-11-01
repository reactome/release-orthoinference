package org.reactome.orthoinference;

import org.gk.model.ReactomeJavaConstants;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.ParseException;
import org.reactome.curation.model.SimpleInstance;
import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

import java.io.IOException;
import java.util.Arrays;
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
	private CuratorToolAPI curatorToolAPI;

	private SimpleInstance ensemblDBInst;
	private SimpleInstance uniprotDbInst;

	private SimpleInstance alternateDbInst;
	private boolean altRefDbExists;
	private String altRefDbId;

	public OrthologyReferenceDatabase(
		ConfigProperties configProperties,
		@Qualifier("targetSpeciesCode") String speciesCode,
		InstanceUtilities instanceUtilities,
		CuratorToolAPI curatorToolAPI
	) throws Exception {
		this.configProperties = configProperties;
		this.speciesCode = speciesCode;
		this.instanceUtilities = instanceUtilities;
		this.curatorToolAPI = curatorToolAPI;

		setEnsEMBLDatabase();
		fetchAndSetUniprotDbInstance();
		setAlternateReferenceDatabase();
	}

	public SimpleInstance getReferenceDatabase(String homologSource) {
		return homologSource.equals("ENSP") ? ensemblDBInst : uniprotDbInst;
	}

	public String getReferenceDatabaseSourceName(String homologSource) {
		return homologSource.equals("ENSP") ? "ENSEMBL" : "UniProt";
	}

	public SimpleInstance getEnsemblDBInst() {
		return this.ensemblDBInst;
	}

	public SimpleInstance getAlternateReferenceDatabase() {
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
	private void fetchAndSetUniprotDbInstance() {
		uniprotDbInst = getCuratorToolAPI().fetchReferenceDatabase("UniProt");
		if (uniprotDbInst == null) {
			throw new IllegalStateException("No UniProt ReferenceDatabase instance found");
		}
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

	private void fetchOrCreateEnsemblDatabase(String dbName, String dbBaseUrl) throws Exception {
		ensemblDBInst = getCuratorToolAPI().fetchReferenceDatabase(dbName);

		if (ensemblDBInst == null) {
			ensemblDBInst = createAndStoreEnsEMBLDatabase(dbName, dbBaseUrl);
		}
	}

	private SimpleInstance createAndStoreEnsEMBLDatabase(String dbName, String dbBaseUrl) {
		SimpleInstance ensEMBLDbInst = new SimpleInstance();
		ensEMBLDbInst.setSchemaClassName(ReactomeJavaConstants.ReferenceDatabase);

		ensEMBLDbInst.setAttribute(ReactomeJavaConstants.name, dbName);
		ensEMBLDbInst.setAttribute(ReactomeJavaConstants.url, dbBaseUrl);
		ensEMBLDbInst.setAttribute(ReactomeJavaConstants.accessUrl, dbBaseUrl + "id/###ID###");
		ensEMBLDbInst.setDisplayName(dbName);
		return getCuratorToolAPI().commit(ensEMBLDbInst);
	}

	private JSONObject getAltRefDBJSON() throws IOException, ParseException {
		return getConfigProperties().getSpeciesConfig().getAltRefDbJSON(getSpeciesCode());
	}

	// Create instance pertaining to any alternative reference DB for the species
	private void createAlternateReferenceDBInstance(JSONObject altRefDbJSON) throws Exception {
		alternateDbInst = new SimpleInstance();
		alternateDbInst.setSchemaClassName(ReactomeJavaConstants.ReferenceDatabase);
			//new GKInstance(getCurrentDBA().getSchema().getClassByName(ReactomeJavaConstants.ReferenceDatabase));

		alternateDbInst.setAttribute(ReactomeJavaConstants.name, ((JSONArray) altRefDbJSON.get("dbname")).get(0));
		alternateDbInst.setAttribute(ReactomeJavaConstants.url, altRefDbJSON.get("url"));
		alternateDbInst.setAttribute(ReactomeJavaConstants.accessUrl, altRefDbJSON.get("access"));
		alternateDbInst.setDisplayName(((JSONArray) altRefDbJSON.get("dbname")).get(0).toString());
		alternateDbInst = instanceUtilities.checkForIdenticalInstances(alternateDbInst, null);
		if (altRefDbJSON.get("alt_id") != null) {
			altRefDbId = (String) altRefDbJSON.get("alt_id");
		}
		altRefDbExists = true;
	}

	private void setAltRefDbToFalse() {
		altRefDbExists = false;
	}

	private ConfigProperties getConfigProperties() {
		return this.configProperties;
	}

	private CuratorToolAPI getCuratorToolAPI() {
		return this.curatorToolAPI;
	}

	private String getSpeciesCode() {
		return this.speciesCode;
	}
}
