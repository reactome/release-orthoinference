package org.reactome.orthoinference;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;

import java.io.FileReader;
import java.io.IOException;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *         Created 1/9/2024
 */
public class SpeciesConfig {
    private String pathToSpeciesConfig;

    public SpeciesConfig(String pathToSpeciesConfig) {
        this.pathToSpeciesConfig = pathToSpeciesConfig;
    }

    public String getSpeciesName(String speciesCode) throws IOException, ParseException {
        JSONArray speciesNames = (JSONArray) getSpeciesObject(speciesCode).get("name");
        return (String) speciesNames.get(0);
    }

    public String getAbbreviation(String speciesCode) throws IOException, ParseException {
        return (String) getSpeciesObject(speciesCode).get("abbreviation");
    }

    public JSONObject getAltRefDbJSON(String speciesCode) throws IOException, ParseException {
        return (JSONObject) getSpeciesObject(speciesCode).get("alt_refdb");
    }

    private JSONObject getSpeciesObject(String speciesCode)
        throws IOException, ParseException {

        JSONParser parser = new JSONParser();
        Object obj = parser.parse(new FileReader(pathToSpeciesConfig));
        JSONObject jsonObject = (JSONObject) obj;

        // Parse Species information (found in Species.json config file)
        return (JSONObject) jsonObject.get(speciesCode);
    }

    private String getPathToSpeciesConfig() {
        return this.pathToSpeciesConfig;
    }
}
