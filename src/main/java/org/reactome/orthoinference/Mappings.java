package org.reactome.orthoinference;

import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.stereotype.Component;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 * Created 5/7/2025
 */
@Component
public class Mappings {
	//private static Mappings mappings;

	private Map<String, String[]> homologueMappings;
	private Map<String, List<String>> ensgMappings;
	private Map<String, String> geneNameMappings;

	public Mappings(
		@Qualifier("sourceSpeciesCode") String sourceSpecies,
		@Qualifier("targetSpeciesCode") String targetSpecies,
		@Qualifier("pathToOrthopairs") String pathToOrthopairs) throws IOException {

		this.homologueMappings = readHomologueMappingFile(targetSpecies, sourceSpecies, pathToOrthopairs);
		this.ensgMappings = readENSGMappingFile(targetSpecies, pathToOrthopairs);
		this.geneNameMappings = readGeneNameMappingFile(targetSpecies, pathToOrthopairs);
	}

//	public static Mappings getInstance() {
//		if (mappings == null) {
//			throw new IllegalStateException("Not initialized yet");
//		}
//		return mappings;
//	}

	public Map<String, String[]> getHomologueMappings() {
		return this.homologueMappings;
	}

	public Map<String, List<String>> getEnsgMappings() {
		return this.ensgMappings;
	}

	public Map<String, String> getGeneNameMappings() {
		return this.geneNameMappings;
	}

	public List<String> getEnsEMBLGeneIdentifiers(String homologueId) {
		return getEnsgMappings().get(homologueId);
	}

	public boolean homologueHasEnsEMBLGeneMapping(String homologueId) {
		return getEnsgMappings().containsKey(homologueId);
	}

	public String getGeneName(String homologueId) {
		return getGeneNameMappings().get(homologueId);
	}

	public boolean homologueHasGeneName(String homologueId) {
		return getGeneNameMappings().containsKey(homologueId);
	}

	public String[] getHomologues(String rgpIdentifier) {
		return getHomologueMappings().get(rgpIdentifier);
	}

	// Read the species-specific orthopair 'mapping' file, and create a HashMap with the contents
	private Map<String, String[]> readHomologueMappingFile(String toSpecies, String fromSpecies, String pathToOrthopairs)
		throws IOException {

		String orthopairsFileName = fromSpecies + "_" + toSpecies + "_mapping.tsv";
		String orthopairsFilePath = Paths.get(pathToOrthopairs, orthopairsFileName).toString();
		FileReader fr = new FileReader(orthopairsFilePath);
		BufferedReader br = new BufferedReader(fr);

		Map<String, String[]> homologueMappings = new HashMap<>();
		String currentLine;
		while ((currentLine = br.readLine()) != null) {
			String[] tabSplit = currentLine.split("\t");
			String mapKey = tabSplit[0];
			String[] spaceSplit = tabSplit[1].split(" ");
			homologueMappings.put(mapKey, spaceSplit);
		}
		br.close();
		fr.close();
		return homologueMappings;
	}

	// Read the species-specific ENSG gene-protein mappings, and create a Hashmap with the contents
	private Map<String, List<String>> readENSGMappingFile(String toSpecies, String pathToOrthopairs) throws IOException {
		String mappingFileName = toSpecies + "_gene_protein_mapping.tsv";
		String mappingFilePath = Paths.get(pathToOrthopairs, mappingFileName).toString();
		FileReader fr = new FileReader(mappingFilePath);
		BufferedReader br = new BufferedReader(fr);

		Map<String, List<String>> ensgMappings = new HashMap<>();
		String currentLine;
		while ((currentLine = br.readLine()) != null) {
			String[] tabSplit = currentLine.split("\t");
			String ensgKey = tabSplit[0];
			String[] proteins = tabSplit[1].split(" ");
			for (String protein : proteins) {
				String proteinId = protein.contains(":") ? protein.split(":")[1] : protein;

				if (ensgMappings.get(proteinId) == null) {
					List<String> singleArray = new ArrayList<>();
					singleArray.add(ensgKey);
					ensgMappings.put(proteinId, singleArray);
				} else {
					ensgMappings.get(proteinId).add(ensgKey);
				}
			}
		}
		br.close();
		fr.close();

		return ensgMappings;
	}

	/**
	 * Read in the {species}_gene_name_mapping.tsv file and create a Map of UniProt identifiers to gene names.
	 * @param species String - 4-letter shortened version of species name (eg: Homo sapiens --> hsap).
	 * @return pathToOrthopairs String - Path to directory containing orthopairs files.
	 * @throws IOException - Thrown if file is not found.
	 */
	private Map<String, String> readGeneNameMappingFile(String species, String pathToOrthopairs)
		throws IOException {

		Path geneNameMappingFilePath = Paths.get(pathToOrthopairs, species + "_gene_name_mapping.tsv");
		Map<String, String> geneNameMappings = new HashMap<>();
		for (String line : Files.readAllLines(geneNameMappingFilePath)) {
			String[] tabSplit = line.split("\t");
			if (tabSplit.length == 2) {
				String uniprotId = tabSplit[0];
				String geneName = tabSplit[1];
				geneNameMappings.put(uniprotId, geneName);
			}
		}

		return geneNameMappings;
	}
}
