package org.reactome.release;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.neo4j.driver.v1.*;

import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

public class Main {
	private static final Logger logger = LogManager.getLogger();

	public static void main( String[] args ) throws IOException {
		logger.info("Beginning Main step...");

		String pathToResources = args.length > 0 ? args[0] : "ncbi/src/main/resources/sample_config.properties";
		Properties props = new Properties();
		try {
			props.load(new FileInputStream(pathToResources));
		} catch (IOException e) {
			e.printStackTrace();
		}
		int version = Integer.parseInt(props.getProperty("reactomeVersion", "67"));
		String outputDir = props.getProperty("outputDir", "archive");
		Files.createDirectories(Paths.get(outputDir));
		logger.info("Files for Reactome version " + version + " will be output to the directory " + outputDir);

		Session graphDBSession = getGraphDBDriver(props).session();

		logger.info("Generating UniProt accession to Main Gene mapping");
		List<NCBIEntry> ncbiEntries = NCBIEntry.getUniProtToNCBIGeneMap(graphDBSession);

		logger.info("Writing proteins_version file");
		NCBIGene.getInstance(ncbiEntries, outputDir, version).writeProteinFile();

		int numGeneXMLFiles = Integer.parseInt(props.getProperty("numGeneXMLFiles", "1"));
		logger.info("Writing gene XML file(s)");
		NCBIGene.getInstance(ncbiEntries, outputDir, version).writeGeneXMLFiles(
			graphDBSession,

			numGeneXMLFiles);

		logger.info("Writing NCBI protein file");
		NCBIProtein.writeNCBIProteinFile(outputDir, version, ncbiEntries);

		logger.info("Writing UCSC files");
		UCSC.writeUCSCFiles(graphDBSession);

		graphDBSession.close();
		logger.info("Finished Main step");
		System.exit(0);
	}

	private static Driver getGraphDBDriver(Properties props) {
		String host = props.getProperty("host","localhost");
		String port = props.getProperty("port", Integer.toString(7687));
		String user = props.getProperty("user", "neo4j");
		String password = props.getProperty("password", "root");

		return GraphDatabase.driver("bolt://" + host + ":" + port, AuthTokens.basic(user, password));
	}
}