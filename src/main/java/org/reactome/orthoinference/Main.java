package org.reactome.orthoinference;

import java.io.FileInputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Properties;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class Main {

	private static final Logger logger = LogManager.getLogger();
	
	public static void main(String[] args) throws Exception {

		String pathToConfig = Paths.get("src", "main", "resources", "config.properties").toString();
		String speciesCode = "";
		if (args.length == 2) {
			pathToConfig = args[0];
			speciesCode = args[1];
		} else if (args.length == 1 && args[0].length() == 4) {
			speciesCode = args[0];
		} else {
			logger.fatal("Please include a 4-letter species code as the first argument (eg: mmus)");
//			System.exit(0);
		}
		String referenceSpeciesCode = "cov1";
		speciesCode="cov2";

		Properties props = new Properties();
		props.load(new FileInputStream(pathToConfig));
		EventsInferrer.inferEvents(props, referenceSpeciesCode, speciesCode);
	}

}
