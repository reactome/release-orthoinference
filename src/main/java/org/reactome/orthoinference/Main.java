package org.reactome.orthoinference;

import java.io.FileInputStream;
import java.nio.file.Paths;
import java.util.Properties;

public class Main {

	public static void main(String[] args) throws Exception {

		String pathToConfig = args.length >= 1 ? args[0] : Paths.get("src", "main", "resources", "config.properties").toString();
		String referenceSpeciesCode = "dvir";
		String speciesCode = "zvir";

		Properties props = new Properties();
		props.load(new FileInputStream(pathToConfig));
		EventsInferrer.inferEvents(props, referenceSpeciesCode, speciesCode);
	}

}
