package org.reactome.orthoinference;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Properties;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class Main {
	private static final Logger logger = LogManager.getLogger();

	private static ConfigProperties configProperties;

	@Parameter(names = {"-c", "-config"}, description = "Path to configuration file")
	private String pathToConfig = Paths.get("src", "main", "resources", "config.properties").toString();

	@Parameter(
		names = {"-s", "-species"},
		description = "Four letter species code as target inference (e.g. mmus)",
		required = true
	)
	private String speciesCode;

	public static void main(String[] args) throws Exception {
		Main main = new Main();
		JCommander.newBuilder()
			.addObject(main)
			.build()
			.parse(args);
		main.run();
	}

	public static ConfigProperties getConfigProperties() {
		return configProperties;
	}

	private void run() throws Exception {
		Properties props = new Properties();
		props.load(Files.newInputStream(Paths.get(getPathToConfig())));
		configProperties = new ConfigProperties(props);

		EventsInferrer eventsInferrer = new EventsInferrer(configProperties, getSpeciesCode());
		eventsInferrer.inferEvents();
	}

	private String getPathToConfig() {
		return this.pathToConfig;
	}

	private String getSpeciesCode() {
		return this.speciesCode;
	}
}
