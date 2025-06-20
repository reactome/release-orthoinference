package org.reactome.orthoinference;

import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.boot.CommandLineRunner;
import org.springframework.stereotype.Component;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 * Created 6/8/2025
 */
@Component
public class Runner implements CommandLineRunner {

	private ConfigProperties configProperties;
	private String speciesCode;
	private EventsInferrer eventsInferrer;

	public Runner(ConfigProperties configProperties, @Qualifier("targetSpeciesCode") String speciesCode, EventsInferrer eventsInferrer) {
		this.configProperties = configProperties;
		this.speciesCode = speciesCode;
		this.eventsInferrer = eventsInferrer;
	}

	@Override
	public void run(String... args) throws Exception {
		System.out.println("Hello from CommandLineRunner!");
		System.out.println("Species code: " + speciesCode);
		System.out.println("Person id: " + configProperties.getPersonId());

		//EventsInferrer eventsInferrer = new EventsInferrer(configProperties, getSpeciesCode());
		eventsInferrer.inferEvents();

	}
}
