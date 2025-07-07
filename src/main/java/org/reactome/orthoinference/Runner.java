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

	private EventsInferrer eventsInferrer;

	public Runner(EventsInferrer eventsInferrer) {

		this.eventsInferrer = eventsInferrer;
	}

	@Override
	public void run(String... args) throws Exception {
		eventsInferrer.inferEvents();
	}
}
