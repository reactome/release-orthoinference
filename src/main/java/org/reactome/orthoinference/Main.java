package org.reactome.orthoinference;

import org.reactome.orthoinference.eventsInferrer.EventsInferrer;
import org.springframework.boot.CommandLineRunner;
import org.springframework.boot.SpringApplication;
import org.springframework.boot.autoconfigure.SpringBootApplication;

@SpringBootApplication
public class Main implements CommandLineRunner {
	private EventsInferrer eventsInferrer;

	public static void main(String[] args) throws Exception {
		SpringApplication.run(Main.class, args);
	}

	public Main(EventsInferrer eventsInferrer) {
		this.eventsInferrer = eventsInferrer;
	}

	@Override
	public void run(String... args) throws Exception {
		this.eventsInferrer.inferEvents();
	}
}
