package org.reactome.orthoinference;

import org.gk.model.ReactomeJavaConstants;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.reactome.curation.model.SimpleInstance;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import static org.reactome.orthoinference.InferenceTestHarness.*;

/**
 * Runs whole inferences through {@link EventsInferrer#inferEvents()} -- every human ReactionlikeEvent in the
 * in-memory database, then the pathways that hold them -- and checks the inferred pathway hierarchy, preceding events,
 * the handling of reactions inferred before, and the report. The homology fixture is described in
 * {@link InferenceTestHarness}; PROTn (n < 98) proteins have homologues and PROT99 has none.
 */
class EventAndPathwayInferenceScenarioTest {

	@TempDir
	Path tempDir;

	private InferenceTestHarness harness;
	private InMemoryCuratorToolAPI api;
	private HumanFixtures fx;

	@BeforeEach
	void setUp() throws Exception {
		harness = new InferenceTestHarness(tempDir, null);
		api = harness.api;
		fx = harness.fx;
	}

	@AfterEach
	void deleteOutputFiles() throws Exception {
		InferenceTestHarness.deleteOutputFiles();
	}

	@Nested
	@DisplayName("Pathway hierarchy")
	class PathwayHierarchy {

		@Test
		void inferredReactionsArePlacedInInferredCopiesOfTheirPathways() throws Exception {
			SimpleInstance r1 = mappedReaction("R1");
			SimpleInstance r2 = mappedReaction("R2");
			SimpleInstance glycolysis = fx.pathway(ReactomeJavaConstants.Pathway, "Glycolysis", r1, r2);
			SimpleInstance metabolism = fx.pathway(TOP_LEVEL_PATHWAY, "Metabolism", glycolysis);

			runInference();

			SimpleInstance inferredGlycolysis = inferred(ReactomeJavaConstants.Pathway, glycolysis);
			assertThat(inferredGlycolysis.getDisplayName(), is("Glycolysis"));
			assertThat(inferredGlycolysis.getAttribute(ReactomeJavaConstants.name),
				is(Collections.singletonList("Glycolysis")));
			assertThat(dbIds(harness.values(inferredGlycolysis, ReactomeJavaConstants.species)),
				contains(harness.targetSpecies.getDbId()));
			assertThat(harness.single(inferredGlycolysis, ReactomeJavaConstants.evidenceType).getDisplayName(),
				is("inferred by electronic annotation"));
			assertThat((String) harness.single(inferredGlycolysis, ReactomeJavaConstants.summation)
				.getAttribute(ReactomeJavaConstants.text), startsWith("This event has been computationally inferred"));
			assertThat(inferredGlycolysis.getAttribute(ReactomeJavaConstants.releaseDate), is(RELEASE_DATE));
			assertThat(harness.stableIdentifier(inferredGlycolysis),
				is("R-" + TARGET_ABBREVIATION + "-" + glycolysis.getDbId()));
			// Linked to the human pathway by inferredFrom only. In the graph model that is the same inferredTo
			// relationship as the human pathway's orthologousEvent, so the human pathway is not committed again to
			// store it; and an orthologousEvent on the inferred pathway would be a second, inferred-to-human one.
			assertThat(dbIds(harness.values(inferredGlycolysis, ReactomeJavaConstants.inferredFrom)),
				contains(glycolysis.getDbId()));
			assertThat(harness.values(inferredGlycolysis, ReactomeJavaConstants.orthologousEvent), is(empty()));
			assertThat(dbIds(api.getCommits()), not(hasItem(glycolysis.getDbId())));

			assertThat(hasEventDbIds(inferredGlycolysis), contains(
				inferred(ReactomeJavaConstants.Reaction, r1).getDbId(),
				inferred(ReactomeJavaConstants.Reaction, r2).getDbId()));
			assertThat(hasEventDbIds(inferred(TOP_LEVEL_PATHWAY, metabolism)), contains(inferredGlycolysis.getDbId()));
		}

		@Test
		void inferredPathwayHoldsOnlyTheReactionsThatWereInferred() throws Exception {
			SimpleInstance r1 = mappedReaction("R1");
			SimpleInstance glycolysis = fx.pathway(ReactomeJavaConstants.Pathway, "Glycolysis", r1, unmappedReaction("R3"));
			fx.pathway(TOP_LEVEL_PATHWAY, "Metabolism", glycolysis);

			runInference();

			assertThat(hasEventDbIds(inferred(ReactomeJavaConstants.Pathway, glycolysis)),
				contains(inferred(ReactomeJavaConstants.Reaction, r1).getDbId()));
		}

		@Test
		void pathwayWithNoInferredReactionIsNotInferred() throws Exception {
			SimpleInstance glycolysis = fx.pathway(ReactomeJavaConstants.Pathway, "Glycolysis", mappedReaction("R1"));
			SimpleInstance fructose = fx.pathway(ReactomeJavaConstants.Pathway, "Fructose", unmappedReaction("R3"));
			SimpleInstance metabolism = fx.pathway(TOP_LEVEL_PATHWAY, "Metabolism", glycolysis, fructose);

			runInference();

			assertThat(hasInferred(ReactomeJavaConstants.Pathway, fructose), is(false));
			assertThat(hasEventDbIds(inferred(TOP_LEVEL_PATHWAY, metabolism)),
				contains(inferred(ReactomeJavaConstants.Pathway, glycolysis).getDbId()));
		}

		@Test
		void reactionInTwoPathwaysIsAddedToBoth() throws Exception {
			SimpleInstance r1 = mappedReaction("R1");
			SimpleInstance glycolysis = fx.pathway(ReactomeJavaConstants.Pathway, "Glycolysis", r1);
			SimpleInstance gluconeogenesis = fx.pathway(ReactomeJavaConstants.Pathway, "Gluconeogenesis", r1);
			SimpleInstance metabolism = fx.pathway(TOP_LEVEL_PATHWAY, "Metabolism", glycolysis, gluconeogenesis);

			runInference();

			long inferredR1 = inferred(ReactomeJavaConstants.Reaction, r1).getDbId();
			SimpleInstance inferredGlycolysis = inferred(ReactomeJavaConstants.Pathway, glycolysis);
			SimpleInstance inferredGluconeogenesis = inferred(ReactomeJavaConstants.Pathway, gluconeogenesis);
			assertThat(hasEventDbIds(inferredGlycolysis), contains(inferredR1));
			assertThat(hasEventDbIds(inferredGluconeogenesis), contains(inferredR1));
			// The shared parent is inferred once, holding both.
			assertThat(api.committedInstancesOfClass(TOP_LEVEL_PATHWAY), hasSize(1));
			assertThat(hasEventDbIds(inferred(TOP_LEVEL_PATHWAY, metabolism)),
				contains(inferredGlycolysis.getDbId(), inferredGluconeogenesis.getDbId()));
		}

		@Test
		void reactionDirectlyUnderATopLevelPathwayIsPlacedInItsInferredCopy() throws Exception {
			SimpleInstance r1 = mappedReaction("R1");
			SimpleInstance metabolism = fx.pathway(TOP_LEVEL_PATHWAY, "Metabolism", r1);

			runInference();

			SimpleInstance inferredMetabolism = inferred(TOP_LEVEL_PATHWAY, metabolism);
			assertThat(inferredMetabolism.getSchemaClassName(), is(TOP_LEVEL_PATHWAY));
			assertThat(hasEventDbIds(inferredMetabolism), contains(inferred(ReactomeJavaConstants.Reaction, r1).getDbId()));
		}

		@Test
		void everyLevelOfANestedHierarchyIsInferred() throws Exception {
			SimpleInstance r1 = mappedReaction("R1");
			SimpleInstance glycolysis = fx.pathway(ReactomeJavaConstants.Pathway, "Glycolysis", r1);
			SimpleInstance carbohydrates = fx.pathway(ReactomeJavaConstants.Pathway, "Carbohydrate metabolism", glycolysis);
			SimpleInstance metabolism = fx.pathway(TOP_LEVEL_PATHWAY, "Metabolism", carbohydrates);

			runInference();

			SimpleInstance inferredGlycolysis = inferred(ReactomeJavaConstants.Pathway, glycolysis);
			SimpleInstance inferredCarbohydrates = inferred(ReactomeJavaConstants.Pathway, carbohydrates);
			assertThat(hasEventDbIds(inferredGlycolysis), contains(inferred(ReactomeJavaConstants.Reaction, r1).getDbId()));
			assertThat(hasEventDbIds(inferredCarbohydrates), contains(inferredGlycolysis.getDbId()));
			assertThat(hasEventDbIds(inferred(TOP_LEVEL_PATHWAY, metabolism)), contains(inferredCarbohydrates.getDbId()));
		}

		@Test
		void pathwaysOnlyUnderDiseaseAreNotInferred() throws Exception {
			SimpleInstance r1 = mappedReaction("R1");
			SimpleInstance glycolysis = fx.pathway(ReactomeJavaConstants.Pathway, "Glycolysis", r1);
			fx.pathway(TOP_LEVEL_PATHWAY, "Metabolism", glycolysis);
			SimpleInstance infectiousDisease = fx.pathway(ReactomeJavaConstants.Pathway, "Infectious disease", r1);
			SimpleInstance disease = fx.pathway(DISEASE_PATHWAY_DB_ID, "Disease", infectiousDisease);

			runInference();

			// The reaction is also under Metabolism, so it is inferred -- but only into the Metabolism branch.
			assertThat(hasEventDbIds(inferred(ReactomeJavaConstants.Pathway, glycolysis)),
				contains(inferred(ReactomeJavaConstants.Reaction, r1).getDbId()));
			assertThat(hasInferred(ReactomeJavaConstants.Pathway, infectiousDisease), is(false));
			assertThat(hasInferred(ReactomeJavaConstants.Pathway, disease), is(false));
		}
	}

	@Nested
	@DisplayName("Preceding events")
	class PrecedingEvents {

		@Test
		void precedingReactionIsLinkedBetweenTheInferredReactions() throws Exception {
			SimpleInstance r1 = mappedReaction("R1");
			SimpleInstance r2 = mappedReaction("R2");
			fx.set(r2, ReactomeJavaConstants.precedingEvent, r1);
			fx.pathway(TOP_LEVEL_PATHWAY, "Metabolism", r1, r2);

			runInference();

			assertThat(dbIds(harness.values(inferred(ReactomeJavaConstants.Reaction, r2), ReactomeJavaConstants.precedingEvent)),
				contains(inferred(ReactomeJavaConstants.Reaction, r1).getDbId()));
		}

		@Test
		void precedingReactionThatWasNotInferredIsLeftOut() throws Exception {
			SimpleInstance r3 = unmappedReaction("R3");
			SimpleInstance r2 = mappedReaction("R2");
			fx.set(r2, ReactomeJavaConstants.precedingEvent, r3);
			fx.pathway(TOP_LEVEL_PATHWAY, "Metabolism", r3, r2);

			runInference();

			assertThat(harness.values(inferred(ReactomeJavaConstants.Reaction, r2), ReactomeJavaConstants.precedingEvent),
				is(empty()));
		}

		@Test
		void precedingPathwayIsLinkedBetweenTheInferredPathways() throws Exception {
			SimpleInstance glycolysis = fx.pathway(ReactomeJavaConstants.Pathway, "Glycolysis", mappedReaction("R1"));
			SimpleInstance citricAcidCycle = fx.pathway(ReactomeJavaConstants.Pathway, "TCA cycle", mappedReaction("R2"));
			fx.set(citricAcidCycle, ReactomeJavaConstants.precedingEvent, glycolysis);
			fx.pathway(TOP_LEVEL_PATHWAY, "Metabolism", glycolysis, citricAcidCycle);

			runInference();

			assertThat(dbIds(harness.values(inferred(ReactomeJavaConstants.Pathway, citricAcidCycle),
					ReactomeJavaConstants.precedingEvent)),
				contains(inferred(ReactomeJavaConstants.Pathway, glycolysis).getDbId()));
		}
	}

	@Nested
	@DisplayName("Reactions already inferred for the species")
	class AlreadyInferred {

		@Test
		void electronicallyInferredCounterpartIsReusedAndPlacedInThePathway() throws Exception {
			SimpleInstance counterpart = electronicCounterpart(harness.targetSpecies);
			SimpleInstance r1 = mappedReaction("R1");
			fx.set(r1, ReactomeJavaConstants.orthologousEvent, counterpart);
			SimpleInstance glycolysis = fx.pathway(ReactomeJavaConstants.Pathway, "Glycolysis", r1);

			runInference();

			assertThat(api.committedInstancesOfClass(ReactomeJavaConstants.Reaction), is(empty()));
			assertThat(hasEventDbIds(inferred(ReactomeJavaConstants.Pathway, glycolysis)), contains(counterpart.getDbId()));
			assertThat(harness.reactionInferrer.getEligibleCount(), is(0));
		}

		@Test
		void counterpartIsAlsoFoundThroughInferredFrom() throws Exception {
			SimpleInstance counterpart = electronicCounterpart(harness.targetSpecies);
			SimpleInstance r1 = mappedReaction("R1");
			fx.set(r1, ReactomeJavaConstants.inferredFrom, counterpart);
			SimpleInstance glycolysis = fx.pathway(ReactomeJavaConstants.Pathway, "Glycolysis", r1);

			runInference();

			assertThat(api.committedInstancesOfClass(ReactomeJavaConstants.Reaction), is(empty()));
			assertThat(hasEventDbIds(inferred(ReactomeJavaConstants.Pathway, glycolysis)), contains(counterpart.getDbId()));
		}

		/**
		 * Pins existing behaviour, as on main: a reaction already inferred by hand is neither inferred again nor
		 * placed in the inferred pathways (EventsInferrer records it in manualHumanEvents, which nothing reads), so a
		 * pathway holding only such reactions gets no inferred copy.
		 */
		@Test
		void manuallyInferredCounterpartIsNeitherReinferredNorPlacedInPathways() throws Exception {
			SimpleInstance counterpart = fx.existingReaction("R1", harness.targetSpecies);
			SimpleInstance r1 = mappedReaction("R1");
			fx.set(r1, ReactomeJavaConstants.orthologousEvent, counterpart);
			fx.pathway(TOP_LEVEL_PATHWAY, "Metabolism", r1);

			runInference();

			assertThat(api.committedInstancesOfClass(ReactomeJavaConstants.Reaction), is(empty()));
			assertThat(api.committedInstancesOfClass(TOP_LEVEL_PATHWAY), is(empty()));
		}

		@Test
		void diseaseCounterpartMeansTheReactionIsSkipped() throws Exception {
			SimpleInstance counterpart = electronicCounterpart(harness.targetSpecies);
			fx.set(counterpart, ReactomeJavaConstants.disease, fx.disease("cancer"));
			SimpleInstance r1 = mappedReaction("R1");
			fx.set(r1, ReactomeJavaConstants.orthologousEvent, counterpart);
			fx.pathway(TOP_LEVEL_PATHWAY, "Metabolism", r1);

			runInference();

			assertThat(api.committedInstancesOfClass(ReactomeJavaConstants.Reaction), is(empty()));
			assertThat(api.committedInstancesOfClass(TOP_LEVEL_PATHWAY), is(empty()));
		}

		@Test
		void chimericCounterpartIsIgnored() throws Exception {
			SimpleInstance counterpart = electronicCounterpart(harness.targetSpecies);
			fx.setSingle(counterpart, ReactomeJavaConstants.isChimeric, true);
			SimpleInstance r1 = mappedReaction("R1");
			fx.set(r1, ReactomeJavaConstants.orthologousEvent, counterpart);

			runInference();

			assertThat(hasInferred(ReactomeJavaConstants.Reaction, r1), is(true));
		}

		@Test
		void counterpartOfAnotherSpeciesIsIgnored() throws Exception {
			SimpleInstance counterpart = electronicCounterpart(fx.species("Mus musculus"));
			SimpleInstance r1 = mappedReaction("R1");
			fx.set(r1, ReactomeJavaConstants.orthologousEvent, counterpart);

			runInference();

			assertThat(hasInferred(ReactomeJavaConstants.Reaction, r1), is(true));
		}

		private SimpleInstance electronicCounterpart(SimpleInstance species) throws Exception {
			SimpleInstance counterpart = fx.existingReaction("R1", species);
			return fx.setSingle(counterpart, ReactomeJavaConstants.evidenceType, harness.utils.getEvidenceType());
		}
	}

	@Nested
	@DisplayName("The run as a whole")
	class WholeRun {

		@Test
		void reportGivesTheInferredOutOfTheEligibleReactions() throws Exception {
			mappedReaction("R1");
			mappedReaction("R2");
			unmappedReaction("R3");
			fx.reaction("ATP => ADP").input(fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build(); // no proteins

			runInference();

			assertThat(Files.readAllLines(reportFile()),
				contains("hsap to " + TARGET_SPECIES_CODE + ":\t2 out of 3 eligible reactions (66.67%)"));
		}

		@Test
		void onlyHumanReactionsAreInferred() throws Exception {
			SimpleInstance r1 = mappedReaction("R1");
			SimpleInstance mouseReaction = fx.existingReaction("mouse reaction", fx.species("Mus musculus"));

			runInference();

			assertThat(hasInferred(ReactomeJavaConstants.Reaction, r1), is(true));
			assertThat(hasInferred(ReactomeJavaConstants.Reaction, mouseReaction), is(false));
			assertThat(harness.reactionInferrer.getEligibleCount(), is(1));
		}
	}

	// ---- Helpers ------------------------------------------------------------------------------------------------

	private static final String TOP_LEVEL_PATHWAY = "TopLevelPathway";

	private void runInference() throws Exception {
		harness.eventsInferrer.inferEvents();
	}

	// PROT1 => PROT2, which has homologues.
	private SimpleInstance mappedReaction(String name) {
		return fx.reaction(name).input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002")).build();
	}

	// PROT99 => PROT1, which cannot be inferred as PROT99 has no homologue.
	private SimpleInstance unmappedReaction(String name) {
		return fx.reaction(name).input(fx.ewas("PROT99", "P00099")).output(fx.ewas("PROT1", "P00001")).build();
	}

	// The instance of the class this run created from the human event, as stored.
	private SimpleInstance inferred(String schemaClass, SimpleInstance humanEvent) {
		List<SimpleInstance> inferredInstances = inferredInstances(schemaClass, humanEvent);
		assertThat("inferred " + schemaClass + " of " + humanEvent, inferredInstances, hasSize(1));
		return inferredInstances.get(0);
	}

	private boolean hasInferred(String schemaClass, SimpleInstance humanEvent) {
		return !inferredInstances(schemaClass, humanEvent).isEmpty();
	}

	private List<SimpleInstance> inferredInstances(String schemaClass, SimpleInstance humanEvent) {
		return api.committedInstancesOfClass(schemaClass).stream()
			.filter(instance -> dbIds(harness.values(instance, ReactomeJavaConstants.inferredFrom))
				.contains(humanEvent.getDbId()))
			.collect(Collectors.toList());
	}

	private List<Long> hasEventDbIds(SimpleInstance pathway) {
		return dbIds(harness.values(pathway, ReactomeJavaConstants.hasEvent));
	}
}
