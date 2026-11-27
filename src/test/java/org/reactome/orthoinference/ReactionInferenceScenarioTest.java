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
import java.util.*;
import java.util.stream.Collectors;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import static org.reactome.orthoinference.InferenceTestHarness.*;

/**
 * Runs representative human reactions through the real inference classes -- ReactionInferrer,
 * OrthologousEntityGenerator, EWASInferrer, ProteinCountUtility, SkipInstanceChecker and the rest -- against an
 * in-memory curator-tool-ws ({@link InMemoryCuratorToolAPI}) and a small orthopairs fixture, and checks what is
 * stored for the target species. The homology fixture is described in {@link InferenceTestHarness}.
 */
class ReactionInferenceScenarioTest {

	@TempDir
	Path tempDir;

	private InferenceTestHarness harness;
	private InMemoryCuratorToolAPI api;
	private HumanFixtures fx;
	private ReactionInferrer reactionInferrer;
	private SimpleInstance targetSpecies;

	@BeforeEach
	void setUp() throws Exception {
		wire(null);
	}

	/**
	 * Builds the inference classes afresh against a new, empty in-memory database.
	 *
	 * @param alternateReferenceDatabaseJson - the species' "alt_refdb" entry in Species.json, or null for none.
	 */
	private void wire(String alternateReferenceDatabaseJson) throws Exception {
		harness = new InferenceTestHarness(tempDir, alternateReferenceDatabaseJson);
		api = harness.api;
		fx = harness.fx;
		reactionInferrer = harness.reactionInferrer;
		targetSpecies = harness.targetSpecies;
	}

	@AfterEach
	void deleteOutputFiles() throws Exception {
		InferenceTestHarness.deleteOutputFiles();
	}

	@Nested
	@DisplayName("Reactions of proteins")
	class Proteins {

		@Test
		void ewasToEwasReactionIsInferredWithItsAnnotations() throws Exception {
			SimpleInstance prot1 = fx.ewas("PROT1", "P00001");
			SimpleInstance prot2 = fx.ewas("PROT2", "P00002");
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2").input(prot1).output(prot2).build();

			SimpleInstance inferred = infer(reaction);

			assertThat(inferred, is(notNullValue()));
			assertThat(inferred.getSchemaClassName(), is(ReactomeJavaConstants.Reaction));
			assertThat(inferred.getDisplayName(), is("PROT1 => PROT2"));
			assertThat(inferred.getAttribute(ReactomeJavaConstants.name), is(Collections.singletonList("PROT1 => PROT2")));
			assertThat(dbIds(values(inferred, ReactomeJavaConstants.species)), contains(targetSpecies.getDbId()));
			assertThat(dbIds(values(inferred, ReactomeJavaConstants.compartment)), contains(fx.cytosol.getDbId()));
			assertThat(single(inferred, ReactomeJavaConstants.evidenceType).getDisplayName(),
				is("inferred by electronic annotation"));
			assertThat((String) single(inferred, ReactomeJavaConstants.summation).getAttribute(ReactomeJavaConstants.text),
				startsWith("This event has been computationally inferred"));
			assertThat(inferred.getAttribute(ReactomeJavaConstants.releaseDate), is(RELEASE_DATE));
			assertThat(single(inferred, ReactomeJavaConstants.stableIdentifier).getAttribute(ReactomeJavaConstants.identifier),
				is("R-" + TARGET_ABBREVIATION + "-" + reaction.getDbId()));

			// Linked to its source by inferredFrom, which in the graph model is the same inferredTo relationship as
			// the source's orthologousEvent. The inferred reaction has no orthologousEvent of its own: that would be a
			// second, inferred-to-human inferredTo.
			assertThat(dbIds(values(inferred, ReactomeJavaConstants.inferredFrom)), contains(reaction.getDbId()));
			assertThat(values(inferred, ReactomeJavaConstants.orthologousEvent), is(empty()));
			assertThat(dbIds(values(reaction, ReactomeJavaConstants.orthologousEvent)), contains(inferred.getDbId()));

			assertInferredEwas(single(inferred, ReactomeJavaConstants.input), prot1, "Q00001", "Tst1");
			assertInferredEwas(single(inferred, ReactomeJavaConstants.output), prot2, "Q00002", "Tst2");

			assertThat(reactionInferrer.getEligibleCount(), is(1));
			assertThat(reactionInferrer.getInferredCount(), is(1));
			assertThat(Files.readAllLines(inferredFile()),
				contains(inferred.getDbId() + "\t" + inferred.getDisplayName()));
		}

		@Test
		void homologueWithoutGeneNameIsNamedAfterItsIdentifier() throws Exception {
			SimpleInstance prot3 = fx.ewas("PROT3", "P00003");
			SimpleInstance reaction = fx.reaction("PROT3 => PROT1").input(prot3).output(fx.ewas("PROT1", "P00001")).build();

			SimpleInstance inferredInput = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(inferredInput.getAttribute(ReactomeJavaConstants.name), is(Collections.singletonList("Q00003")));
			assertThat(inferredInput.getDisplayName(), is("Q00003 [cytosol]"));
		}

		@Test
		void ewasWithTwoHomologuesBecomesDefinedSetOfBoth() throws Exception {
			SimpleInstance prot10 = fx.ewas("PROT10", "P00010");
			SimpleInstance reaction = fx.reaction("PROT10 => PROT1").input(prot10).output(fx.ewas("PROT1", "P00001")).build();

			SimpleInstance inferredInput = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(inferredInput.getSchemaClassName(), is(ReactomeJavaConstants.DefinedSet));
			assertThat(inferredInput.getDisplayName(), is("Homologues of PROT10 [cytosol]"));
			assertThat(dbIds(values(inferredInput, ReactomeJavaConstants.species)), contains(targetSpecies.getDbId()));
			assertThat(dbIds(values(inferredInput, ReactomeJavaConstants.inferredFrom)), contains(prot10.getDbId()));

			List<SimpleInstance> members = inflated(inferredInput, ReactomeJavaConstants.hasMember);
			assertThat(members.stream().map(ReactionInferenceScenarioTest.this::referenceIdentifier).collect(Collectors.toList()),
				contains("Q00010", "Q00011"));
			// Both homologues come from the same human EWAS, so the second is given a paralog suffix.
			assertThat(members.stream().map(ReactionInferenceScenarioTest.this::stableIdentifier).collect(Collectors.toList()),
				contains("R-TST-" + prot10.getDbId(), "R-TST-" + prot10.getDbId() + "-2"));
		}

		@Test
		void phosphorylatedEwasIsRenamedAndKeepsItsResidue() throws Exception {
			SimpleInstance residue = fx.phosphoSerine(15);
			SimpleInstance phosphoProt1 = fx.ewas("p-S15-PROT1", "P00001", residue);
			SimpleInstance reaction = fx.reaction("PROT2 phosphorylates PROT1")
				.input(fx.ewas("PROT1", "P00001")).output(phosphoProt1).build();

			SimpleInstance inferredOutput = single(infer(reaction), ReactomeJavaConstants.output);

			assertThat(inferredOutput.getAttribute(ReactomeJavaConstants.name), is(Arrays.asList("phospho-Tst1", "Q00001")));
			assertThat(inferredOutput.getDisplayName(), is("phospho-Tst1 [cytosol]"));

			SimpleInstance inferredResidue = single(inferredOutput, ReactomeJavaConstants.hasModifiedResidue);
			assertThat(inferredResidue.getAttribute(ReactomeJavaConstants.coordinate), is(15));
			assertThat(inferredResidue.getDisplayName(), is("O-phospho-L-serine at 15 (in Homo sapiens)"));
			assertThat(single(inferredResidue, ReactomeJavaConstants.referenceSequence).getDbId(),
				is(single(inferredOutput, ReactomeJavaConstants.referenceEntity).getDbId()));
			assertThat(single(inferredResidue, ReactomeJavaConstants.psiMod).getDbId(),
				is(single(residue, ReactomeJavaConstants.psiMod).getDbId()));
		}

		@Test
		void unmappedInputStopsInference() throws Exception {
			SimpleInstance reaction = fx.reaction("PROT99 => PROT1")
				.input(fx.ewas("PROT99", "P00099")).output(fx.ewas("PROT1", "P00001")).build();

			assertNotInferred(reaction, true);
		}

		@Test
		void unmappedOutputStopsInference() throws Exception {
			SimpleInstance reaction = fx.reaction("PROT1 => PROT99")
				.input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT99", "P00099")).build();

			assertNotInferred(reaction, true);
		}
	}

	@Nested
	@DisplayName("Small molecules")
	class SmallMolecules {

		@Test
		void simpleEntityIsUsedAsItIs() throws Exception {
			SimpleInstance atp = fx.simpleEntity("ATP");
			SimpleInstance adp = fx.simpleEntity("ADP");
			SimpleInstance reaction = fx.reaction("PROT1 hydrolyses ATP")
				.input(atp, fx.ewas("PROT1", "P00001")).output(adp, fx.ewas("PROT2", "P00002")).build();

			SimpleInstance inferred = infer(reaction);

			assertThat(dbIds(values(inferred, ReactomeJavaConstants.input)), hasItem(atp.getDbId()));
			assertThat(dbIds(values(inferred, ReactomeJavaConstants.output)), hasItem(adp.getDbId()));
			assertThat(values(inferred, ReactomeJavaConstants.input), hasSize(2));
		}

		@Test
		void reactionWithoutProteinsIsNotEligible() throws Exception {
			SimpleInstance reaction = fx.reaction("ATP => ADP")
				.input(fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build();

			assertNotInferred(reaction, false);
		}
	}

	@Nested
	@DisplayName("Complexes")
	class Complexes {

		@Test
		void complexWithThreeOfFourProteinsMappedIsInferredWithGhostForTheFourth() throws Exception {
			SimpleInstance prot99 = fx.ewas("PROT99", "P00099");
			SimpleInstance complex = fx.complex("PROT1:PROT2:PROT3:PROT99",
				fx.ewas("PROT1", "P00001"), fx.ewas("PROT2", "P00002"), fx.ewas("PROT3", "P00003"), prot99);
			SimpleInstance reaction = fx.reaction("complex binds ATP")
				.input(complex, fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build();

			SimpleInstance inferredComplex = inflated(infer(reaction), ReactomeJavaConstants.input).stream()
				.filter(input -> input.getSchemaClassName().equals(ReactomeJavaConstants.Complex))
				.findFirst().orElseThrow(AssertionError::new);

			assertThat(inferredComplex.getDisplayName(), is(complex.getDisplayName()));
			assertThat(dbIds(values(inferredComplex, ReactomeJavaConstants.species)), contains(targetSpecies.getDbId()));
			assertThat(dbIds(values(inferredComplex, ReactomeJavaConstants.inferredFrom)), contains(complex.getDbId()));
			assertThat((String) single(inferredComplex, ReactomeJavaConstants.summation).getAttribute(ReactomeJavaConstants.text),
				startsWith("This complex/polymer has been computationally inferred"));

			List<SimpleInstance> components = inflated(inferredComplex, ReactomeJavaConstants.hasComponent);
			assertThat(components, hasSize(4));
			assertThat(components.subList(0, 3).stream().map(ReactionInferenceScenarioTest.this::referenceIdentifier).collect(Collectors.toList()),
				contains("Q00001", "Q00002", "Q00003"));

			SimpleInstance ghost = components.get(3);
			assertThat(ghost.getSchemaClassName(), is(ReactomeJavaConstants.GenomeEncodedEntity));
			assertThat(ghost.getDisplayName(), is("Ghost homologue of " + prot99.getDisplayName()));
			assertThat(dbIds(values(ghost, ReactomeJavaConstants.species)), contains(targetSpecies.getDbId()));
			assertThat(dbIds(values(ghost, ReactomeJavaConstants.inferredFrom)), contains(prot99.getDbId()));
		}

		@Test
		void complexWithHalfItsProteinsMappedStopsInference() throws Exception {
			SimpleInstance complex = fx.complex("PROT1:PROT2:PROT98:PROT99",
				fx.ewas("PROT1", "P00001"), fx.ewas("PROT2", "P00002"),
				fx.ewas("PROT98", "P00098"), fx.ewas("PROT99", "P00099"));
			SimpleInstance reaction = fx.reaction("complex binds ATP")
				.input(complex, fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build();

			assertNotInferred(reaction, true);
		}
	}

	@Nested
	@DisplayName("Entity sets")
	class EntitySets {

		@Test
		void definedSetKeepsOnlyItsInferredMembers() throws Exception {
			SimpleInstance set = fx.definedSet("PROT1,2,99",
				fx.ewas("PROT1", "P00001"), fx.ewas("PROT2", "P00002"), fx.ewas("PROT99", "P00099"));
			SimpleInstance reaction = fx.reaction("set => PROT3").input(set).output(fx.ewas("PROT3", "P00003")).build();

			SimpleInstance inferredSet = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(inferredSet.getSchemaClassName(), is(ReactomeJavaConstants.DefinedSet));
			assertThat(inferredSet.getDisplayName(), is(set.getDisplayName()));
			assertThat(dbIds(values(inferredSet, ReactomeJavaConstants.species)), contains(targetSpecies.getDbId()));
			assertThat(inflated(inferredSet, ReactomeJavaConstants.hasMember).stream()
				.map(ReactionInferenceScenarioTest.this::referenceIdentifier).collect(Collectors.toList()), contains("Q00001", "Q00002"));
		}

		@Test
		void definedSetWithOneInferredMemberIsReplacedByThatMember() throws Exception {
			SimpleInstance prot1 = fx.ewas("PROT1", "P00001");
			SimpleInstance set = fx.definedSet("PROT1,99", prot1, fx.ewas("PROT99", "P00099"));
			SimpleInstance reaction = fx.reaction("set => PROT3").input(set).output(fx.ewas("PROT3", "P00003")).build();

			SimpleInstance inferredInput = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(inferredInput.getSchemaClassName(), is(ReactomeJavaConstants.EntityWithAccessionedSequence));
			assertThat(referenceIdentifier(inferredInput), is("Q00001"));
			assertThat(dbIds(values(inferredInput, ReactomeJavaConstants.inferredFrom)), hasItem(prot1.getDbId()));
		}

		@Test
		void definedSetWithNoInferredMemberStopsInference() throws Exception {
			SimpleInstance set = fx.definedSet("PROT98,99", fx.ewas("PROT98", "P00098"), fx.ewas("PROT99", "P00099"));
			SimpleInstance reaction = fx.reaction("set => PROT3").input(set).output(fx.ewas("PROT3", "P00003")).build();

			assertNotInferred(reaction, true);
		}

		@Test
		void candidateSetKeepsItsInferredMembersAndCandidates() throws Exception {
			SimpleInstance set = fx.candidateSet("PROT1 or PROT2",
				Collections.singletonList(fx.ewas("PROT1", "P00001")),
				Collections.singletonList(fx.ewas("PROT2", "P00002")));
			SimpleInstance reaction = fx.reaction("set => PROT3").input(set).output(fx.ewas("PROT3", "P00003")).build();

			SimpleInstance inferredSet = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(inferredSet.getSchemaClassName(), is(ReactomeJavaConstants.CandidateSet));
			assertThat(referenceIdentifier(single(inferredSet, ReactomeJavaConstants.hasMember)), is("Q00001"));
			assertThat(referenceIdentifier(single(inferredSet, ReactomeJavaConstants.hasCandidate)), is("Q00002"));
		}

		@Test
		void candidateSetWithNoInferredCandidateBecomesDefinedSetOfItsMembers() throws Exception {
			SimpleInstance set = fx.candidateSet("PROT1, PROT2 or PROT99",
				Arrays.asList(fx.ewas("PROT1", "P00001"), fx.ewas("PROT2", "P00002")),
				Collections.singletonList(fx.ewas("PROT99", "P00099")));
			SimpleInstance reaction = fx.reaction("set => PROT3").input(set).output(fx.ewas("PROT3", "P00003")).build();

			SimpleInstance inferredSet = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(inferredSet.getSchemaClassName(), is(ReactomeJavaConstants.DefinedSet));
			assertThat(dbIds(values(inferredSet, ReactomeJavaConstants.species)), contains(targetSpecies.getDbId()));
			assertThat(inflated(inferredSet, ReactomeJavaConstants.hasMember).stream()
				.map(ReactionInferenceScenarioTest.this::referenceIdentifier).collect(Collectors.toList()), contains("Q00001", "Q00002"));
		}
	}

	@Nested
	@DisplayName("Catalysts")
	class Catalysts {

		@Test
		void catalystActivityIsInferredWithItsEntityAndActiveUnit() throws Exception {
			SimpleInstance enzyme = fx.ewas("PROT3", "P00003");
			SimpleInstance catalyst = fx.catalystActivity(enzyme, enzyme);
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2")
				.input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002")).catalyst(catalyst).build();

			SimpleInstance inferredCatalyst = single(infer(reaction), ReactomeJavaConstants.catalystActivity);

			assertThat(inferredCatalyst.getSchemaClassName(), is(ReactomeJavaConstants.CatalystActivity));
			assertThat(inferredCatalyst.getDisplayName(), is(catalyst.getDisplayName()));
			// The GO molecular function is species-independent, so the human one is used.
			assertThat(single(inferredCatalyst, ReactomeJavaConstants.activity).getDbId(), is(fx.catalyticActivity.getDbId()));
			SimpleInstance inferredEnzyme = single(inferredCatalyst, ReactomeJavaConstants.physicalEntity);
			assertInferredEwas(inferredEnzyme, enzyme, "Q00003", null);
			assertThat(dbIds(values(inferredCatalyst, ReactomeJavaConstants.activeUnit)), contains(inferredEnzyme.getDbId()));
		}

		@Test
		void everyCatalystActivityIsKept() throws Exception {
			SimpleInstance catalyst3 = fx.catalystActivity(fx.ewas("PROT3", "P00003"));
			SimpleInstance catalyst4 = fx.catalystActivity(fx.ewas("PROT4", "P00004"));
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2")
				.input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002"))
				.catalyst(catalyst3, catalyst4).build();

			SimpleInstance inferred = infer(reaction);

			assertThat(inflated(inferred, ReactomeJavaConstants.catalystActivity).stream()
					.map(inferredCatalyst -> referenceIdentifier(single(inferredCatalyst, ReactomeJavaConstants.physicalEntity)))
					.collect(Collectors.toList()),
				contains("Q00003", "Q00004"));
		}

		@Test
		void unmappedCatalystStopsInference() throws Exception {
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2")
				.input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002"))
				.catalyst(fx.catalystActivity(fx.ewas("PROT99", "P00099"))).build();

			assertNotInferred(reaction, true);
		}
	}

	@Nested
	@DisplayName("Regulations")
	class Regulations {

		@Test
		void regulationWithInferredRegulatorIsAdded() throws Exception {
			SimpleInstance regulator = fx.ewas("PROT4", "P00004");
			SimpleInstance regulation = fx.regulation(ReactomeJavaConstants.PositiveRegulation, regulator);
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2")
				.input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002")).regulatedBy(regulation).build();

			SimpleInstance inferredRegulation = single(infer(reaction), ReactomeJavaConstants.regulatedBy);

			assertThat(inferredRegulation.getSchemaClassName(), is(ReactomeJavaConstants.PositiveRegulation));
			assertThat(inferredRegulation.getDisplayName(), is(regulation.getDisplayName()));
			assertInferredEwas(single(inferredRegulation, ReactomeJavaConstants.regulator), regulator, "Q00004", null);
		}

		@Test
		void regulationWithUnmappedRegulatorIsDropped() throws Exception {
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2")
				.input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002"))
				.regulatedBy(fx.regulation(ReactomeJavaConstants.NegativeRegulation, fx.ewas("PROT99", "P00099")))
				.build();

			SimpleInstance inferred = infer(reaction);

			assertThat(inferred, is(notNullValue()));
			assertThat(values(inferred, ReactomeJavaConstants.regulatedBy), is(empty()));
		}

		@Test
		void requirementWithUnmappedRegulatorStopsInference() throws Exception {
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2")
				.input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002"))
				.regulatedBy(fx.regulation(ReactomeJavaConstants.Requirement, fx.ewas("PROT99", "P00099")))
				.build();

			assertNotInferred(reaction, true);
		}
	}

	@Nested
	@DisplayName("Reactions that are skipped")
	class Skipped {

		@Test
		void chimericReactionIsSkipped() throws Exception {
			assertNotInferred(proteinReaction().with(ReactomeJavaConstants.isChimeric, true).build(), false);
		}

		@Test
		void diseaseReactionIsSkipped() throws Exception {
			SimpleInstance reaction = proteinReaction()
				.with(ReactomeJavaConstants.disease, new ArrayList<>(Collections.singletonList(fx.disease("cancer"))))
				.build();

			assertNotInferred(reaction, false);
		}

		@Test
		void manuallyInferredReactionIsSkipped() throws Exception {
			SimpleInstance otherReaction = proteinReaction().build();
			SimpleInstance reaction = proteinReaction()
				.with(ReactomeJavaConstants.inferredFrom, new ArrayList<>(Collections.singletonList(otherReaction)))
				.build();

			assertNotInferred(reaction, false);
		}

		@Test
		void reactionOnlyUnderDiseaseIsSkipped() throws Exception {
			SimpleInstance reaction = proteinReaction().build();
			fx.pathway(DISEASE_PATHWAY_DB_ID, "Disease", reaction);

			assertNotInferred(reaction, false);
		}

		@Test
		void reactionUnderDiseaseAndAnotherTopLevelPathwayIsInferred() throws Exception {
			SimpleInstance reaction = proteinReaction().build();
			fx.pathway(DISEASE_PATHWAY_DB_ID, "Disease", reaction);
			fx.pathway(1_430_728L, "Metabolism", reaction);

			assertThat(infer(reaction), is(notNullValue()));
		}

		private HumanFixtures.ReactionBuilder proteinReaction() {
			return fx.reaction("PROT1 => PROT2").input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002"));
		}
	}

	@Nested
	@DisplayName("Reuse across reactions")
	class Reuse {

		@Test
		void reactionIsInferredOnlyOnce() throws Exception {
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2")
				.input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002")).build();

			reactionInferrer.inferReaction(reaction);
			reactionInferrer.inferReaction(reaction);

			assertThat(api.committedInstancesOfClass(ReactomeJavaConstants.Reaction), hasSize(1));
			assertThat(reactionInferrer.getEligibleCount(), is(1));
			assertThat(reactionInferrer.getInferredCount(), is(1));
		}

		@Test
		void entitySharedByTwoReactionsIsInferredOnce() throws Exception {
			SimpleInstance prot2 = fx.ewas("PROT2", "P00002");
			SimpleInstance first = fx.reaction("PROT1 => PROT2").input(fx.ewas("PROT1", "P00001")).output(prot2).build();
			SimpleInstance second = fx.reaction("PROT2 => PROT3").input(prot2).output(fx.ewas("PROT3", "P00003")).build();

			SimpleInstance firstOutput = single(infer(first), ReactomeJavaConstants.output);
			SimpleInstance secondInput = single(infer(second), ReactomeJavaConstants.input);

			assertThat(secondInput.getDbId(), is(firstOutput.getDbId()));
			assertThat(api.committedInstancesOfClass(ReactomeJavaConstants.EntityWithAccessionedSequence).stream()
				.filter(ewas -> "Q00002".equals(referenceIdentifier(ewas))).count(), is(1L));
		}
	}

	@Nested
	@DisplayName("Instances already in the database")
	class ExistingInstances {

		@Test
		void previouslyInferredEwasIsReusedRatherThanCreated() throws Exception {
			SimpleInstance existingRgp =
				fx.existingReferenceGeneProduct(ReactomeJavaConstants.ReferenceGeneProduct, "Q00001", targetSpecies);
			SimpleInstance existingEwas = fx.existingEwas(existingRgp, targetSpecies, "R-TST-1");
			SimpleInstance prot1 = fx.ewas("PROT1", "P00001");
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2").input(prot1).output(fx.ewas("PROT2", "P00002")).build();

			SimpleInstance inferredInput = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(inferredInput.getDbId(), is(existingEwas.getDbId()));
			assertThat(stableIdentifier(inferredInput), is("R-TST-1"));
			// Now also inferred from this release's human EWAS.
			assertThat(dbIds(values(inferredInput, ReactomeJavaConstants.inferredFrom)), hasItem(prot1.getDbId()));
			assertThat(newInstancesFor(ReactomeJavaConstants.EntityWithAccessionedSequence, "Q00001"), is(0L));
			assertThat(newInstancesFor(ReactomeJavaConstants.ReferenceGeneProduct, "Q00001"), is(0L));
		}

		@Test
		void previouslyInferredComplexIsReusedWhenItsComponentsAre() throws Exception {
			SimpleInstance existingEwas1 = fx.existingEwas(
				fx.existingReferenceGeneProduct(ReactomeJavaConstants.ReferenceGeneProduct, "Q00001", targetSpecies),
				targetSpecies, "R-TST-1");
			SimpleInstance existingEwas2 = fx.existingEwas(
				fx.existingReferenceGeneProduct(ReactomeJavaConstants.ReferenceGeneProduct, "Q00002", targetSpecies),
				targetSpecies, "R-TST-2");
			SimpleInstance existingComplex =
				fx.existingComplex("PROT1:PROT2", targetSpecies, "R-TST-3", existingEwas1, existingEwas2);
			SimpleInstance complex = fx.complex("PROT1:PROT2", fx.ewas("PROT1", "P00001"), fx.ewas("PROT2", "P00002"));
			SimpleInstance reaction = fx.reaction("complex binds ATP")
				.input(complex, fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build();

			SimpleInstance inferredComplex = inputOfClass(infer(reaction), ReactomeJavaConstants.Complex);

			assertThat(inferredComplex.getDbId(), is(existingComplex.getDbId()));
			assertThat(api.committedInstancesOfClass(ReactomeJavaConstants.Complex), is(empty()));
		}

		@Test
		void storedInstanceOfExactlyTheInferredClassIsPreferredOverASubclass() throws Exception {
			// The isoform has the lower dbId, so only the class preference picks the ReferenceGeneProduct.
			fx.existingReferenceGeneProduct(ReactomeJavaConstants.ReferenceIsoform, "Q00001", targetSpecies);
			SimpleInstance existingRgp =
				fx.existingReferenceGeneProduct(ReactomeJavaConstants.ReferenceGeneProduct, "Q00001", targetSpecies);
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2")
				.input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002")).build();

			SimpleInstance inferredInput = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(single(inferredInput, ReactomeJavaConstants.referenceEntity).getDbId(), is(existingRgp.getDbId()));
		}

		@Test
		void lowestDbIdIsPreferredAmongEquallyGoodMatches() throws Exception {
			SimpleInstance firstRgp =
				fx.existingReferenceGeneProduct(ReactomeJavaConstants.ReferenceGeneProduct, "Q00001", targetSpecies);
			fx.existingReferenceGeneProduct(ReactomeJavaConstants.ReferenceGeneProduct, "Q00001", targetSpecies);
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2")
				.input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002")).build();

			SimpleInstance inferredInput = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(single(inferredInput, ReactomeJavaConstants.referenceEntity).getDbId(), is(firstRgp.getDbId()));
		}
	}

	@Nested
	@DisplayName("Protein counts of nested sets and complexes")
	class NestedProteinCounts {

		@Test
		void candidateSetOfOnlyCandidatesIsInferredWithItsMappedCandidates() throws Exception {
			SimpleInstance prot1 = fx.ewas("PROT1", "P00001");
			SimpleInstance prot99 = fx.ewas("PROT99", "P00099");
			SimpleInstance set = fx.candidateSet("PROT1 or PROT99", HumanFixtures.none(), Arrays.asList(prot1, prot99));
			SimpleInstance reaction = fx.reaction("set => PROT3").input(set).output(fx.ewas("PROT3", "P00003")).build();

			SimpleInstance inferredSet = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(inferredSet.getSchemaClassName(), is(ReactomeJavaConstants.CandidateSet));
			assertThat(referenceIdentifier(single(inferredSet, ReactomeJavaConstants.hasCandidate)), is("Q00001"));
			assertThat(values(inferredSet, ReactomeJavaConstants.hasMember), is(empty()));
		}

		@Test
		void candidateSetOfOnlyUnmappedCandidatesStopsInference() throws Exception {
			SimpleInstance set = fx.candidateSet("PROT98 or PROT99", HumanFixtures.none(),
				Arrays.asList(fx.ewas("PROT98", "P00098"), fx.ewas("PROT99", "P00099")));
			SimpleInstance reaction = fx.reaction("set => PROT3").input(set).output(fx.ewas("PROT3", "P00003")).build();

			assertNotInferred(reaction, true);
		}

		/**
		 * Pins existing behaviour rather than endorsing it. ProteinCountUtility.getCandidateProteinCounts emulates the
		 * Perl version: candidates are counted in dbId order, and if the first one counted has no homologue the whole
		 * CandidateSet counts as having no inferrable protein -- even though a later candidate is mapped. The same set
		 * with its candidates' dbIds the other way round is inferred (see
		 * candidateSetOfOnlyCandidatesIsInferredWithItsMappedCandidates). If this ordering quirk is fixed, this test
		 * should change with it.
		 */
		@Test
		void candidateSetIsRejectedWhenItsLowestDbIdCandidateIsUnmapped() throws Exception {
			SimpleInstance prot99 = fx.ewas("PROT99", "P00099"); // created first, so its reference has the lower dbId
			SimpleInstance prot1 = fx.ewas("PROT1", "P00001");
			SimpleInstance set = fx.candidateSet("PROT1 or PROT99", HumanFixtures.none(), Arrays.asList(prot1, prot99));
			SimpleInstance reaction = fx.reaction("set => PROT3").input(set).output(fx.ewas("PROT3", "P00003")).build();

			assertNotInferred(reaction, true);
		}

		@Test
		void definedSetOfComplexesKeepsTheComplexesThatPassTheThreshold() throws Exception {
			SimpleInstance complex12 = fx.complex("PROT1:PROT2", fx.ewas("PROT1", "P00001"), fx.ewas("PROT2", "P00002"));
			SimpleInstance complex34 = fx.complex("PROT3:PROT4", fx.ewas("PROT3", "P00003"), fx.ewas("PROT4", "P00004"));
			SimpleInstance complex9899 =
				fx.complex("PROT98:PROT99", fx.ewas("PROT98", "P00098"), fx.ewas("PROT99", "P00099"));
			SimpleInstance set = fx.definedSet("dimers", complex12, complex34, complex9899);
			SimpleInstance reaction = fx.reaction("dimers => PROT5").input(set).output(fx.ewas("PROT5", "P00005")).build();

			SimpleInstance inferredSet = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(inferredSet.getSchemaClassName(), is(ReactomeJavaConstants.DefinedSet));
			List<SimpleInstance> members = inflated(inferredSet, ReactomeJavaConstants.hasMember);
			assertThat(members.stream().map(SimpleInstance::getSchemaClassName).collect(Collectors.toList()),
				contains(ReactomeJavaConstants.Complex, ReactomeJavaConstants.Complex));
			assertThat(members.stream().map(member -> dbIds(values(member, ReactomeJavaConstants.inferredFrom)))
					.collect(Collectors.toList()),
				contains(contains(complex12.getDbId()), contains(complex34.getDbId())));
		}

		// Counted protein by protein, PROT1, PROT2 and PROT99 would be 2 of 3 (66%) and fail the 75% threshold; the
		// set counts as one protein, which is inferrable, so the complex is 2 of 2.
		@Test
		void setInsideAComplexCountsAsOneProteinTowardTheThreshold() throws Exception {
			SimpleInstance set = fx.definedSet("PROT2 or PROT99", fx.ewas("PROT2", "P00002"), fx.ewas("PROT99", "P00099"));
			SimpleInstance complex = fx.complex("PROT1:set", fx.ewas("PROT1", "P00001"), set);
			SimpleInstance reaction = fx.reaction("complex binds ATP")
				.input(complex, fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build();

			SimpleInstance inferredComplex = inputOfClass(infer(reaction), ReactomeJavaConstants.Complex);

			// The set inferred only PROT2's homologue, so that homologue stands in for it.
			assertThat(inflated(inferredComplex, ReactomeJavaConstants.hasComponent).stream()
				.map(ReactionInferenceScenarioTest.this::referenceIdentifier).collect(Collectors.toList()),
				contains("Q00001", "Q00002"));
		}
	}

	@Nested
	@DisplayName("Other kinds of physical entity")
	class OtherEntityKinds {

		@Test
		void polymerIsInferredWithItsRepeatedUnit() throws Exception {
			SimpleInstance polymer = fx.polymer("PROT1 polymer", fx.ewas("PROT1", "P00001"));
			SimpleInstance reaction = fx.reaction("polymer => PROT2").input(polymer).output(fx.ewas("PROT2", "P00002")).build();

			SimpleInstance inferredPolymer = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(inferredPolymer.getSchemaClassName(), is(ReactomeJavaConstants.Polymer));
			assertThat(inferredPolymer.getDisplayName(), is(polymer.getDisplayName()));
			assertThat(dbIds(values(inferredPolymer, ReactomeJavaConstants.inferredFrom)), contains(polymer.getDbId()));
			assertThat(referenceIdentifier(single(inferredPolymer, ReactomeJavaConstants.repeatedUnit)), is("Q00001"));
			assertThat((String) single(inferredPolymer, ReactomeJavaConstants.summation).getAttribute(ReactomeJavaConstants.text),
				startsWith("This complex/polymer has been computationally inferred"));
		}

		@Test
		void drugAndOtherEntityAreUsedAsTheyAre() throws Exception {
			SimpleInstance drug = fx.chemicalDrug("aspirin");
			SimpleInstance photon = fx.otherEntity("photon");
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2")
				.input(drug, photon, fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002")).build();

			SimpleInstance inferred = infer(reaction);

			assertThat(dbIds(values(inferred, ReactomeJavaConstants.input)), hasItems(drug.getDbId(), photon.getDbId()));
		}

		@Test
		void genomeEncodedEntityOnItsOwnStopsInference() throws Exception {
			SimpleInstance reaction = fx.reaction("PROTX => PROT2")
				.input(fx.genomeEncodedEntity("PROTX")).output(fx.ewas("PROT2", "P00002")).build();

			assertNotInferred(reaction, true);
		}

		@Test
		void genomeEncodedEntityInAComplexBecomesAGhost() throws Exception {
			SimpleInstance protX = fx.genomeEncodedEntity("PROTX");
			SimpleInstance complex = fx.complex("PROT1:PROT2:PROT3:PROTX",
				fx.ewas("PROT1", "P00001"), fx.ewas("PROT2", "P00002"), fx.ewas("PROT3", "P00003"), protX);
			SimpleInstance reaction = fx.reaction("complex binds ATP")
				.input(complex, fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build();

			List<SimpleInstance> components =
				inflated(inputOfClass(infer(reaction), ReactomeJavaConstants.Complex), ReactomeJavaConstants.hasComponent);

			assertGhostOf(components.get(3), protX);
		}

		/**
		 * Pins existing behaviour: a set is inferred only if it has a species of its own, so a set without one is used
		 * as it is -- the human set, in the inferred reaction -- even though its members could be inferred.
		 */
		@Test
		void setWithoutItsOwnSpeciesIsUsedAsItIs() throws Exception {
			SimpleInstance set = fx.definedSet("PROT1 or PROT2", fx.ewas("PROT1", "P00001"), fx.ewas("PROT2", "P00002"));
			set.setAttribute(ReactomeJavaConstants.species, null);
			api.add(set);
			SimpleInstance reaction = fx.reaction("set => PROT3").input(set).output(fx.ewas("PROT3", "P00003")).build();

			assertThat(single(infer(reaction), ReactomeJavaConstants.input).getDbId(), is(set.getDbId()));
		}

		@Test
		void complexOfSmallMoleculesIsUsedAsItIs() throws Exception {
			SimpleInstance atpMg = fx.complex("ATP:Mg2+", fx.simpleEntity("ATP"), fx.simpleEntity("Mg2+"));
			SimpleInstance reaction = fx.reaction("PROT1 binds ATP:Mg2+")
				.input(atpMg, fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002")).build();

			assertThat(dbIds(values(infer(reaction), ReactomeJavaConstants.input)), hasItem(atpMg.getDbId()));
		}
	}

	@Nested
	@DisplayName("Ghosts and override caches inside complexes")
	class InsideComplexes {

		@Test
		void candidateSetInferringNothingInsideAComplexBecomesAGhost() throws Exception {
			SimpleInstance set = fx.candidateSet("PROT99 candidates", HumanFixtures.none(),
				Collections.singletonList(fx.ewas("PROT99", "P00099")));
			SimpleInstance complex = fx.complex("PROT1:PROT2:PROT3:set",
				fx.ewas("PROT1", "P00001"), fx.ewas("PROT2", "P00002"), fx.ewas("PROT3", "P00003"), set);
			SimpleInstance reaction = fx.reaction("complex binds ATP")
				.input(complex, fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build();

			List<SimpleInstance> components =
				inflated(inputOfClass(infer(reaction), ReactomeJavaConstants.Complex), ReactomeJavaConstants.hasComponent);

			assertGhostOf(components.get(3), set);
		}

		@Test
		void definedSetInferringNothingInsideAComplexBecomesAGhost() throws Exception {
			SimpleInstance set = fx.definedSet("PROT99 set", fx.ewas("PROT99", "P00099"));
			SimpleInstance complex = fx.complex("PROT1:PROT2:PROT3:set",
				fx.ewas("PROT1", "P00001"), fx.ewas("PROT2", "P00002"), fx.ewas("PROT3", "P00003"), set);
			SimpleInstance reaction = fx.reaction("complex binds ATP")
				.input(complex, fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build();

			List<SimpleInstance> components =
				inflated(inputOfClass(infer(reaction), ReactomeJavaConstants.Complex), ReactomeJavaConstants.hasComponent);

			assertGhostOf(components.get(3), set);
		}

		@Test
		void complexListedTwiceInAComplexIsInferredOnce() throws Exception {
			SimpleInstance dimer = fx.complex("PROT1:PROT2", fx.ewas("PROT1", "P00001"), fx.ewas("PROT2", "P00002"));
			SimpleInstance tetramer = fx.complex("(PROT1:PROT2)2", dimer, dimer);
			SimpleInstance reaction = fx.reaction("tetramer binds ATP")
				.input(tetramer, fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build();

			List<SimpleInstance> components = values(inputOfClass(infer(reaction), ReactomeJavaConstants.Complex),
				ReactomeJavaConstants.hasComponent);

			assertThat(components, hasSize(2));
			assertThat(components.get(1).getDbId(), is(components.get(0).getDbId()));
			assertThat(api.committedInstancesOfClass(ReactomeJavaConstants.Complex), hasSize(2));
		}

		@Test
		void setListedTwiceInAComplexBecomesTheSameGhost() throws Exception {
			SimpleInstance set = fx.definedSet("PROT99 set", fx.ewas("PROT99", "P00099"));
			SimpleInstance complex = fx.complex("set:set:PROT1:PROT2:PROT3",
				set, set, fx.ewas("PROT1", "P00001"), fx.ewas("PROT2", "P00002"), fx.ewas("PROT3", "P00003"));
			SimpleInstance reaction = fx.reaction("complex binds ATP")
				.input(complex, fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build();

			List<SimpleInstance> components =
				inflated(inputOfClass(infer(reaction), ReactomeJavaConstants.Complex), ReactomeJavaConstants.hasComponent);

			assertGhostOf(components.get(0), set);
			assertThat(components.get(1).getDbId(), is(components.get(0).getDbId()));
		}

		// Inside a complex, a component is inferred regardless of the threshold (with a ghost for what has no
		// homologue). That result is cached apart, so the same complex met on its own is still held to the threshold.
		@Test
		void complexInferredInsideALargerComplexIsStillRejectedOnItsOwn() throws Exception {
			SimpleInstance prot99 = fx.ewas("PROT99", "P00099");
			SimpleInstance halfMapped = fx.complex("PROT1:PROT99", fx.ewas("PROT1", "P00001"), prot99);
			SimpleInstance large = fx.complex("large", halfMapped,
				fx.ewas("PROT2", "P00002"), fx.ewas("PROT3", "P00003"), fx.ewas("PROT4", "P00004"));
			SimpleInstance largeReaction = fx.reaction("large complex binds ATP")
				.input(large, fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build();
			SimpleInstance halfMappedReaction = fx.reaction("PROT1:PROT99 binds ATP")
				.input(halfMapped, fx.simpleEntity("ATP")).output(fx.simpleEntity("ADP")).build();

			SimpleInstance inferredLarge = inputOfClass(infer(largeReaction), ReactomeJavaConstants.Complex);
			SimpleInstance inferredHalfMapped = inflated(inferredLarge, ReactomeJavaConstants.hasComponent).get(0);
			List<SimpleInstance> halfMappedComponents = inflated(inferredHalfMapped, ReactomeJavaConstants.hasComponent);
			assertThat(referenceIdentifier(halfMappedComponents.get(0)), is("Q00001"));
			assertGhostOf(halfMappedComponents.get(1), prot99);

			assertThat(infer(halfMappedReaction), is(nullValue()));
			assertThat(reactionInferrer.getEligibleCount(), is(2));
			assertThat(reactionInferrer.getInferredCount(), is(1));
		}
	}

	@Nested
	@DisplayName("EWAS details")
	class EwasDetails {

		@Test
		void homologueWithoutEnsemblGeneIsLeftOut() throws Exception {
			// P00020 has two homologues, but only Q00020 has an Ensembl gene, so no DefinedSet is made.
			SimpleInstance reaction = fx.reaction("PROT20 => PROT1")
				.input(fx.ewas("PROT20", "P00020")).output(fx.ewas("PROT1", "P00001")).build();

			SimpleInstance inferredInput = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(inferredInput.getSchemaClassName(), is(ReactomeJavaConstants.EntityWithAccessionedSequence));
			assertThat(referenceIdentifier(inferredInput), is("Q00020"));
		}

		@Test
		void ewasWhoseOnlyHomologueHasNoEnsemblGeneStopsInference() throws Exception {
			SimpleInstance reaction = fx.reaction("PROT22 => PROT1")
				.input(fx.ewas("PROT22", "P00022")).output(fx.ewas("PROT1", "P00001")).build();

			assertNotInferred(reaction, true);
		}

		@Test
		void coordinatesAreCopiedAndAddTheSourceNameToTheNames() throws Exception {
			SimpleInstance reaction = fx.reaction("PROT1 fragment => PROT2")
				.input(fx.ewasWithCoordinates("PROT1", "P00001", 10, 200)).output(fx.ewas("PROT2", "P00002")).build();

			SimpleInstance inferredInput = single(infer(reaction), ReactomeJavaConstants.input);

			assertThat(inferredInput.getAttribute(ReactomeJavaConstants.startCoordinate), is(10));
			assertThat(inferredInput.getAttribute(ReactomeJavaConstants.endCoordinate), is(200));
			assertThat(inferredInput.getAttribute(ReactomeJavaConstants.name), is(Arrays.asList("Tst1", "PROT1", "Q00001")));
		}

		@Test
		void interChainCrosslinkIsInferredOntoTheHomologue() throws Exception {
			SimpleInstance ubiquitin = fx.simpleEntity("ubiquitin");
			SimpleInstance ubiquitinReference = fx.referenceGeneProduct("P0CG48");
			SimpleInstance crosslink = fx.interChainCrosslink(ubiquitin, ubiquitinReference);
			SimpleInstance reaction = fx.reaction("Ub-PROT1 => PROT2")
				.input(fx.ewas("Ub-PROT1", "P00001", crosslink)).output(fx.ewas("PROT2", "P00002")).build();

			SimpleInstance inferredInput = single(infer(reaction), ReactomeJavaConstants.input);
			SimpleInstance inferredCrosslink = single(inferredInput, ReactomeJavaConstants.hasModifiedResidue);

			assertThat(inferredCrosslink.getSchemaClassName(), is(ReactomeJavaConstants.InterChainCrosslinkedResidue));
			assertThat(single(inferredCrosslink, ReactomeJavaConstants.referenceSequence).getDbId(),
				is(single(inferredInput, ReactomeJavaConstants.referenceEntity).getDbId()));
			assertThat(single(inferredCrosslink, ReactomeJavaConstants.modification).getDbId(), is(ubiquitin.getDbId()));
			assertThat(single(inferredCrosslink, ReactomeJavaConstants.psiMod).getDbId(),
				is(single(crosslink, ReactomeJavaConstants.psiMod).getDbId()));
			assertThat(dbIds(values(inferredCrosslink, ReactomeJavaConstants.secondReferenceSequence)),
				contains(ubiquitinReference.getDbId()));
			// With no coordinate, it is named after what it links.
			assertThat(inferredCrosslink.getDisplayName(), is("UniProt:Q00001 ubiquitin [cytosol] N6-glycyl-L-lysine"));
		}

		@Test
		void alternateReferenceDatabaseAddsASecondGeneReference() throws Exception {
			wire("{\"dbname\": [\"TestBase\"], \"url\": \"https://testbase.example\", " +
				"\"access\": \"https://testbase.example/###ID###\", \"alt_id\": \"ENSTSTG\"}");
			SimpleInstance reaction = fx.reaction("PROT1 => PROT2")
				.input(fx.ewas("PROT1", "P00001")).output(fx.ewas("PROT2", "P00002")).build();

			SimpleInstance inferredInput = single(infer(reaction), ReactomeJavaConstants.input);
			List<SimpleInstance> genes =
				inflated(single(inferredInput, ReactomeJavaConstants.referenceEntity), ReactomeJavaConstants.referenceGene);

			assertThat(genes.stream().map(gene -> gene.getAttribute(ReactomeJavaConstants.identifier))
				.collect(Collectors.toList()), contains("ENSTSTG00001", "00001"));
			SimpleInstance alternateGene = genes.get(1);
			assertThat(alternateGene.getDisplayName(), is("TestBase:ENSTSTG00001"));
			assertThat(single(alternateGene, ReactomeJavaConstants.referenceDatabase).getDisplayName(), is("TestBase"));
			assertThat(dbIds(values(alternateGene, ReactomeJavaConstants.species)), contains(targetSpecies.getDbId()));
		}
	}

	// ---- Helpers ------------------------------------------------------------------------------------------------

	private SimpleInstance inputOfClass(SimpleInstance inferredReaction, String schemaClass) {
		return inflated(inferredReaction, ReactomeJavaConstants.input).stream()
			.filter(input -> input.getSchemaClassName().equals(schemaClass))
			.findFirst()
			.orElseThrow(() -> new AssertionError("No " + schemaClass + " input in " + inferredReaction));
	}

	private void assertGhostOf(SimpleInstance ghost, SimpleInstance humanEntity) {
		assertThat(ghost.getSchemaClassName(), is(ReactomeJavaConstants.GenomeEncodedEntity));
		assertThat(ghost.getDisplayName(), is("Ghost homologue of " + humanEntity.getDisplayName()));
		assertThat(dbIds(values(ghost, ReactomeJavaConstants.species)), contains(targetSpecies.getDbId()));
		assertThat(dbIds(values(ghost, ReactomeJavaConstants.inferredFrom)), hasItem(humanEntity.getDbId()));
	}

	// The number of instances of the class this run created whose reference sequence (or own identifier) is the one
	// given.
	private long newInstancesFor(String schemaClass, String identifier) {
		return api.committedInstancesOfClass(schemaClass).stream()
			.filter(instance -> identifier.equals(instance.getAttribute(ReactomeJavaConstants.identifier)) ||
				(instance.getAttribute(ReactomeJavaConstants.referenceEntity) != null &&
					identifier.equals(referenceIdentifier(instance))))
			.count();
	}


	// Infers the reaction and returns the inferred reaction as stored, or null if none was inferred.
	private SimpleInstance infer(SimpleInstance reaction) throws Exception {
		reactionInferrer.inferReaction(reaction);
		SimpleInstance inferred = reactionInferrer.getInferredEvent().get(reaction.getDbId());
		return inferred != null ? api.stored(inferred.getDbId()) : null;
	}

	private void assertNotInferred(SimpleInstance reaction, boolean eligible) throws Exception {
		assertThat(infer(reaction), is(nullValue()));
		assertThat(api.committedInstancesOfClass(ReactomeJavaConstants.Reaction), is(empty()));
		assertThat(values(reaction, ReactomeJavaConstants.orthologousEvent), is(empty()));
		assertThat(reactionInferrer.getEligibleCount(), is(eligible ? 1 : 0));
		assertThat(reactionInferrer.getInferredCount(), is(0));
	}

	private void assertInferredEwas(SimpleInstance inferredEwas, SimpleInstance humanEwas, String homologueId,
									String geneName) {
		assertThat(inferredEwas.getSchemaClassName(), is(ReactomeJavaConstants.EntityWithAccessionedSequence));
		assertThat(dbIds(values(inferredEwas, ReactomeJavaConstants.species)), contains(targetSpecies.getDbId()));
		assertThat(dbIds(values(inferredEwas, ReactomeJavaConstants.compartment)), contains(fx.cytosol.getDbId()));
		assertThat(dbIds(values(inferredEwas, ReactomeJavaConstants.inferredFrom)), hasItem(humanEwas.getDbId()));
		assertThat(stableIdentifier(inferredEwas), is("R-" + TARGET_ABBREVIATION + "-" + humanEwas.getDbId()));

		List<String> expectedNames = geneName != null ? Arrays.asList(geneName, homologueId) :
			Collections.singletonList(homologueId);
		assertThat(inferredEwas.getAttribute(ReactomeJavaConstants.name), is(expectedNames));
		assertThat(inferredEwas.getDisplayName(), is(expectedNames.get(0) + " [cytosol]"));

		SimpleInstance referenceGeneProduct = single(inferredEwas, ReactomeJavaConstants.referenceEntity);
		assertThat(referenceGeneProduct.getSchemaClassName(), is(ReactomeJavaConstants.ReferenceGeneProduct));
		assertThat(referenceGeneProduct.getAttribute(ReactomeJavaConstants.identifier), is(homologueId));
		assertThat(referenceGeneProduct.getDisplayName(), is("UniProt:" + homologueId));
		assertThat(single(referenceGeneProduct, ReactomeJavaConstants.referenceDatabase).getDbId(), is(fx.uniProt.getDbId()));
		assertThat(dbIds(values(referenceGeneProduct, ReactomeJavaConstants.species)), contains(targetSpecies.getDbId()));
		assertThat(inflated(referenceGeneProduct, ReactomeJavaConstants.referenceGene).stream()
				.map(gene -> gene.getAttribute(ReactomeJavaConstants.identifier)).collect(Collectors.toList()),
			contains("ENSTSTG" + homologueId.substring(1)));
	}

	private String referenceIdentifier(SimpleInstance ewas) {
		return harness.referenceIdentifier(ewas);
	}

	private String stableIdentifier(SimpleInstance instance) {
		return harness.stableIdentifier(instance);
	}

	private SimpleInstance single(SimpleInstance instance, String attribute) {
		return harness.single(instance, attribute);
	}

	private List<SimpleInstance> inflated(SimpleInstance instance, String attribute) {
		return harness.inflated(instance, attribute);
	}

	private List<SimpleInstance> values(SimpleInstance instance, String attribute) {
		return harness.values(instance, attribute);
	}
}
