package org.reactome.orthoinference;

import org.gk.model.ReactomeJavaConstants;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Nested;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.reactome.curation.model.SimpleInstance;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

/**
 * Runs representative human reactions through the real inference classes -- ReactionInferrer,
 * OrthologousEntityGenerator, EWASInferrer, ProteinCountUtility, SkipInstanceChecker and the rest -- against an
 * in-memory curator-tool-ws ({@link InMemoryCuratorToolAPI}) and a small orthopairs fixture, and checks what is
 * stored for the target species.
 *
 * <p>Homology fixture (human UniProt -> target UniProt):
 * <pre>
 *   P00001..P00005 -> Q00001..Q00005   (one homologue each; Q00001 and Q00002 have gene names)
 *   P00010         -> Q00010, Q00011   (two homologues)
 *   P00098, P00099 -> none
 * </pre>
 */
class ReactionInferenceScenarioTest {

	private static final String TARGET_SPECIES_CODE = "ttst";
	private static final String TARGET_SPECIES_NAME = "Testus testus";
	private static final String TARGET_ABBREVIATION = "TST";
	private static final String RELEASE_DATE = "2026-12-01";
	private static final long DISEASE_PATHWAY_DB_ID = 1643685L;

	@TempDir
	Path tempDir;

	private InMemoryCuratorToolAPI api;
	private HumanFixtures fx;
	private InstanceUtilities utils;
	private ReactionInferrer reactionInferrer;
	private SimpleInstance targetSpecies;

	@BeforeEach
	void setUp() throws Exception {
		// The inference classes cache what they infer in static fields, which would carry over between tests.
		clearStaticState(ReactionInferrer.class, OrthologousEntityGenerator.class, EWASInferrer.class,
			StableIdentifierGenerator.class);

		ConfigProperties config = new ConfigProperties(writeConfigFixtures());
		api = new InMemoryCuratorToolAPI(config);
		fx = new HumanFixtures(api);

		// Wired as Spring wires them.
		StableIdentifierGenerator stableIdGenerator = new StableIdentifierGenerator(config, TARGET_SPECIES_CODE, api);
		utils = new InstanceUtilities(config, TARGET_SPECIES_CODE, stableIdGenerator, api);
		// A @Bean in production, so created at start-up; createNewInferredGKInstance relies on it being set.
		targetSpecies = utils.getSpeciesInstance();
		Mappings mappings = new Mappings(config, TARGET_SPECIES_CODE);
		OrthologyReferenceDatabase referenceDatabase =
			new OrthologyReferenceDatabase(config, TARGET_SPECIES_CODE, utils, api);
		EWASInferrer ewasInferrer = new EWASInferrer(config, referenceDatabase, mappings, utils);
		ProteinCountUtility proteinCountUtility = new ProteinCountUtility(mappings, api, utils);
		OrthologousEntityGenerator entityGenerator = new OrthologousEntityGenerator(
			ewasInferrer, new SpeciesCheckUtility(utils), utils, proteinCountUtility);
		reactionInferrer = new ReactionInferrer(config, TARGET_SPECIES_CODE, entityGenerator, utils,
			proteinCountUtility, stableIdGenerator, new SkipInstanceChecker(utils, api), api);
	}

	@AfterEach
	void deleteEligibleAndInferredFiles() throws Exception {
		Files.deleteIfExists(eligibleFile());
		Files.deleteIfExists(inferredFile());
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

			// Linked to its source both ways.
			assertThat(dbIds(values(inferred, ReactomeJavaConstants.inferredFrom)), contains(reaction.getDbId()));
			assertThat(dbIds(values(inferred, ReactomeJavaConstants.orthologousEvent)), contains(reaction.getDbId()));
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

	// ---- Helpers ------------------------------------------------------------------------------------------------

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
		return (String) single(ewas, ReactomeJavaConstants.referenceEntity).getAttribute(ReactomeJavaConstants.identifier);
	}

	private String stableIdentifier(SimpleInstance instance) {
		return (String) single(instance, ReactomeJavaConstants.stableIdentifier).getAttribute(ReactomeJavaConstants.identifier);
	}

	// The attribute's only value, as stored.
	private SimpleInstance single(SimpleInstance instance, String attribute) {
		List<SimpleInstance> attributeValues = inflated(instance, attribute);
		assertThat(attribute + " of " + instance, attributeValues, hasSize(1));
		return attributeValues.get(0);
	}

	private List<SimpleInstance> inflated(SimpleInstance instance, String attribute) {
		return values(instance, attribute).stream().map(api::inflate).collect(Collectors.toList());
	}

	private List<SimpleInstance> values(SimpleInstance instance, String attribute) {
		return utils.getAttributeValues(instance, attribute);
	}

	private static List<Long> dbIds(List<SimpleInstance> instances) {
		return instances.stream().map(SimpleInstance::getDbId).collect(Collectors.toList());
	}

	// ---- Fixtures -----------------------------------------------------------------------------------------------

	private ConfigProperties.OrthoinferenceProps writeConfigFixtures() throws Exception {
		Path speciesConfig = tempDir.resolve("Species.json");
		Files.write(speciesConfig, Collections.singletonList(
			"{\"" + TARGET_SPECIES_CODE + "\": {\"name\": [\"" + TARGET_SPECIES_NAME + "\"], " +
				"\"abbreviation\": \"" + TARGET_ABBREVIATION + "\"}}"));

		Path orthopairs = Files.createDirectory(tempDir.resolve("orthopairs"));
		Files.write(orthopairs.resolve("hsap_" + TARGET_SPECIES_CODE + "_mapping.tsv"), Arrays.asList(
			"P00001\tUniProt:Q00001",
			"P00002\tUniProt:Q00002",
			"P00003\tUniProt:Q00003",
			"P00004\tUniProt:Q00004",
			"P00005\tUniProt:Q00005",
			"P00010\tUniProt:Q00010 UniProt:Q00011"));
		Files.write(orthopairs.resolve(TARGET_SPECIES_CODE + "_gene_protein_mapping.tsv"),
			Arrays.asList("00001", "00002", "00003", "00004", "00005", "00010", "00011").stream()
				.map(id -> "ENSTSTG" + id + "\tUniProt:Q" + id)
				.collect(Collectors.toList()));
		Files.write(orthopairs.resolve(TARGET_SPECIES_CODE + "_gene_name_mapping.tsv"),
			Arrays.asList("Q00001\tTst1", "Q00002\tTst2"));

		ConfigProperties.OrthoinferenceProps props = new ConfigProperties.OrthoinferenceProps();
		props.setPathToSpeciesConfig(speciesConfig.toString());
		props.setPathToOrthopairs(orthopairs.toString());
		props.setSourceSpeciesCode("hsap");
		props.setDateOfRelease(RELEASE_DATE);
		props.setPersonId(8_888_888L);
		return props;
	}

	private static Path eligibleFile() {
		return Paths.get("eligible_" + TARGET_SPECIES_CODE + "_75.txt");
	}

	private static Path inferredFile() {
		return Paths.get("inferred_" + TARGET_SPECIES_CODE + "_75.txt");
	}

	// Empties the static collections and zeroes the static counters of each class.
	private static void clearStaticState(Class<?>... classes) throws IllegalAccessException {
		for (Class<?> clazz : classes) {
			for (Field field : clazz.getDeclaredFields()) {
				if (!Modifier.isStatic(field.getModifiers())) {
					continue;
				}
				field.setAccessible(true);
				Object value = field.get(null);
				if (value instanceof Map) {
					((Map<?, ?>) value).clear();
				} else if (value instanceof Collection) {
					((Collection<?>) value).clear();
				} else if (field.getType() == Integer.class && !Modifier.isFinal(field.getModifiers())) {
					field.set(null, 0);
				}
			}
		}
	}
}
