package org.reactome.orthoinference;

import org.gk.model.ReactomeJavaConstants;
import org.reactome.curation.model.SimpleInstance;

import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.hasSize;

/**
 * The inference classes wired together as Spring wires them, against an {@link InMemoryCuratorToolAPI} and a small
 * species configuration and orthopairs fixture written to a temporary directory, plus helpers for reading what was
 * stored.
 *
 * <p>Homology fixture (human UniProt -> target UniProt):
 * <pre>
 *   P00001..P00005 -> Q00001..Q00005   (one homologue each; Q00001 and Q00002 have gene names)
 *   P00010         -> Q00010, Q00011   (two homologues)
 *   P00020         -> Q00020, Q00021   (Q00021 has no Ensembl gene in the gene-protein mapping)
 *   P00022         -> Q00022           (which has no Ensembl gene)
 *   P00098, P00099 -> none
 * </pre>
 * Each target protein Qnnnnn is encoded by the Ensembl gene ENSTSTGnnnnn.
 *
 * <p>A made-up target species code is used because the inference writes eligible_/inferred_ files and a report
 * named after it, and its release number, into the working directory; {@link #deleteOutputFiles()} removes them.
 */
class InferenceTestHarness {

	static final String TARGET_SPECIES_CODE = "ttst";
	static final String TARGET_SPECIES_NAME = "Testus testus";
	static final String TARGET_ABBREVIATION = "TST";
	static final String RELEASE_DATE = "2026-12-01";
	static final int RELEASE_NUMBER = 9999;
	static final long DISEASE_PATHWAY_DB_ID = 1643685L;

	final InMemoryCuratorToolAPI api;
	final HumanFixtures fx;
	final InstanceUtilities utils;
	final ReactionInferrer reactionInferrer;
	final PathwaysInferrer pathwaysInferrer;
	final EventsInferrer eventsInferrer;
	final SimpleInstance targetSpecies;

	/**
	 * Builds the inference classes afresh against a new, empty in-memory database.
	 *
	 * @param tempDir - where the species configuration and orthopairs files are written.
	 * @param alternateReferenceDatabaseJson - the species' "alt_refdb" entry in Species.json, or null for none.
	 */
	InferenceTestHarness(Path tempDir, String alternateReferenceDatabaseJson) throws Exception {
		// The inference classes cache what they infer in static fields, which would carry over between tests.
		clearStaticState(ReactionInferrer.class, OrthologousEntityGenerator.class, EWASInferrer.class,
			StableIdentifierGenerator.class, PathwaysInferrer.class, EventsInferrer.class);

		ConfigProperties config = new ConfigProperties(writeConfigFixtures(tempDir, alternateReferenceDatabaseJson));
		api = new InMemoryCuratorToolAPI(config);
		fx = new HumanFixtures(api);

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
		pathwaysInferrer = new PathwaysInferrer(config, stableIdGenerator, utils, api);
		eventsInferrer = new EventsInferrer(config, TARGET_SPECIES_CODE, reactionInferrer, pathwaysInferrer, utils, api);
	}

	static void deleteOutputFiles() throws Exception {
		Files.deleteIfExists(eligibleFile());
		Files.deleteIfExists(inferredFile());
		Files.deleteIfExists(reportFile());
	}

	static Path eligibleFile() {
		return Paths.get("eligible_" + TARGET_SPECIES_CODE + "_75.txt");
	}

	static Path inferredFile() {
		return Paths.get("inferred_" + TARGET_SPECIES_CODE + "_75.txt");
	}

	static Path reportFile() {
		return Paths.get("report_ortho_inference_test_reactome_" + RELEASE_NUMBER + ".txt");
	}

	// ---- Reading what was stored --------------------------------------------------------------------------------

	// The attribute's only value, as stored.
	SimpleInstance single(SimpleInstance instance, String attribute) {
		List<SimpleInstance> attributeValues = inflated(instance, attribute);
		assertThat(attribute + " of " + instance, attributeValues, hasSize(1));
		return attributeValues.get(0);
	}

	List<SimpleInstance> inflated(SimpleInstance instance, String attribute) {
		return values(instance, attribute).stream().map(api::inflate).collect(Collectors.toList());
	}

	List<SimpleInstance> values(SimpleInstance instance, String attribute) {
		return utils.getAttributeValues(instance, attribute);
	}

	static List<Long> dbIds(List<SimpleInstance> instances) {
		return instances.stream().map(SimpleInstance::getDbId).collect(Collectors.toList());
	}

	String referenceIdentifier(SimpleInstance ewas) {
		return (String) single(ewas, ReactomeJavaConstants.referenceEntity).getAttribute(ReactomeJavaConstants.identifier);
	}

	String stableIdentifier(SimpleInstance instance) {
		return (String) single(instance, ReactomeJavaConstants.stableIdentifier).getAttribute(ReactomeJavaConstants.identifier);
	}

	// ---- Fixtures -----------------------------------------------------------------------------------------------

	private static ConfigProperties.OrthoinferenceProps writeConfigFixtures(
		Path tempDir, String alternateReferenceDatabaseJson) throws Exception {

		Path speciesConfig = tempDir.resolve("Species.json");
		Files.write(speciesConfig, Collections.singletonList(
			"{\"" + TARGET_SPECIES_CODE + "\": {\"name\": [\"" + TARGET_SPECIES_NAME + "\"], " +
				"\"abbreviation\": \"" + TARGET_ABBREVIATION + "\"" +
				(alternateReferenceDatabaseJson != null ? ", \"alt_refdb\": " + alternateReferenceDatabaseJson : "") +
				"}}"));

		Path orthopairs = Files.createDirectories(tempDir.resolve("orthopairs"));
		Files.write(orthopairs.resolve("hsap_" + TARGET_SPECIES_CODE + "_mapping.tsv"), Arrays.asList(
			"P00001\tUniProt:Q00001",
			"P00002\tUniProt:Q00002",
			"P00003\tUniProt:Q00003",
			"P00004\tUniProt:Q00004",
			"P00005\tUniProt:Q00005",
			"P00010\tUniProt:Q00010 UniProt:Q00011",
			"P00020\tUniProt:Q00020 UniProt:Q00021",
			"P00022\tUniProt:Q00022"));
		// Q00021 and Q00022 are left out on purpose: a homologue without an Ensembl gene is not inferred.
		Files.write(orthopairs.resolve(TARGET_SPECIES_CODE + "_gene_protein_mapping.tsv"),
			Arrays.asList("00001", "00002", "00003", "00004", "00005", "00010", "00011", "00020").stream()
				.map(id -> "ENSTSTG" + id + "\tUniProt:Q" + id)
				.collect(Collectors.toList()));
		Files.write(orthopairs.resolve(TARGET_SPECIES_CODE + "_gene_name_mapping.tsv"),
			Arrays.asList("Q00001\tTst1", "Q00002\tTst2"));

		ConfigProperties.OrthoinferenceProps props = new ConfigProperties.OrthoinferenceProps();
		props.setPathToSpeciesConfig(speciesConfig.toString());
		props.setPathToOrthopairs(orthopairs.toString());
		props.setSourceSpeciesCode("hsap");
		props.setDateOfRelease(RELEASE_DATE);
		props.setReleaseNumber(RELEASE_NUMBER);
		props.setPersonId(8_888_888L);
		return props;
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
