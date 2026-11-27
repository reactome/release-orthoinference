package org.reactome.orthoinference;

import org.gk.model.ReactomeJavaConstants;
import org.reactome.curation.model.SimpleInstance;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Builds human (source) instances and stores them in an {@link InMemoryCuratorToolAPI}, shaped as curator-tool-ws
 * returns them: multi-valued attributes as Lists, and every PhysicalEntity and Event with a R-HSA StableIdentifier.
 */
class HumanFixtures {

	private final InMemoryCuratorToolAPI api;
	private long nextDbId = 100_000L;

	final SimpleInstance uniProt;
	final SimpleInstance humanSpecies;
	final SimpleInstance cytosol;
	final SimpleInstance catalyticActivity;

	HumanFixtures(InMemoryCuratorToolAPI api) {
		this.api = api;

		uniProt = instance(ReactomeJavaConstants.ReferenceDatabase, "UniProt");
		uniProt.setAttribute(ReactomeJavaConstants.name, list("UniProt"));
		store(uniProt);

		humanSpecies = instance(ReactomeJavaConstants.Species, "Homo sapiens");
		humanSpecies.setAttribute(ReactomeJavaConstants.name, list("Homo sapiens"));
		store(humanSpecies);

		cytosol = instance(ReactomeJavaConstants.Compartment, "cytosol");
		cytosol.setAttribute(ReactomeJavaConstants.name, list("cytosol"));
		cytosol.setAttribute("accession", "0005829");
		store(cytosol);

		catalyticActivity = instance("GO_MolecularFunction", "catalytic activity");
		catalyticActivity.setAttribute(ReactomeJavaConstants.name, list("catalytic activity"));
		catalyticActivity.setAttribute("accession", "0003824");
		store(catalyticActivity);
	}

	// ---- Physical entities --------------------------------------------------------------------------------------

	SimpleInstance referenceGeneProduct(String uniProtId) {
		SimpleInstance rgp = instance(ReactomeJavaConstants.ReferenceGeneProduct, "UniProt:" + uniProtId);
		rgp.setAttribute(ReactomeJavaConstants.identifier, uniProtId);
		rgp.setAttribute(ReactomeJavaConstants.referenceDatabase, uniProt);
		rgp.setAttribute(ReactomeJavaConstants.species, list(humanSpecies));
		return store(rgp);
	}

	SimpleInstance ewas(String name, String uniProtId, SimpleInstance... modifiedResidues) {
		SimpleInstance ewas = physicalEntity(ReactomeJavaConstants.EntityWithAccessionedSequence, name);
		ewas.setAttribute(ReactomeJavaConstants.species, list(humanSpecies));
		ewas.setAttribute(ReactomeJavaConstants.referenceEntity, referenceGeneProduct(uniProtId));
		if (modifiedResidues.length > 0) {
			ewas.setAttribute(ReactomeJavaConstants.hasModifiedResidue, list(modifiedResidues));
		}
		return storeWithStableId(ewas);
	}

	SimpleInstance phosphoSerine(int coordinate) {
		SimpleInstance psiMod = instance("PsiMod", "O-phospho-L-serine");
		psiMod.setAttribute(ReactomeJavaConstants.name, list("O-phospho-L-serine"));
		psiMod.setAttribute(ReactomeJavaConstants.identifier, "00046");
		store(psiMod);

		SimpleInstance residue = instance(ReactomeJavaConstants.ModifiedResidue, "O-phospho-L-serine at " + coordinate);
		residue.setAttribute(ReactomeJavaConstants.coordinate, coordinate);
		residue.setAttribute(ReactomeJavaConstants.psiMod, psiMod);
		return store(residue);
	}

	// A small molecule: no species, so it is used as it is rather than inferred.
	SimpleInstance simpleEntity(String name) {
		return storeWithStableId(physicalEntity(ReactomeJavaConstants.SimpleEntity, name));
	}

	SimpleInstance complex(String name, SimpleInstance... components) {
		SimpleInstance complex = physicalEntity(ReactomeJavaConstants.Complex, name);
		complex.setAttribute(ReactomeJavaConstants.species, list(humanSpecies));
		complex.setAttribute(ReactomeJavaConstants.hasComponent, list(components));
		return storeWithStableId(complex);
	}

	SimpleInstance definedSet(String name, SimpleInstance... members) {
		SimpleInstance set = physicalEntity(ReactomeJavaConstants.DefinedSet, name);
		set.setAttribute(ReactomeJavaConstants.species, list(humanSpecies));
		set.setAttribute(ReactomeJavaConstants.hasMember, list(members));
		return storeWithStableId(set);
	}

	SimpleInstance candidateSet(String name, List<SimpleInstance> members, List<SimpleInstance> candidates) {
		SimpleInstance set = physicalEntity(ReactomeJavaConstants.CandidateSet, name);
		set.setAttribute(ReactomeJavaConstants.species, list(humanSpecies));
		if (!members.isEmpty()) {
			set.setAttribute(ReactomeJavaConstants.hasMember, new ArrayList<>(members));
		}
		set.setAttribute(ReactomeJavaConstants.hasCandidate, new ArrayList<>(candidates));
		return storeWithStableId(set);
	}

	// ---- Catalysts and regulations ------------------------------------------------------------------------------

	SimpleInstance catalystActivity(SimpleInstance physicalEntity, SimpleInstance... activeUnits) {
		SimpleInstance catalyst = instance(ReactomeJavaConstants.CatalystActivity,
			"catalytic activity of " + physicalEntity.getDisplayName());
		catalyst.setAttribute(ReactomeJavaConstants.activity, catalyticActivity);
		catalyst.setAttribute(ReactomeJavaConstants.physicalEntity, physicalEntity);
		if (activeUnits.length > 0) {
			catalyst.setAttribute(ReactomeJavaConstants.activeUnit, list(activeUnits));
		}
		return store(catalyst);
	}

	SimpleInstance regulation(String regulationClass, SimpleInstance regulator) {
		SimpleInstance regulation = instance(regulationClass,
			regulationClass + " by " + regulator.getDisplayName());
		regulation.setAttribute(ReactomeJavaConstants.regulator, regulator);
		return store(regulation);
	}

	// ---- Events -------------------------------------------------------------------------------------------------

	ReactionBuilder reaction(String name) {
		return new ReactionBuilder(name);
	}

	SimpleInstance pathway(long dbId, String name, SimpleInstance... events) {
		SimpleInstance pathway = instance(ReactomeJavaConstants.Pathway, name);
		pathway.setDbId(dbId);
		pathway.setAttribute(ReactomeJavaConstants.species, list(humanSpecies));
		pathway.setAttribute(ReactomeJavaConstants.hasEvent, list(events));
		return api.add(pathway);
	}

	SimpleInstance disease(String name) {
		SimpleInstance disease = instance(ReactomeJavaConstants.Disease, name);
		disease.setAttribute(ReactomeJavaConstants.name, list(name));
		return store(disease);
	}

	class ReactionBuilder {
		private final SimpleInstance reaction;

		private ReactionBuilder(String name) {
			reaction = instance(ReactomeJavaConstants.Reaction, name);
			reaction.setAttribute(ReactomeJavaConstants.name, list(name));
			reaction.setAttribute(ReactomeJavaConstants.species, list(humanSpecies));
			reaction.setAttribute(ReactomeJavaConstants.compartment, list(cytosol));
		}

		ReactionBuilder input(SimpleInstance... inputs) {
			return add(ReactomeJavaConstants.input, inputs);
		}

		ReactionBuilder output(SimpleInstance... outputs) {
			return add(ReactomeJavaConstants.output, outputs);
		}

		ReactionBuilder catalyst(SimpleInstance... catalystActivities) {
			return add(ReactomeJavaConstants.catalystActivity, catalystActivities);
		}

		ReactionBuilder regulatedBy(SimpleInstance... regulations) {
			return add(ReactomeJavaConstants.regulatedBy, regulations);
		}

		ReactionBuilder with(String attribute, Object value) {
			reaction.setAttribute(attribute, value);
			return this;
		}

		/**
		 * Stores the reaction and returns it as the inference receives it: inflated, referring to everything else
		 * through shells.
		 */
		SimpleInstance build() {
			storeWithStableId(reaction);
			return api.inflate(reaction);
		}

		@SuppressWarnings("unchecked")
		private ReactionBuilder add(String attribute, SimpleInstance... values) {
			List<SimpleInstance> current = (List<SimpleInstance>) reaction.getAttribute(attribute);
			List<SimpleInstance> updated = current != null ? new ArrayList<>(current) : new ArrayList<>();
			updated.addAll(Arrays.asList(values));
			reaction.setAttribute(attribute, updated);
			return this;
		}
	}

	// ---- Helpers ------------------------------------------------------------------------------------------------

	private SimpleInstance physicalEntity(String schemaClass, String name) {
		SimpleInstance entity = instance(schemaClass, name + " [" + cytosol.getDisplayName() + "]");
		entity.setAttribute(ReactomeJavaConstants.name, list(name));
		entity.setAttribute(ReactomeJavaConstants.compartment, list(cytosol));
		return entity;
	}

	private SimpleInstance instance(String schemaClass, String displayName) {
		SimpleInstance instance = new SimpleInstance();
		instance.setDbId(nextDbId++);
		instance.setSchemaClassName(schemaClass);
		instance.setDisplayName(displayName);
		return instance;
	}

	private SimpleInstance storeWithStableId(SimpleInstance instance) {
		String identifier = "R-HSA-" + instance.getDbId();
		SimpleInstance stableIdentifier = instance(ReactomeJavaConstants.StableIdentifier, identifier + ".1");
		stableIdentifier.setAttribute(ReactomeJavaConstants.identifier, identifier);
		stableIdentifier.setAttribute(ReactomeJavaConstants.identifierVersion, "1");
		store(stableIdentifier);
		instance.setAttribute(ReactomeJavaConstants.stableIdentifier, stableIdentifier);
		return store(instance);
	}

	private SimpleInstance store(SimpleInstance instance) {
		return api.add(instance);
	}

	@SafeVarargs
	private static <T> List<T> list(T... values) {
		return new ArrayList<>(Arrays.asList(values));
	}

	static List<SimpleInstance> none() {
		return Collections.emptyList();
	}
}
