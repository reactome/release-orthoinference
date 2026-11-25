package org.reactome.orthoinference;

import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;

import org.gk.model.ReactomeJavaConstants;
import org.reactome.curation.model.SimpleInstance;
import org.springframework.stereotype.Component;

@Component
public class SpeciesCheckUtility {

	private InstanceUtilities instanceUtilities;

	public SpeciesCheckUtility(InstanceUtilities instanceUtilities) {
		this.instanceUtilities = instanceUtilities;
	}

	// Determines if there is a species attribute in any constituent instances of entityInst.
	// Unless it's an 'OtherEntity' (which will return false), the function will check the instance or iterate on its
	// sub-instances until it finds an existing 'species' attribute, or else it will return false.
	@SuppressWarnings("unchecked")
	// As in the GKInstance version, an EntitySet, Complex or Polymer is judged only by its constituents, not by its own
	// species: e.g. a Polymer of a SimpleEntity has species Homo sapiens but no constituent with a species, so it is
	// used as it is rather than inferred (which would fail the protein-count threshold).
	public boolean hasOrContainsSpeciesAttribute(SimpleInstance entity) throws Exception {
		if (isOtherEntity(entity)) {
			return false;
		}
		if (isEntitySet(entity)) {
			return anyPhysicalEntityInAttributeHasSpecies(entity, ReactomeJavaConstants.hasMember) ||
				(isCandidateSet(entity) &&
					anyPhysicalEntityInAttributeHasSpecies(entity, ReactomeJavaConstants.hasCandidate));
		}
		if (isComplex(entity)) {
			return anyPhysicalEntityInAttributeHasSpecies(entity, ReactomeJavaConstants.hasComponent);
		}
		if (isPolymer(entity)) {
			return anyPhysicalEntityInAttributeHasSpecies(entity, ReactomeJavaConstants.repeatedUnit);
		}
		return hasSpecies(entity);
	}

	private boolean isOtherEntity(SimpleInstance entity) {
		return instanceUtilities.isSchemaClassType(entity, ReactomeJavaConstants.OtherEntity);
	}

	private boolean isEntitySet(SimpleInstance entity) {
		return instanceUtilities.isSchemaClassType(entity, ReactomeJavaConstants.EntitySet);
	}

	private boolean isCandidateSet(SimpleInstance entity) {
		return instanceUtilities.isSchemaClassType(entity, ReactomeJavaConstants.CandidateSet);
	}

	private boolean isComplex(SimpleInstance entity) {
		return instanceUtilities.isSchemaClassType(entity, ReactomeJavaConstants.Complex);
	}

	private boolean isPolymer(SimpleInstance entity) {
		return instanceUtilities. isSchemaClassType(entity, ReactomeJavaConstants.Polymer);
	}

	private boolean hasSpecies(SimpleInstance entity) {
		return entity.getAttribute(ReactomeJavaConstants.species) != null;
	}

	private boolean anyPhysicalEntityInAttributeHasSpecies(SimpleInstance entity, String attribute)
		throws Exception {

		for (SimpleInstance physicalEntity : getPhysicalEntitiesByAttribute(entity, attribute)) {
			if (hasOrContainsSpeciesAttribute(physicalEntity)) {
				return true;
			}
		}
		return false;
	}

	// The attribute's values are shells (dbId, displayName and schemaClass only), which never carry a species or
	// constituents of their own, so each is inflated before it is checked.
	private List<SimpleInstance> getPhysicalEntitiesByAttribute(SimpleInstance entity, String attribute) {
		return instanceUtilities.getInflatedInstances(entity, attribute).stream()
			.filter(Objects::nonNull)
			.collect(Collectors.toList());
	}
}
