package org.reactome.orthoinference;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

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
	public boolean hasOrContainsSpeciesAttribute(SimpleInstance entity) throws Exception {
		if (isOtherEntity(entity)) {
			return false;
		}

		return isEntitySetContainingSpecies(entity) ||
			isComplexContainingSpecies(entity) ||
			isPolymerContainingSpecies(entity) ||
			hasSpecies(entity);
	}

	private boolean isEntitySetContainingSpecies(SimpleInstance entity) throws Exception {
		if (isEntitySet(entity)) {
			boolean membersHaveSpecies =
				anyPhysicalEntityInAttributeHasSpecies(entity, ReactomeJavaConstants.hasMember);

			return isCandidateSet(entity) ?
				anyPhysicalEntityInAttributeHasSpecies(entity, ReactomeJavaConstants.hasCandidate) || membersHaveSpecies :
				membersHaveSpecies;
		}
		return false;
	}

	private boolean isComplexContainingSpecies(SimpleInstance entity) throws Exception {
		return isComplex(entity) && anyPhysicalEntityInAttributeHasSpecies(entity, ReactomeJavaConstants.hasComponent);
	}

	private boolean isPolymerContainingSpecies(SimpleInstance entity) throws Exception {
		return isPolymer(entity) && anyPhysicalEntityInAttributeHasSpecies(entity, ReactomeJavaConstants.repeatedUnit);
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

	private List<SimpleInstance> getPhysicalEntitiesByAttribute(SimpleInstance entity, String attribute) {

		Collection<SimpleInstance> attributeInstanceList = (Collection<SimpleInstance>) entity.getAttribute(attribute);

		if (attributeInstanceList == null) {
			return new ArrayList<>();
		}
		return new ArrayList<>(attributeInstanceList);
	}
}
