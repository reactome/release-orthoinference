package org.reactome.orthoinference;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;

public class SpeciesCheckUtility {

	// Determines if there is a species attribute in any constituent instances of entityInst.
	// Unless it's an 'OtherEntity' (which will return false), the function will check the instance or iterate on its
	// sub-instances until it finds an existing 'species' attribute, or else it will return false.
	@SuppressWarnings("unchecked")
	public static boolean hasOrContainsSpeciesAttribute(GKInstance entity) throws Exception {
		if (isOtherEntity(entity)) {
			return false;
		}

		return isEntitySetContainingSpecies(entity) ||
			isComplexContainingSpecies(entity) ||
			isPolymerContainingSpecies(entity) ||
			hasSpecies(entity);
	}

	private static boolean isEntitySetContainingSpecies(GKInstance entity) throws Exception {
		if (isEntitySet(entity)) {
			boolean membersHaveSpecies =
				anyPhysicalEntityInAttributeHasSpecies(entity, ReactomeJavaConstants.hasMember);

			return isCandidateSet(entity) ?
				anyPhysicalEntityInAttributeHasSpecies(entity, ReactomeJavaConstants.hasCandidate) || membersHaveSpecies :
				membersHaveSpecies;
		}
		return false;
	}

	private static boolean isComplexContainingSpecies(GKInstance entity) throws Exception {
		return isComplex(entity) && anyPhysicalEntityInAttributeHasSpecies(entity, ReactomeJavaConstants.hasComponent);
	}

	private static boolean isPolymerContainingSpecies(GKInstance entity) throws Exception {
		return isPolymer(entity) && anyPhysicalEntityInAttributeHasSpecies(entity, ReactomeJavaConstants.repeatedUnit);
	}

	private static boolean isOtherEntity(GKInstance entity) {
		return isSchemaClassType(entity, ReactomeJavaConstants.OtherEntity);
	}

	private static boolean isEntitySet(GKInstance entity) {
		return isSchemaClassType(entity, ReactomeJavaConstants.EntitySet);
	}

	private static boolean isCandidateSet(GKInstance entity) {
		return isSchemaClassType(entity, ReactomeJavaConstants.CandidateSet);
	}

	private static boolean isComplex(GKInstance entity) {
		return isSchemaClassType(entity, ReactomeJavaConstants.Complex);
	}

	private static boolean isPolymer(GKInstance entity) {
		return isSchemaClassType(entity, ReactomeJavaConstants.Polymer);
	}

	private static boolean isSchemaClassType(GKInstance entity, String schemaClassType) {
		return entity.getSchemClass().isa(schemaClassType);
	}

	private static boolean hasSpecies(GKInstance entity) throws Exception {
		return entity.getSchemClass().isValidAttribute(ReactomeJavaConstants.species) &&
			entity.getAttributeValue(ReactomeJavaConstants.species) != null;
	}

	private static boolean anyPhysicalEntityInAttributeHasSpecies(GKInstance entity, String attribute)
		throws Exception {

		for (GKInstance physicalEntity : getPhysicalEntitiesByAttribute(entity, attribute)) {
			if (hasOrContainsSpeciesAttribute(physicalEntity)) {
				return true;
			}
		}
		return false;
	}

	private static List<GKInstance> getPhysicalEntitiesByAttribute(GKInstance entity, String attribute)
		throws Exception {

		Collection<GKInstance> attributeInstanceList = entity.getAttributeValuesList(attribute);

		if (attributeInstanceList == null) {
			return new ArrayList<>();
		}
		return new ArrayList<>(attributeInstanceList);
	}
}
