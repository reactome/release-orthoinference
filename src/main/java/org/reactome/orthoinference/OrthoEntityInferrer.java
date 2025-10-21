package org.reactome.orthoinference;

import org.gk.model.GKInstance;

public interface OrthoEntityInferrer {

    GKInstance createOrthoEntity(GKInstance entityInst, boolean override) throws Exception;
}