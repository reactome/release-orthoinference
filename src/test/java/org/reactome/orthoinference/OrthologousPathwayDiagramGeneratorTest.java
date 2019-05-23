package org.reactome.orthoinference;


import org.gk.model.GKInstance;
import org.gk.model.ReactomeJavaConstants;
import org.gk.pathwaylayout.PredictedPathwayDiagramGeneratorFromDB;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.mockito.Mock;
import org.mockito.Mockito;
import org.mockito.MockitoAnnotations;
import org.powermock.api.mockito.PowerMockito;
import org.powermock.core.classloader.annotations.PowerMockIgnore;
import org.powermock.core.classloader.annotations.PrepareForTest;
import org.powermock.modules.junit4.PowerMockRunner;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

@RunWith(PowerMockRunner.class)
@PrepareForTest({OrthologousPathwayDiagramGenerator.class, PredictedPathwayDiagramGeneratorFromDB.class})
@PowerMockIgnore({"org.apache.logging.log4j.*", "javax.management.*", "javax.script.*",
        "javax.xml.*", "com.sun.org.apache.xerces.*", "org.xml.sax.*", "com.sun.xml.*", "org.w3c.dom.*", "org.mockito.*"})
public class OrthologousPathwayDiagramGeneratorTest {

    @Mock
    MySQLAdaptor mockAdaptor;
    @Mock
    MySQLAdaptor mockPrevAdaptor;
    @Mock
    GKInstance mockSpeciesInst;
    @Mock
    GKInstance mockAlternateSpeciesInst;
    @Mock
    GKInstance mockDiagramInst;
    @Mock
    GKInstance mockPathwayInst;
    @Mock
    GKInstance mockOrthoEventInst;
    @Mock
    GKInstance mockEvidenceTypeInst;
    @Mock
    GKInstance mockOrthoDiagramInst;
    @Mock
    PredictedPathwayDiagramGeneratorFromDB mockDiagramGenerator;

    long mockId = 12345L;

    Collection<GKInstance> mockDiagramInstances = new ArrayList<>();
    Collection<GKInstance> mockPrevDiagramInstances = new ArrayList<>();
    List<GKInstance> mockOrthoEventInstances = new ArrayList<>();


    OrthologousPathwayDiagramGenerator testDiagramGenerator;


    @Before
    public void setUp() {
        MockitoAnnotations.initMocks(this);
        testDiagramGenerator = new OrthologousPathwayDiagramGenerator(mockAdaptor, mockPrevAdaptor, mockSpeciesInst, mockId, mockId);
    }
    @Test
    public void generateOrthologousPathwayDiagramsTest() throws Exception {

        mockDiagramInstances.add(mockDiagramInst);
        mockOrthoEventInstances.add(mockOrthoEventInst);

        PowerMockito.whenNew(PredictedPathwayDiagramGeneratorFromDB.class).withNoArguments().thenReturn(mockDiagramGenerator);
        Mockito.when(mockAdaptor.fetchInstance(mockId)).thenReturn(mockSpeciesInst);
        Mockito.when(mockSpeciesInst.getDisplayName()).thenReturn("Species");
        Mockito.when(mockAdaptor.fetchInstancesByClass(ReactomeJavaConstants.PathwayDiagram)).thenReturn(mockDiagramInstances);
        Mockito.when(mockDiagramInst.getAttributeValue(ReactomeJavaConstants.representedPathway)).thenReturn(mockPathwayInst);
        Mockito.when(mockPathwayInst.getAttributeValue(ReactomeJavaConstants.species)).thenReturn(mockSpeciesInst);
        Mockito.when(mockPathwayInst.getAttributeValuesList(ReactomeJavaConstants.orthologousEvent)).thenReturn(mockOrthoEventInstances);
        Mockito.when(mockOrthoEventInst.getAttributeValue(ReactomeJavaConstants.species)).thenReturn(mockSpeciesInst);
        Mockito.when(mockOrthoEventInst.getAttributeValue(ReactomeJavaConstants.evidenceType)).thenReturn(mockEvidenceTypeInst);
        Mockito.when(mockDiagramGenerator.generatePredictedDiagram(mockOrthoEventInst, mockPathwayInst, mockDiagramInst)).thenReturn(mockOrthoDiagramInst);
        testDiagramGenerator.generateOrthologousPathwayDiagrams();
    }

    @Test
    public void fewerPathwayDiagramsTest() throws Exception {
        mockDiagramInstances.add(mockDiagramInst);
        mockPrevDiagramInstances.addAll(Arrays.asList(mockDiagramInst, mockDiagramInst));

        Mockito.when(mockAdaptor.fetchInstancesByClass(ReactomeJavaConstants.PathwayDiagram)).thenReturn(mockDiagramInstances);
        Mockito.when(mockPrevAdaptor.fetchInstancesByClass(ReactomeJavaConstants.PathwayDiagram)).thenReturn(mockPrevDiagramInstances);
        Mockito.when(mockDiagramInst.getAttributeValue(ReactomeJavaConstants.representedPathway)).thenReturn(mockPathwayInst);
        Mockito.when(mockPathwayInst.getAttributeValue(ReactomeJavaConstants.species)).thenReturn(mockAlternateSpeciesInst);
        testDiagramGenerator.generateOrthologousPathwayDiagrams();

    }
}