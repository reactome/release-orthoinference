package org.reactome.orthoinference;

import org.gk.model.GKInstance;
import org.gk.persistence.MySQLAdaptor;
import org.junit.Before;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;
import org.junit.runner.RunWith;
import org.mockito.Mock;
import org.mockito.Mockito;
import org.mockito.MockitoAnnotations;
import org.powermock.api.mockito.PowerMockito;
import org.powermock.core.classloader.annotations.PowerMockIgnore;
import org.powermock.core.classloader.annotations.PrepareForTest;
import org.powermock.modules.junit4.PowerMockRunner;
import static org.junit.Assert.*;

import java.util.Collection;

@RunWith(PowerMockRunner.class)
@PrepareForTest({StableIdentifierGenerator.class, InstanceUtilities.class})
@PowerMockIgnore({"org.apache.logging.log4j.*", "javax.management.*", "javax.script.*",
        "javax.xml.*", "com.sun.org.apache.xerces.*", "org.xml.sax.*", "com.sun.xml.*", "org.w3c.dom.*", "org.mockito.*"})
public class StableIdentifierGeneratorTest {

    private static Object identifier = "R-HSA-123456";

    @Mock
    MySQLAdaptor mockAdaptor;

    @Mock
    GKInstance mockInferredInst;

    @Mock
    Collection<GKInstance> mockInstanceCollection;
    @Mock
    GKInstance mockOriginalInst;

    @Mock
    GKInstance mockStableIdentifierInst;

    @Mock
    GKInstance mockOrthoStableIdentifierInst;

    StableIdentifierGenerator stIdGenerator;

    @Before
    public void setUp() {
        MockitoAnnotations.initMocks(this);
        stIdGenerator = new StableIdentifierGenerator(mockAdaptor, "ABC");
    }

    @Test
    public void generateOrthologousStableIdReturnsStableIdentifierInst() throws Exception {

        PowerMockito.mockStatic(InstanceUtilities.class);
        Mockito.when(mockOriginalInst.getAttributeValue("stableIdentifier")).thenReturn(mockStableIdentifierInst);
        Mockito.when(mockStableIdentifierInst.getAttributeValue("identifier")).thenReturn(identifier);
        Mockito.when(mockAdaptor.fetchInstanceByAttribute("StableIdentifier", "identifier", "=", "R-ABC-123456")).thenReturn(mockInstanceCollection);
        PowerMockito.when(InstanceUtilities.createNewInferredGKInstance(mockStableIdentifierInst)).thenReturn(mockOrthoStableIdentifierInst);
        assertEquals(stIdGenerator.generateOrthologousStableId(mockInferredInst, mockOriginalInst), mockOrthoStableIdentifierInst);
    }

    @Test
    public void nullStableIdentifierReturnsRuntimeException() throws Exception {
        Mockito.when(mockOriginalInst.getAttributeValue("stableIdentifier")).thenReturn(null);
        try {
            stIdGenerator.generateOrthologousStableId(mockInferredInst, mockOriginalInst);
        } catch (RuntimeException e) {
            e.printStackTrace();
            assertTrue(e.getMessage().contains("No stable identifier instance found for " + mockOriginalInst));
        }
    }
}