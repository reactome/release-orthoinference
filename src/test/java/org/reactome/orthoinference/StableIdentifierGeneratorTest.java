package org.reactome.orthoinference;

import org.gk.model.ReactomeJavaConstants;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.extension.ExtendWith;
import org.mockito.Mock;
import org.mockito.Mockito;
import org.mockito.junit.jupiter.MockitoExtension;
import org.reactome.curation.model.SimpleInstance;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;
import static org.junit.jupiter.api.Assertions.assertThrows;

@ExtendWith(MockitoExtension.class)
public class StableIdentifierGeneratorTest {

    private static final String SPECIES_CODE = "abcd";
    private static final String SPECIES_ABBREVIATION = "ABC";

    @Mock
    ConfigProperties mockConfigProperties;

    @Mock
    SpeciesConfig mockSpeciesConfig;

    @Mock
    CuratorToolAPI mockCuratorToolAPI;

    StableIdentifierGenerator stIdGenerator;

    @BeforeEach
    public void setUp() throws Exception {
        // Lenient, as the tests that fail before an identifier is built never ask for the abbreviation.
        Mockito.lenient().when(mockConfigProperties.getSpeciesConfig()).thenReturn(mockSpeciesConfig);
        Mockito.lenient().when(mockSpeciesConfig.getAbbreviation(SPECIES_CODE)).thenReturn(SPECIES_ABBREVIATION);
        stIdGenerator = new StableIdentifierGenerator(mockConfigProperties, SPECIES_CODE, mockCuratorToolAPI);
    }

    // StableIdentifierGenerator counts paralogs in a static map that is kept across tests, so each test uses its
    // own source identifier to get a predictable target identifier.

    @Test
    public void generateOrthologousStableIdReturnsExistingStableIdentifierInst() throws Exception {
        SimpleInstance existingStableIdentifierInst = createStableIdentifierInst("R-ABC-100001");
        Mockito.when(mockCuratorToolAPI.fetchStableIdentifier("R-ABC-100001")).thenReturn(existingStableIdentifierInst);

        SimpleInstance orthoStableIdentifierInst = stIdGenerator.generateOrthologousStableId(
            new SimpleInstance(), createOriginalInst(createStableIdentifierInst("R-HSA-100001")));

        assertThat(orthoStableIdentifierInst, is(sameInstance(existingStableIdentifierInst)));
        Mockito.verify(mockCuratorToolAPI, Mockito.never()).commit(Mockito.any(SimpleInstance.class));
    }

    @Test
    public void generateOrthologousStableIdCreatesAndCommitsNewStableIdentifierInst() throws Exception {
        SimpleInstance orthoStableIdentifierInst = stIdGenerator.generateOrthologousStableId(
            new SimpleInstance(), createOriginalInst(createStableIdentifierInst("R-HSA-100002")));

        assertThat(orthoStableIdentifierInst.getSchemaClassName(), is(equalTo(ReactomeJavaConstants.StableIdentifier)));
        assertThat(orthoStableIdentifierInst.getAttribute(ReactomeJavaConstants.identifier),
            is(equalTo("R-ABC-100002")));
        assertThat(orthoStableIdentifierInst.getAttribute(ReactomeJavaConstants.identifierVersion), is(equalTo("1")));
        assertThat(orthoStableIdentifierInst.getDisplayName(), is(equalTo("R-ABC-100002.1")));
        Mockito.verify(mockCuratorToolAPI).commit(orthoStableIdentifierInst);
    }

    @Test
    public void paralogStableIdentifierIsGivenCountSuffix() throws Exception {
        SimpleInstance originalInst = createOriginalInst(createStableIdentifierInst("R-HSA-100003"));

        stIdGenerator.generateOrthologousStableId(new SimpleInstance(), originalInst);
        SimpleInstance paralogStableIdentifierInst =
            stIdGenerator.generateOrthologousStableId(new SimpleInstance(), originalInst);

        assertThat(paralogStableIdentifierInst.getAttribute(ReactomeJavaConstants.identifier),
            is(equalTo("R-ABC-100003-2")));
        Mockito.verify(mockCuratorToolAPI).fetchStableIdentifier("R-ABC-100003");
        Mockito.verify(mockCuratorToolAPI).fetchStableIdentifier("R-ABC-100003-2");
    }

    @Test
    public void shellStableIdentifierInstIsInflatedForItsIdentifier() throws Exception {
        SimpleInstance shellStableIdentifierInst = createShellStableIdentifierInst();
        Mockito.when(mockCuratorToolAPI.inflate(shellStableIdentifierInst))
            .thenReturn(createStableIdentifierInst("R-HSA-100004"));

        SimpleInstance orthoStableIdentifierInst = stIdGenerator.generateOrthologousStableId(
            new SimpleInstance(), createOriginalInst(shellStableIdentifierInst));

        assertThat(orthoStableIdentifierInst.getAttribute(ReactomeJavaConstants.identifier),
            is(equalTo("R-ABC-100004")));
    }

    @Test
    public void inferredInstWithStableIdentifierReturnsNull() throws Exception {
        SimpleInstance inferredInst = new SimpleInstance();
        inferredInst.setAttribute(ReactomeJavaConstants.stableIdentifier, createStableIdentifierInst("R-ABC-100005"));

        assertThat(stIdGenerator.generateOrthologousStableId(
            inferredInst, createOriginalInst(createStableIdentifierInst("R-HSA-100005"))), is(nullValue()));
        Mockito.verifyNoInteractions(mockCuratorToolAPI);
    }

    @Test
    public void nullStableIdentifierThrowsRuntimeException() {
        SimpleInstance originalInst = createOriginalInst(null);

        RuntimeException exception = assertThrows(RuntimeException.class,
            () -> stIdGenerator.generateOrthologousStableId(new SimpleInstance(), originalInst));
        assertThat(exception.getMessage(), is(equalTo("No stable identifier instance found for " + originalInst)));
    }

    @Test
    public void stableIdentifierWithoutIdentifierThrowsRuntimeException() {
        SimpleInstance shellStableIdentifierInst = createShellStableIdentifierInst();
        Mockito.when(mockCuratorToolAPI.inflate(shellStableIdentifierInst)).thenReturn(shellStableIdentifierInst);

        RuntimeException exception = assertThrows(RuntimeException.class,
            () -> stIdGenerator.generateOrthologousStableId(
                new SimpleInstance(), createOriginalInst(shellStableIdentifierInst)));
        assertThat(exception.getMessage(), is(equalTo(
            "No identifier found for stable identifier instance " + shellStableIdentifierInst)));
    }

    private SimpleInstance createOriginalInst(SimpleInstance stableIdentifierInst) {
        SimpleInstance originalInst = new SimpleInstance();
        originalInst.setDbId(1L);
        originalInst.setDisplayName("Original instance");
        originalInst.setSchemaClassName(ReactomeJavaConstants.EntityWithAccessionedSequence);
        if (stableIdentifierInst != null) {
            originalInst.setAttribute(ReactomeJavaConstants.stableIdentifier, stableIdentifierInst);
        }
        return originalInst;
    }

    private SimpleInstance createStableIdentifierInst(String identifier) {
        SimpleInstance stableIdentifierInst = new SimpleInstance();
        stableIdentifierInst.setSchemaClassName(ReactomeJavaConstants.StableIdentifier);
        stableIdentifierInst.setDisplayName(identifier + ".1");
        stableIdentifierInst.setAttribute(ReactomeJavaConstants.identifier, identifier);
        stableIdentifierInst.setAttribute(ReactomeJavaConstants.identifierVersion, "1");
        return stableIdentifierInst;
    }

    // A referred instance as curator-tool-ws returns it: dbId, displayName and schemaClass only.
    private SimpleInstance createShellStableIdentifierInst() {
        SimpleInstance shellStableIdentifierInst = new SimpleInstance();
        shellStableIdentifierInst.setDbId(2L);
        shellStableIdentifierInst.setDisplayName("R-HSA-100004.1");
        shellStableIdentifierInst.setSchemaClassName(ReactomeJavaConstants.StableIdentifier);
        return shellStableIdentifierInst;
    }
}