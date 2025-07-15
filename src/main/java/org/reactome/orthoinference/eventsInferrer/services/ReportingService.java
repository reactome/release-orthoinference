package org.reactome.orthoinference.eventsInferrer.services;

import org.reactome.orthoinference.ConfigProperties;
import org.springframework.stereotype.Component;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

@Component
public class ReportingService {
    private final ConfigProperties configProperties;

    public ReportingService(ConfigProperties configProperties) {
        this.configProperties = configProperties;
    }

    public void outputReport(String species, int eligibleCount, int inferredCount) throws IOException {
        float percentInferred = (float) 100 * inferredCount / eligibleCount;
        // Create file if it doesn't exist
        String reportFilename = "report_ortho_inference_test_reactome_" + configProperties.getReleaseVersion() + ".txt";
        if (!Files.exists(Paths.get(reportFilename))) {
            createNewFile(reportFilename);
        }
        String results = "hsap to " + species + ":\t" + inferredCount + " out of " + eligibleCount +
                " eligible reactions (" + String.format("%.2f", percentInferred) + "%)\n";
        Files.write(Paths.get(reportFilename), results.getBytes(), StandardOpenOption.APPEND);
    }

    private void createNewFile(String filename) throws IOException {
        File file = new File(filename);
        if (file.exists()) {
            file.delete();
        }
        file.createNewFile();
    }

}

