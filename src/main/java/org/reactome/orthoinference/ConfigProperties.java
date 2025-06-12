package org.reactome.orthoinference;

import org.gk.persistence.MySQLAdaptor;
import org.springframework.beans.factory.annotation.Value;
import org.springframework.boot.ApplicationArguments;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.stereotype.Component;

import java.sql.SQLException;
import java.util.Properties;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *         Created 1/9/2024
 */
@Configuration
public class ConfigProperties {
    @Value("${user}")
    private String user;
    @Value("${password}")
    private String password;
    @Value("${currentDbName}")
    private String currentDbName;
    @Value("${previousDbName}")
    private String previousDbName;
    @Value("${host}")
    private String host;
    @Value("${port}")
    private int port;

    @Value("${pathToSpeciesConfig}")
    private String pathToSpeciesConfig;
    @Value("${pathToOrthopairs}")
    private String pathToOrthopairs;

    @Value("${dateOfRelease}")
    private String dateOfRelease;
    @Value("${releaseNumber}")
    private int releaseNumber;
    @Value("${personId}")
    private long personId;

    @Value("${sourceSpeciesCode}")
    private String sourceSpeciesCode;

    public ConfigProperties() {}

    @Bean(name = "currentDBA")
    public MySQLAdaptor getCurrentDBA() throws SQLException {
        return getDBA(currentDbName);
    }

    @Bean(name = "previousDBA")
    public MySQLAdaptor getPreviousDBA() throws SQLException {
        return getDBA(previousDbName);
    }

    @Bean(name = "releaseVersion")
    public int getReleaseVersion() {
        return this.releaseNumber;
    }

    @Bean(name = "pathToSpeciesConfig")
    public String getPathToSpeciesConfig() {
        return this.pathToSpeciesConfig;
    }

    @Bean(name = "pathToOrthopairs")
    public String getPathToOrthopairs() {
        return this.pathToOrthopairs != null ? this.pathToOrthopairs : "orthopairs";
    }

    @Bean
    public SpeciesConfig getSpeciesConfig() {
        String pathToSpeciesConfig = this.pathToSpeciesConfig != null ?
            this.pathToSpeciesConfig : "src/main/resources/Species.json";
        return new SpeciesConfig(pathToSpeciesConfig);
    }

    @Bean(name = "dateOfRelease")
    public String getDateOfRelease() {
        return this.dateOfRelease;
    }

    @Bean(name = "personId")
    public long getPersonId() {
        return this.personId;
    }

    @Bean(name = "sourceSpeciesCode")
    public String getSourceSpeciesCode() {
        return this.sourceSpeciesCode;
    }

    @Bean(name = "targetSpeciesCode")
    public String speciesCode(ApplicationArguments args) {
        final String speciesCodeArgName = "speciesCode";
        if (!args.containsOption(speciesCodeArgName) || args.getOptionValues(speciesCodeArgName).isEmpty()) {
            throw new IllegalArgumentException("Missing required command-line argument: --" + speciesCodeArgName);
        }

        return args.getOptionValues(speciesCodeArgName).get(0);
    }

    private MySQLAdaptor getDBA(String dbName) throws SQLException {
        return new MySQLAdaptor(this.host,dbName,this.user,this.password,this.port);
    }
}
