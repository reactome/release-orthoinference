package org.reactome.orthoinference;

import org.gk.persistence.MySQLAdaptor;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.boot.context.properties.EnableConfigurationProperties;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.Lazy;
import org.springframework.boot.ApplicationArguments;

import java.sql.SQLException;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *      Created 1/9/2024
 */
@Configuration
@EnableConfigurationProperties(ConfigProperties.OrthoinferenceProps.class)
public class ConfigProperties {

    private final OrthoinferenceProps props;

    public ConfigProperties(OrthoinferenceProps props) {
        this.props = props;
    }

    @Lazy
    @Bean(name = "currentDBA")
    public MySQLAdaptor getCurrentDBA() throws SQLException {
        return getDBA(props.getCurrentDbName());
    }

    public MySQLAdaptor getPreviousDBA() throws SQLException {
        return getDBA(props.getPreviousDbName());
    }

    public int getReleaseVersion() {
        return props.getReleaseNumber();
    }

    @Bean(name = "pathToSpeciesConfig")
    public String getPathToSpeciesConfig() {
        return props.getPathToSpeciesConfig() != null ? props.getPathToSpeciesConfig() : "src/main/resources/Species.json";
    }

    public String getPathToOrthopairs() {
        return props.getPathToOrthopairs() != null ? props.getPathToOrthopairs() : "orthopairs";
    }

    public SpeciesConfig getSpeciesConfig() {
        return new SpeciesConfig(getPathToSpeciesConfig());
    }

    public String getDateOfRelease() {
        return props.getDateOfRelease();
    }

    public long getPersonId() {
        return props.getPersonId();
    }

    public String getSourceSpeciesCode() {
        return props.getSourceSpeciesCode();
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
        return new MySQLAdaptor(props.getHost(), dbName, props.getUser(), props.getPassword(), props.getPort());
    }

    @ConfigurationProperties(prefix = "orthoinference")
    public static class OrthoinferenceProps {
        private String user;
        private String password;
        private String currentDbName;
        private String previousDbName;
        private String host;
        private int port;
        private String pathToSpeciesConfig;
        private String pathToOrthopairs;
        private String dateOfRelease;
        private int releaseNumber;
        private long personId;
        private String sourceSpeciesCode;

        public String getUser() {
            return this.user;
        }

        public String getPassword() {
            return this.password;
        }

        public String getCurrentDbName() {
            return this.currentDbName;
        }

        public String getPreviousDbName() {
            return this.previousDbName;
        }

        public String getHost() {
            return this.host;
        }

        public void setHost(String host) {
            this.host = host;
        }

        public int getPort() {
            return this.port;
        }

        public String getPathToSpeciesConfig() {
            return this.pathToSpeciesConfig;
        }

        public String getPathToOrthopairs() {
            return this.pathToOrthopairs;
        }

        public String getDateOfRelease() {
            return this.dateOfRelease;
        }

        public int getReleaseNumber() {
            return this.releaseNumber;
        }

        public long getPersonId() {
            return this.personId;
        }

        public String getSourceSpeciesCode() {
            return this.sourceSpeciesCode;
        }
    }
}
