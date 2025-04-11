package org.reactome.orthoinference;

import org.gk.persistence.MySQLAdaptor;

import java.sql.SQLException;
import java.util.Properties;

/**
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 *         Created 1/9/2024
 */
public class ConfigProperties {
    private Properties props;

    public ConfigProperties(Properties props) {
        this.props = props;
    }

    public MySQLAdaptor getCurrentDBA() throws SQLException {
        String dbName = props.getProperty("release_current.name");
        return getDBA(props, dbName);
    }

    public MySQLAdaptor getPreviousDBA() throws SQLException {
        String previousDbName = props.getProperty("release_previous.name");
        return getDBA(props, previousDbName);
    }

    public String getReleaseVersion() {
        return props.getProperty("releaseNumber");
    }

    public String getPathToOrthopairs() {
        return props.getProperty("pathToOrthopairs", "orthopairs");
    }

    public SpeciesConfig getSpeciesConfig() {
        String pathToSpeciesConfig = props.getProperty("pathToSpeciesConfig", "src/main/resources/Species.json");
        return new SpeciesConfig(pathToSpeciesConfig);
    }

    public String getDateOfRelease() {
        return props.getProperty("dateOfRelease");
    }

    public int getPersonId() {
        return Integer.valueOf(props.getProperty("personId"));
    }

    private MySQLAdaptor getDBA(Properties props, String dbName) throws SQLException {
        String username = props.getProperty("release.database.user");
        String password = props.getProperty("release.database.password");
        String host = props.getProperty("release.database.host");
        int port = Integer.valueOf(props.getProperty("release.database.port"));

        return new MySQLAdaptor(host,username,dbName,password,port);
    }

    private Properties getProps() {
        return this.props;
    }
}
