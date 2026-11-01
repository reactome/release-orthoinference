package org.reactome.orthoinference;

import org.reactome.curation.CuratorToolWsApplication;
import org.reactome.curation.controller.CurationController;
import org.springframework.boot.WebApplicationType;
import org.springframework.boot.builder.SpringApplicationBuilder;
import org.springframework.context.ConfigurableApplicationContext;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.core.io.Resource;
import org.springframework.core.io.UrlResource;
import org.springframework.core.env.MapPropertySource;
import org.springframework.core.io.support.ResourcePropertySource;
import org.springframework.web.context.support.StandardServletEnvironment;

import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

/**
 * Starts curator-tool-ws in an application context of its own, separate from orthoinference's, and exposes its
 * CurationController for CuratorToolAPI. curator-tool-ws is a full servlet application with its own properties,
 * Neo4j, security and AspectJ set-up, so it is kept apart rather than component-scanned into this context.
 *
 * @author Joel Weiser (joel.weiser@oicr.on.ca)
 * Created 10/14/2026
 */
@Configuration
public class CuratorToolConfig {

    // The following code is copied directly from the slicing tool project.
    @Bean(name = "curatorToolContext", destroyMethod = "close")
    public ConfigurableApplicationContext curatorToolContext() throws IOException {
        // curator-tool-ws's bundled application.properties forces DEBUG for these loggers and binds the
        // HTTP connector to 9090. System properties outrank a classpath application.properties in Spring
        // Boot's precedence order, so these settings take hold for the batch run without editing
        // curator-tool-ws. (SpringApplicationBuilder.properties(...) are default/lowest precedence and
        // would NOT override application.properties.)
        System.setProperty("logging.level.org.springframework.data.neo4j", "WARN");
        System.setProperty("logging.level.org.springframework.security", "WARN");
        // Disable the HTTP server; the full servlet context is kept for correct AspectJ wiring.
        System.setProperty("server.port", "-1");

        // Orthoinference's own application.properties comes first on the classpath and hides curator-tool-ws's file
        // of the same name, so Spring Boot never loads curator-tool-ws's settings (its Neo4j credentials among
        // them). The file is read from curator-tool-ws's own jar instead and added below the system properties
        // above, so those still take precedence.
        StandardServletEnvironment environment = new StandardServletEnvironment();
        // This context shares the JVM's MBean server with orthoinference's, so if JMX is switched on (as an IDE's
        // Spring Boot launcher does with -D flags) both register the same SpringApplication admin MBean and the
        // second fails with InstanceAlreadyExistsException. JMX is turned off for this context alone, ahead of the
        // system properties so those flags cannot turn it back on.
        Map<String, Object> jmxDisabled = new HashMap<>();
        jmxDisabled.put("spring.application.admin.enabled", "false");
        jmxDisabled.put("spring.jmx.enabled", "false");
        environment.getPropertySources().addFirst(new MapPropertySource("curator-tool-ws JMX disabled", jmxDisabled));
        environment.getPropertySources().addLast(
            new ResourcePropertySource("curator-tool-ws application.properties", curatorToolProperties()));

        return new SpringApplicationBuilder(CuratorToolWsApplication.class)
            .environment(environment)
            .web(WebApplicationType.SERVLET)
            .run();
    }

    private Resource curatorToolProperties() throws MalformedURLException {
        // A jar when curator-tool-ws comes from the Maven repository; a classes directory when it is opened as a
        // module alongside this project in an IDE.
        URL curatorToolLocation = CuratorToolWsApplication.class.getProtectionDomain().getCodeSource().getLocation();
        return curatorToolLocation.getPath().endsWith(".jar") ?
            new UrlResource("jar:" + curatorToolLocation + "!/application.properties") :
            new UrlResource(new URL(curatorToolLocation, "application.properties"));
    }

    // The CurationController is handed to CuratorToolAPI rather than registered as a bean of its own: this context
    // would re-run @Autowired injection on it, and its dependencies (CurationService, QAService, ...) exist only in
    // curator-tool-ws's context. curatorToolContext() is called directly, not injected, because a
    // ConfigurableApplicationContext parameter resolves to this application's own context ahead of any bean.
    @Bean
    public CuratorToolAPI curatorToolAPI(ConfigProperties configProperties) throws IOException {
        return new CuratorToolAPI(curatorToolContext().getBean(CurationController.class), configProperties);
    }
}