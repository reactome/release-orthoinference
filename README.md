<h2> Orthoinference: Generating Reactome's Computationally Inferred Reactions and Pathways</h2>

This module has been rewritten from Perl to Java.

Additionally, Orthoinference now generates orthologous Stable Identifiers and PathwayDiagrams. <b>StableIdentifier</b> generation represents the updated version of <a href="https://github.com/reactome/Release/blob/master/scripts/release/generate_stable_ids_orthoinference/add_ortho_stable_ids.pl">add_ortho_stable_ids.pl</a> (from the GenerateOrthoStableIds release <a href="https://github.com/reactome/Release/tree/master/scripts/release/generate_stable_ids_orthoinference">step</a>). It now happens during orthoinference whenever a Pathway, ReactionlikeEvent or PhysicalEntity instance is inferred. Orthologous <b>PathwayDiagram</b> generation is completed after all ReactionlikeEvents and Pathways for a species have been inferred, and is based off the <a href="https://github.com/reactome/Release/tree/master/scripts/release/orthodiagrams">OrthoDiagrams</a> release step, which itself utilized the <a href="https://github.com/reactome/CuratorTool/blob/master/src/org/gk/pathwaylayout/PredictedPathwayDiagramGeneratorFromDB.java">PredictedPathwayDiagramGeneratorFromDB</a> method from the CuratorTool repository.

<h3> The Inference Process </h3>

In a nutshell, the inference process follows this workflow:

![Orthoinference Overview Image](https://github.com/reactome/release-orthoinference/blob/develop/OrthoinferenceOverview.png)

For each species, we take all Human <b>ReactionlikeEvents</b> (RlE) instances (<i>Reaction, BlackBoxEvent, Polymerisation, Depolymerisation, FailedReaction</i>) in the `release_current` database that is a copy of the the `slice_current` database after <a href="https://github.com/reactome/release-update-stable-ids">updateStableIds</a> has been run. For each of these RlE instances, there are a few basic rules that must be followed for an inference to be attempted. It must pass a series of <a href="https://github.com/reactome/release-orthoinference/tree/develop/src/main/java/org/reactome/orthoinference/SkipInstanceChecker.java">filters</a> and have <b>at least 1</b> protein instance, determined using the <a href="https://github.com/reactome/release-orthoinference/tree/develop/src/main/java/org/reactome/orthoinference/ProteinCountUtility.java">ProteinCountUtility</a>. 

If the RlE passes all these tests, it is considered <i>eligible</i> for inference. Inference is first attempted on the RlE's <b>input</b> and <b>output</b> attributes, followed by <b>catalyst</b> and <b>regulation</b> inference attempts. <u>If the input or output inferences fail, then the process is terminated for that RlE since they are required components of any ReactionlikeEvent.</u> 
  
During inference, the attribute (input/output/catalyst/regulation) is broken down into its individual <b>PhysicalEntity</b> (PE) instances. These PE instances are each run through the <a href="https://github.com/reactome/release-orthoinference/tree/develop/src/main/java/org/reactome/orthoinference/OrthologousEntityGenerator.java">createOrthoEntity</a> method. This method infers according to PE type: <i>GenomeEncodedEntity/EntityWithAccessionedSequence, Complex/Polymer, EntitySet or SimpleEntity</i>. In cases where the PE itself contains multiple <i>GKInstance</i> attributes (eg: Complexes <i>hasComponent</i> attribute, EntitySets <i>hasMember</i> attribute), these too are run through the createOrthoEntity method and inferred. Through this system, PE instances will be recursively broken down until they reach the <b>EntityWithAccessionedSequence</b> (EWAS) level and are inferred in the <a href="https://github.com/reactome/release-orthoinference/tree/develop/src/main/java/org/reactome/orthoinference/EWASInferrer.java">inferEWAS</a> method.

After all valid ReactionlikeEvents instances have been inferred for a species, the final step is to populate the <b>Pathway</b> instances these RlEs are found in. Pathway inference takes place in the <a href="https://github.com/reactome/release-orthoinference/tree/develop/src/main/java/org/reactome/orthoinference/PathwaysInferrer.java">PathwaysInferrer</a> class, and involves creating the hierarchy structure from the Human Pathway equivalent.

Once the database has been updated with inferred Pathways, orthologous PathwayDiagram generation takes place. This generates the diagrams that are visible in Reactome's Pathway Browser for each species. These diagrams are based off the existing Human diagrams in the database.

<h3> Preparing Orthoinference </h3>

Orthoinference can be run once the <a href="https://github.com/reactome/release-update-stable-ids">UpdateStableIds</a> step has been run. Historically, it had been run following the now-deprecated MyISAM step. Before running the new Java Orthoinference code, there are a few requirements:<br>

- <a href="https://github.com/reactome/release-orthopairs">Orthopair</a> file generation must have been completed.
- The `slice_current` database will need to be dumped and restored to a new `release_current` database.
- Set all values in the `config.properties` file
- `normal_event_skip_list.txt` needs to be placed in `src/main/resources/` folder.

<h4> Setting config.properties </h4>
  
  Create or update the `config.properties` file in the`orthoinference/src/main/resources/` folder, setting the following properties to match the current release:
  
  ```
  ### Sample config.properties file for Orthoinference
  username=mysqlUsername
  password=mysqlPassword
  database=test_reactome_##
  host=localhost
  port=3306
  pathToSpeciesConfig=src/main/resources/Species.json
  pathToOrthopairs=path/to/orthopairs/
  releaseNumber=releaseNumber
  dateOfRelease=yyyy-mm-dd
  personId=reactomePersonInstanceId
  ```
  
  <h4> Orthoinference skiplists </h4>
  
  Historically, the list of ReactionlikeEvents that are skipped during Orthoinference have been manually created by a Curator. As of August 2020, this process has been automated in two ways: 1) A static skiplist that is hard-coded into the orthoinference code and 2) through automated enforcement based on the membership of the ReactionlikeEvent in the Disease TopLevelPathway.
  
 The static skiplist currently consists of the HIV Infection (<b>162906</b>), Influenza Infection (DbId: <b>168255</b>) and Amyloid Fiber Formation (<b>977225</b>; only non-Disease skiplist instance) Pathways. Orthoinference creates a skiplist of all ReactionlikeEvents contained within these Pathways.
  
  The automated skiplist is only focused on skipping inference for instances that are children of the Disease TopLevelPathway. If a ReactionlikeEvent only exists as a child of the Disease pathway, than inference will be skipped. In cases where a ReactionlikeEvent is a member of the Disease AND another TopLevelPathway, then <b>Reaction</b> inference proceeds as normal. During <b>Pathway</b> inference, the inference of the non-Disease pathways are allowed while the Disease Pathway inference (and its children) are suppressed for the instance. 
  
  <h3> Running Orthoinference </h3>
  
  <b>Note</b>: For the orthoinference steps and QA, it is recommended that <b>Java 8</b> be used and <b>mySQL 5.5</b> or <b>5.7</b> be used.
  
  Once all <a href="#-preparing-orthoinference-">prerequisites</a> have been completed, running the <a href="https://github.com/reactome/release-orthoinference/blob/develop/runOrthoinference.sh">runOrthoinference.sh</a> script will begin the process. This bash script performs a git pull to update the repo with any changes that may have happened between releases. It then builds the orthoinference jar file with all dependencies and then executes it for each species that will be projected to.
  
  <b>Note</b>: To run orthoinference on particular species, modify the 'allSpecies' array in <a href="https://github.com/reactome/release-orthoinference/blob/develop/runOrthoinference.sh">runOrthoinference.sh</a> so that it only contains the species you wish to project too. Alternatively, if the jar file has been built and only one species needs to be inferred, run the following command:<br> 
`java -jar target/orthoinference-[version]-jar-with-dependencies.jar [speciesCode]`
- Replace '[speciesCode]' with the 4 letter species code corresponding to the species you wish to infer too
- Replace '[version]' with the build version for the jar (e.g. 0.0.1-SNAPSHOT)
- Orthoinference benefits from an increased memory heap, which can be modified with the `-Xmx####m` tag before `-jar`.
  
 During orthoinference, many files are produced:
 
 - Log files in the `logs/` folder provide information pertaining to each inference attempt and is useful for tracing errors.
   - They are organized by time stamp.
 - `eligible_(speciesCode)_75.txt` lists all ReactionlikeEvents that can be inferred. This should be the same for all species.
   - The 75 refers to the percent of distinct proteins that must exist in <b>Complex/Polymer</b> instances for an inference attempt to continue. It is a holdover name from Perl Orthoinference.
 - `inferred_(speciesCode)_75.txt` lists all ReactionlikeEvents that were successfully inferred for the species.
 - `report_ortho_inference_test_reactome_##.txt` shows the percentage of inferences that were successful for each species.
 
 Once the Java code has been finished, verify that all `eligible_(speciesCode)_75.txt` files have the same number of lines. If the line counts are different, something likely went wrong during inference and will need to be investigated.
 
 Finally, once you have verified that orthoinference seems to have run correctly, run <a href="https://github.com/reactome/Release/blob/master/scripts/release/updateDisplayName.pl">updateDisplayName.pl</a>. 
 
<h3> Verifying Orthoinference </h3>

Once all scripts in the previous step have been run there is a QA process that should be followed. Orthoinference is a foundational step for the Reactome release pipeline, and ensuring that this process worked as expected will save much time later in the Release process if anything erroneous happened. 

<b> Recommended QA </b><br>

1) Compare line counts of the `eligible_(speciesCode)_75.txt` and `inferred_(speciesCode)_75.txt` files to the previous release. Are they relatively close to each other? If any are significantly smaller for the eligible files of a particular species, perhaps check the <b>Orthopair</b> files that correspond to the species of interest. Are the eligible and/or inferred line counts considerably smaller for all species? Something may have gone wrong during the inference process itself. Check log files to see if anything obvious jumps out. Otherwise, more extensive troubleshooting with the script itself will be required.

Next, we want to make sure that the new `release_current` database can be imported to the <b>graphDb</b> in Neo4j and that it reports acceptable graph QA numbers. <br><br><b>It is recommended that the following steps be run on your workstation.</b>

2) Run the <a href="https://github.com/reactome/graph-importer">graph-importer</a> module. This should take some time (~30 minutes) and will output the results from <a href="https://github.com/reactome/database-checker">database-checker</a> as well as the imported graphDb in `target/graph.db/`.

    -  The <b>database-checker</b> results should resemble the following:
    ```
    The consistency check finished reporting 13 as follows:
             1 entries for (Taxon, crossReference)

    Valid schema class check finish reporting 1 entry
    ```
    The database-checker  module just looks for any attributes of an instance that are <i>required</i> (as determined by the data model) and are not filled. Small numbers reported are OK but any newly reported entries should be investigated. 
    
3) Finally, running the <a href="https://github.com/reactome/graph-qa">graph-qa</a> step will check a series of graphDB QA items and rank them by urgency. To run graph-qa, you will need to have an instance of neo4j running with the imported graph DB. To quickly get neo4j installed and running, a docker container is recommended. The following command can be used to get it running quickly and painlessly:<br><br>
` docker run -p 7687:7687 -p 7474:7474 -e NEO4J_dbms_allow__upgrade=true -e NEO4J_dbms_allowFormatMigration=true -e NEO4J_dbms_memory_heap_max__size=XXG -v $(pwd)/target/graph.db:/var/lib/neo4j/data/databases/graph.db neo4j:3.4.9` 
    -  Make sure that the location of the `graph.db/` folder is properly specified in the last argument 
    -  Adjust the `NEO4J_dbms_memory_heap_max__size` argument so that it is  appropriate for your computer/server. 
  
  To verify that the graphDb has been properly populated, open `localhost:7474` (username and password default is `neo4j`), and click on the database icon at the top left. A panel titled <i>Database Information</i> should open up and display all nodes in the <a href="https://reactome.org/content/schema/">Data Model</a>. If none of this appears, chances are the neo4j instance is not using the imported graphDB. Verify the `graph.db/` folder is in fact populated. 

To verify that the graphDb has been properly populated, open `localhost:7474` in your browser once the docker container is built (username and password default is <b>neo4j</b>), and click on the <i>database icon</i> at the top left. A panel titled <i>Database Information</i> should open up and display all nodes in the <a href="https://reactome.org/content/schema">Data Model</a>. If none of this appears, chances are the neo4j instance is not using the imported graphDB. Verify the `graph.db/` folder is in fact populated and make sure the location of the `graph.db/` folder is properly specified in the docker command's last argument.

Once you have verified that the graphDb is populated, graph-qa can be run. Build the jar file using `mvn clean compile package` and then follow the instructions at the <a href="https://github.com/reactome/graph-qa">repo</a>. After `graph-qa` has been run, it will output the number of offending instances for each QA category, ranked by <b>urgency/priority</b>. These results can be compared with the QA results from the previous release.  QA reports can be found <a href="https://drive.google.com/drive/folders/0Byf6VJ-jol9OVnZiWjlzVlJOSGc">here</a>.

4) Once you are satisfied with the results from each of these steps, send the <b>graph-qa</b> and <b>database-checker</b> results to the Curator overseeing release. <a href="https://github.com/reactome/database-checker">Database-checker</a> can be re-run fairly painlessly by following the instructions on its GitHub page. If the Curator is satisfied with the QA results, you can move onto one of the next steps of release. At this point in the release pipeline that could be <a href="https://github.com/reactome/release-update-dois">UpdateDOIs</a> or <a href="https://github.com/reactome/add-links">AddLinks</a>.

<b> Additional Orthoinference troubleshooting</b>

Below are suggestions for further checking that the Orthoinference process was correctly run. They don't provide as much information or have a specific SOP though, making them optional QA processes. 

-  The <b>WebELVTool</b> found in the <a href="https://github.com/reactome/CuratorTool">CuratorTool</a> repo can be used to check the validity of the `release_current` database. The WebELVTool jar is built by running an <b>ant</b> build on the <a href="https://github.com/reactome/CuratorTool/blob/master/ant/WebELVDiagram.xml">WebELVDiagram.xml</a> file. The jar file will appear in the `WebELVTool/` folder above the `CuratorTool/` folder. To run the program, navigate to the directory with the jar file and run:<br> `java -jar webELVTool.jar (host) (release_current) (mysql_username) (mysql_password) (mysql_port) (reactome_author_ID)`<br>The output should be many lines of 'Predicting pathway' or 'Working on pathway'. If the script runs successfully, then it is implied `release_current` is appropriately populated.
-  The <b>CuratorTool</b> program (different from the CuratorTool repository mentioned above), can be downloaded from the Reactome website <a href="https://reactome.org/download-data/reactome-curator-tool">here</a> can also be used to explore the inferred data. There isn't a recommended set of tests or checks, but familiarity with the CuratorTool can be quite useful for determining if instances are populated correctly. 
