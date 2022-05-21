// This Jenkinsfile is used by Jenkins to run the Orthoinference step of Reactome's release.
// It requires that the Orthopairs and UpdateStableIdentifiers steps have been run successfully before it can be run.

import org.reactome.release.jenkins.utilities.Utilities

// Shared library maintained at 'release-jenkins-utils' repository.
def utils = new Utilities()
pipeline{
	agent any

	stages{
		// This stage checks that upstream projects Orthopairs and UpdateStableIdentifier, were run successfully for their last build.
		stage('Check if Orthopairs and UpdateStableIdentifiers builds succeeded'){
			steps{
				script{
					utils.checkUpstreamBuildsSucceeded("Orthopairs")
					utils.checkUpstreamBuildsSucceeded("ConfirmReleaseConfigs")
				}
			}
		}
		stage('Setup: Download Orthopairs files from S3 bucket'){
			steps{
				script{
					def releaseVersion = utils.getReleaseVersion()
					sh "mkdir -p orthopairs"
					sh "aws s3 --no-progress cp --recursive ${env.S3_RELEASE_DIRECTORY_URL}/${releaseVersion}/orthopairs/data/orthopairs/ ./orthopairs/"
					sh "gunzip orthopairs/*"
				}
			}
		}
		// This stage backs up the release current database before it is modified.
		stage('Setup: Backup release_current'){
			steps{
				script{
					withCredentials([usernamePassword(credentialsId: 'mySQLUsernamePassword', passwordVariable: 'pass', usernameVariable: 'user')]){
						utils.takeDatabaseDumpAndGzip("${env.RELEASE_CURRENT_DB}", "orthoinference", "before", "${env.RELEASE_SERVER}")
					}
				}
			}
		}
		// This stage builds the jar file using maven. It also runs the Main orthoinference process as a 'sub-stage'.
		// This was due to a restriction in iterating over a list of species names. To iterate, you need to first have a 'stage > steps > script' hierarchy.
		// At the script level, you can iterate over a list and then create new stages from this iteration. The choice was between an empty stage or to do a sub-stage.
		stage('Setup: Run Orthoinference on species list'){
			steps{
				// This script block executes the main orthoinference code one species at a time.
				// It takes all Human Reaction instances in the database and attempts to project each Reaction to each species by
				// stripping them down to the reaction's constituent proteins, checks if the protein homolog exists for that species, and infers it in Reactome's data model.
				// If enough proteins (>= 75%) are inferrable in a Reaction, then it is created and stored in the database for this release. This is done from scratch each time.
				script{
					speciesList = ['mmus', 'rnor', 'cfam', 'btau', 'sscr', 'drer', 'xtro', 'ggal', 'dmel', 'cele', 'ddis', 'spom', 'scer', 'pfal']
					for (species in speciesList) {
						stage("Main: Infer ${species}"){
							script{
								withCredentials([file(credentialsId: 'Config', variable: 'ConfigFile')]){
									// Changes name of output log files to include 4-letter species name, for easier file management.
									sh "git checkout src/main/resources/log4j2.xml"
									sh "sed -i -e 's/OrthoInference/${species}-OrthoInference/g' src/main/resources/log4j2.xml"
									utils.buildJarFile()
									sh "java -Xmx${env.JAVA_MEM_MAX}m -jar target/orthoinference-*-jar-with-dependencies.jar $ConfigFile ${species}"
								}
							}
						}
					}
				}
			}
		}
		// This stage sorts the order of the species in the output report file and then symlinks it to the website_files_update directory.
		stage('Post: Sort output report & create symlink'){
			steps{
				script{
					def releaseVersion = utils.getReleaseVersion()
					sh "./formatOrthoinferenceReport.sh --release ${releaseVersion}"
					def inferenceReportFilename = "report_ortho_inference_test_reactome_${releaseVersion}_sorted.txt"
					sh "cp ${inferenceReportFilename} ${env.WEBSITE_FILES_UPDATE_ABS_PATH}/"
					dir("${env.WEBSITE_FILES_UPDATE_ABS_PATH}"){
						// Creates hard-link of sorted report_ortho_inference file, in website_files_update folder.
						// This allows the generically named file to be committed to git so that we can track it over releases.
						sh "ln -f ${inferenceReportFilename} report_ortho_inference.txt"
					}
				}
			}
		}
		// This stage downloads the previous releases orthoinference files (eligible, inferred), and outputs line count differences between them.
		stage('Post: Orthoinference file line counts') {
		    steps{
		        script{
		            def releaseVersion = utils.getReleaseVersion()
		            def previousReleaseVersion = utils.getPreviousReleaseVersion()
		            def orthoinferencesDir = "orthoinferences"
		            def currentDir = pwd()
		            
			    // Create the 'orthoinferences' and 'previousReleaseVersion' directories
		            sh "mkdir -p ${orthoinferencesDir} ${previousReleaseVersion}"
		            sh "mv eligible* inferred* ${orthoinferencesDir}/"
		            sh "aws s3 --recursive --no-progress cp s3://reactome/private/releases/${previousReleaseVersion}/orthoinference/data/orthoinferences/ ${previousReleaseVersion}/"
		            sh "gunzip ${previousReleaseVersion}/*"
	                    utils.outputLineCountsOfFilesBetweenFolders("$orthoinferencesDir", "$previousReleaseVersion", "$currentDir")
	                    sh "rm -r ${previousReleaseVersion}"
		        }
		    }
		}
		// This stage backs up the release_current database after it is modified.
		stage('Post: Backup DB'){
			steps{
				script{
					withCredentials([usernamePassword(credentialsId: 'mySQLUsernamePassword', passwordVariable: 'pass', usernameVariable: 'user')]){
						utils.takeDatabaseDumpAndGzip("${env.RELEASE_CURRENT_DB}", "orthoinference", "after", "${env.RELEASE_SERVER}")
					}
				}
			}
		}
		// This stage generates the graph database using the graph-importer module, and replaces the current graph db with it.
		stage('Post: Generate Graph Database'){
			steps{
				script{
					// Gets a copy of 'changeGraphDatabase', which Jenkins can execute as sudo. Changes permissions of file to user read/write only.
					utils.cloneOrUpdateLocalRepo("release-jenkins-utils")
					sh "cp -f release-jenkins-utils/scripts/changeGraphDatabase.sh ${env.JENKINS_HOME_PATH}"
					sh "chmod 700 ${env.JENKINS_HOME_PATH}/changeGraphDatabase.sh"
					utils.cloneOrUpdateLocalRepo("graph-importer")
					
					dir("graph-importer"){
						utils.buildJarFileWithPackage()
						// This generates the graph database.
						withCredentials([usernamePassword(credentialsId: 'mySQLUsernamePassword', passwordVariable: 'pass', usernameVariable: 'user')]){
							sh "java -jar target/GraphImporter-exec.jar --name ${env.RELEASE_CURRENT_DB} --user $user --password $pass --neo4j /tmp/graph.db"
							sh "sudo service tomcat9 stop"
							sh "sudo service neo4j stop"
							// This static script adjusts permissions of the graph.db folder and moves it to /var/lib/neo4j/data/databases/.
							sh "sudo bash ${env.JENKINS_HOME_PATH}/changeGraphDatabase.sh"
							sh "sudo service neo4j start"
							sh "sudo service tomcat9 start"
							sh "rm ${env.JENKINS_HOME_PATH}/changeGraphDatabase.sh"
						}
					}
				}
			}			
		}
		// This stage runs the graph-qa script that will be emailed to Curation.
		stage('Post: Run graph-qa'){
			steps{
				script{
					utils.cloneOrUpdateLocalRepo("graph-qa")
					dir("graph-qa"){
						utils.buildJarFile()
						withCredentials([usernamePassword(credentialsId: 'neo4jUsernamePassword', passwordVariable: 'pass', usernameVariable: 'user')]){
							sh "java -jar target/graph-qa-jar-with-dependencies.jar -u $user -p  $pass --verbose"
						}
					}
				}
			}
		}
		// This stage emails the contents of the GraphQA_Summary_vXX.csv that is generated by graph-qa for review, as well as the GraphQA_Summary_vXX-1 (previous release). 
		stage('Post: Email graph-qa output'){
			steps{
				script{
					def releaseVersion = utils.getReleaseVersion()
					def previousReleaseVersion = utils.getPreviousReleaseVersion()
					def prevGraphQAFileName = "GraphQA_Summary_v${previousReleaseVersion}.csv"
					def currentGraphQAFileName = "GraphQA_Summary_v${releaseVersion}.csv"
					def s3PathPrevGraphQA = "${env.S3_RELEASE_DIRECTORY_URL}/${previousReleaseVersion}/orthoinference/logs/${prevGraphQAFileName}.gz"
					// Get previous release graph-qa output
					sh "aws s3 cp ${s3PathPrevGraphQA} ."
					sh "gunzip ${prevGraphQAFileName}.gz"
					// Email the graph-qa outputs from the current and previous releases
					def emailSubject = "Orthoinference graph-qa for v${releaseVersion}"
					def emailBody = "Hello,\n\nThis is an automated message from Jenkins regarding an update for v${releaseVersion}. The Orthoinference step has finished running. Attached to this email should be the summary report output by graph-qa for both v${releaseVersion} and v${previousReleaseVersion}. Please compare these and confirm if they look appropriate with the developer running Release. \n\nThanks!"
					def emailAttachments = "graph-qa/reports/${currentGraphQAFileName}, ${prevGraphQAFileName}"
					utils.sendEmailWithAttachment("$emailSubject", "$emailBody", "$emailAttachments")
					
					sh "rm ${prevGraphQAFileName}"
				}
			}
		}
		// All databases, logs, reports, and data files generated by this step are compressed before moving them to the Reactome S3 bucket. 
		// All remaining files/folders are then deleted that are not a part of the release-orthoinference repository.
		stage('Post: Archive Outputs'){
			steps{
				script{
				    def releaseVersion = utils.getReleaseVersion()
				    def dataFiles = ["orthoinferences", "report_ortho_inference_test_reactome_${releaseVersion}*.txt"]
					// Additional log files from post-step QA need to be pulled in
					def logFiles = ["graph-importer/logs/*", "graph-qa/logs/*", "graph-qa/reports/*"]
					// This folder is utilized for post-step QA. Jenkins creates multiple temporary directories
					// cloning and checking out repositories, which is why the wildcard is added.
					def foldersToDelete = ["orthopairs", "release-jenkins-utils*", "graph-importer*", "graph-qa*"]
					utils.cleanUpAndArchiveBuildFiles("orthoinference", dataFiles, logFiles, foldersToDelete)
				}
			}
		}
	}
}
