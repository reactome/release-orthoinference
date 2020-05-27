import groovy.json.JsonSlurper
// This Jenkinsfile is used by Jenkins to run the Orthoinference step of Reactome's release.
// It requires that the Orthopairs and UpdateStableIdentifiers steps have been run successfully before it can be run.
def currentRelease
def previousRelease
pipeline{
	agent any

	stages{
		// This stage checks that upstream projects Orthopairs and UpdateStableIdentifier, were run successfully for their last build.
		stage('Check if Orthopairs and UpdateStableIdentifiers builds succeeded'){
			steps{
				script{
					currentRelease = (pwd() =~ /(\d+)\//)[0][1];
					previousRelease = (pwd() =~ /(\d+)\//)[0][1].toInteger() - 1;
					// This queries the Jenkins API to confirm that the most recent builds of Orthopairs and UpdateStableIdentifiers were successful.
					checkUpstreamBuildsSucceeded("Orthopairs", "$currentRelease")
					checkUpstreamBuildsSucceeded("UpdateStableIdentifiers", "$currentRelease")
				}
			}
		}
		/*
		// Orthoinference utilizes a skiplist of Reaction DbIds to prevent particular reactions from being inferred.
		stage('User Input Required: Confirm skiplist uploaded'){
			steps{
				script{
					def userInput = input(
						id: 'userInput', message: "Has the Orthoinference skiplist been uploaded as a locally scoped credential? (yes/no)",
						parameters: [
							[$class: 'TextParameterDefinition', defaultValue: '', description: 'Confirmation of orthoinference skiplist', name: 'response']
						])

					if (userInput.toLowerCase().startsWith("y")) {
						echo("Proceeding with Orthoinference step.")
					} else {
						error("Please upload the skiplist to Jenkins>Releases>${currentRelease}>Credentials>orthoinferenceSkipList. You should have received this skiplist from Curation.")
					}
				}
			}
		}
		stage('Setup: Download Orthopairs files from S3 bucket'){
			steps{
				script{
					sh "mkdir -p orthopairs"
					sh "aws s3 --no-progress cp  --recursive ${env.S3_RELEASE_DIRECTORY_URL}/${currentRelease}/orthopairs/data/orthopairs/ ./orthopairs/"
					sh "gunzip orthopairs/*"
				}
			}
		}
		// This stage backs up the release current database before it is modified.
		stage('Setup: Backup release_current'){
			steps{
				script{
					withCredentials([usernamePassword(credentialsId: 'mySQLUsernamePassword', passwordVariable: 'pass', usernameVariable: 'user')]){
						def release_current_before_orthoinference_dump = "${env.RELEASE_CURRENT}_${currentRelease}_before_orthoinference.dump"
						sh "mysqldump -u$user -p$pass ${env.RELEASE_CURRENT} > $release_current_before_orthoinference_dump"
						sh "gzip -f $release_current_before_orthoinference_dump"
					}
				}
			}
		}
		// This stage builds the jar file using maven. It also runs the Main orthoinference process as a 'sub-stage'.
		// This was due to a restriction in iterating over a list of species names. To iterate, you need to first have a 'stage > steps > script' hierarchy.
		// At the script level, you can iterate over a list and then create new stages from this iteration. The choice was between an empty stage or to do a sub-stage.
		stage('Setup: Build jar file'){
			steps{
				script{
					sh "mvn clean compile assembly:single"
				}
				// This script block executes the main orthoinference code one species at a time.
				// It takes all Human Reaction instances in the database and attempts to project each Reaction to each species by
				// stripping them down to the reaction's constituent proteins, checks if the protein homolog exists for that species, and infers it in Reactome's data model.
				// If enough proteins (>= 75%) are inferrable in a Reaction, then it is created and stored in the database for this release. This is done from scratch each time.
				script{
					withCredentials([file(credentialsId: 'orthoinferenceSkipList', variable: 'skipListFile')]) {
						sh "cp -f $skipListFile normal_event_skip_list.txt"
						speciesList = ['mmus', 'rnor', 'cfam', 'btau', 'sscr', 'drer', 'xtro', 'ggal', 'dmel', 'cele', 'ddis', 'spom', 'scer', 'pfal']
						for (species in speciesList) {
							stage("Main: Infer ${species}"){
								script{
									withCredentials([file(credentialsId: 'Config', variable: 'ConfigFile')]){
										sh "java -Xmx${env.JAVA_MEM_MAX}m -jar target/orthoinference-*-jar-with-dependencies.jar $ConfigFile ${species}"
									}
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
					sh "./formatOrthoinferenceReport.sh --release ${currentRelease}"
					def cwd = pwd()
					sh "ln -sf ${cwd}/report_ortho_inference_test_reactome_${currentRelease}_sorted.txt ${env.WEBSITE_FILES_UPDATE_ABS_PATH}/report_ortho_inference.txt"
				}
			}
		}
		// This stage backs up the release_current database after it is modified.
		stage('Post: Backup DB'){
			steps{
				script{
					withCredentials([usernamePassword(credentialsId: 'mySQLUsernamePassword', passwordVariable: 'pass', usernameVariable: 'user')]){
						def release_current_after_orthoinference_dump = "${env.RELEASE_CURRENT}_${currentRelease}_after_orthoinference.dump"
						sh "mysqldump -u$user -p$pass ${env.RELEASE_CURRENT} > $release_current_after_orthoinference_dump"
						sh "gzip -f $release_current_after_orthoinference_dump"
					}
				}
			}
		}
		*/
		stage('Post: Generate Graph Database'){
			steps{
				script{
					sh "git clone https://github.com/reactome/graph-importer"
					dir("graph-importer"){
						sh "mvn clean compile assembly:single"
						withCredentials([usernamePassword(credentialsId: 'mySQLUsernamePassword', passwordVariable: 'pass', usernameVariable: 'user')]){
							sh "sudo su neo4j"
							sh "java -jar target/GraphImporter-jar-with-dependencies.jar --name ${env.RELEASE_CURRENT} --user $user --password $pass --neo4j ${env.NEO4J_DIR}/${currentRelease}.graph.db"
						}
					}
				}
			}			
		}
		/*
		// This stage archives all logs and database backups produced by Orthoinference. It also archives the eligible/inferred files produced by orthoinference.
		stage('Post: Archive logs and backups'){
			steps{
				script{
					sh "mkdir -p archive/${currentRelease}/logs"
					sh "mv --backup=numbered *_${currentRelease}_*.dump.gz archive/${currentRelease}/"
					sh "gzip logs/*"
					sh "mv logs/* archive/${currentRelease}/logs/"
					sh "mkdir -p ${currentRelease}"
					sh "gzip -f *.txt"
					sh "mv *.txt.gz ${currentRelease}/"
				}
			}
		}
		*/
	}
}

// Utility function that checks upstream builds of this project were successfully built.
def checkUpstreamBuildsSucceeded(String stepName, String currentRelease) {
	def statusUrl = httpRequest authentication: 'jenkinsKey', validResponseCodes: "${env.VALID_RESPONSE_CODES}", url: "${env.JENKINS_JOB_URL}/job/$currentRelease/job/$stepName/lastBuild/api/json"
	if (statusUrl.getStatus() == 404) {
		error("$stepName has not yet been run. Please complete a successful build.")
	} else {
		def statusJson = new JsonSlurper().parseText(statusUrl.getContent())
		if(statusJson['result'] != "SUCCESS"){
			error("Most recent $stepName build status: " + statusJson['result'] + ". Please complete a successful build.")
		}
	}
}
