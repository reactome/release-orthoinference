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
					currentRelease = (pwd() =~ /Releases\/(\d+)\//)[0][1];
					previousRelease = (pwd() =~ /Releases\/(\d+)\//)[0][1].toInteger() - 1;
					// This queries the Jenkins API to confirm that the most recent builds of Orthopairs and UpdateStableIdentifiers were successful.
//					checkUpstreamBuildsSucceeded("Orthopairs", "$currentRelease")
//					checkUpstreamBuildsSucceeded("Relational-Database-Updates/job/UpdateStableIdentifiers", "$currentRelease")
				}
			}
		}
		/*
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
		// This stage generates the graph database using the graph-importer module, and replaces the current graph db with it.
		stage('Post: Generate Graph Database'){
			steps{
				script{
					cloneOrPullGitRepo("release-jenkins-utils")
					sh "cp -f release-jenkins-utils/scripts/changeGraphDatabase.sh ${env.JENKINS_HOME_PATH}"
					sh "chmod 700 ${env.JENKINS_HOME_PATH}changeGraphDatabase.sh"
					cloneOrPullGitRepo("graph-importer")
					dir("graph-importer"){
						sh "mvn clean compile assembly:single"
						withCredentials([usernamePassword(credentialsId: 'mySQLUsernamePassword', passwordVariable: 'pass', usernameVariable: 'user')]){
							sh "java -jar target/GraphImporter-jar-with-dependencies.jar --name ${env.RELEASE_CURRENT} --user $user --password $pass --neo4j /tmp/graph.db"
							sh "sudo service tomcat7 stop"
							sh "sudo service neo4j stop"
							// This static script adjusts permissions of the graph.db folder and moves it to /var/lib/neo4j/data/databases/.
							sh "sudo bash ${env.JENKINS_HOME_PATH}changeGraphDatabase.sh"
							sh "sudo service neo4j start"
							sh "sudo service tomcat7 start"
							sh "rm ${env.JENKINS_HOME_PATH}changeGraphDatabase.sh"
						}
					}
				}
			}			
		}
		// This stage runs the graph-qa script that will be emailed to Curation.
		stage('Post: Run graph-qa'){
			steps{
				script{
					cloneOrPullGitRepo("graph-qa")
					dir("graph-qa"){
						sh "mvn clean compile assembly:single"
						withCredentials([usernamePassword(credentialsId: 'neo4jUsernamePassword', passwordVariable: 'pass', usernameVariable: 'user')]){
							sh "java -jar target/graph-qa-jar-with-dependencies.jar -u $user -p  $pass --verbose"
						}
					}
				}
			}
		}
		// This stage emails the contents of the GraphQA_Summary_vXX.csv that is generated by graph-qa for review. 
		stage('Post: Email graph-qa output'){
			steps{
				script{
					emailext (
						body: "Hello,\n\nThis is an automated message from Jenkins regarding an update for v${currentRelease}. The Orthoinference step has finished running. Attached to this email should be the summary report output by graph-qa. Please compare this with the graph-qa output from the previous release and confirm if they look appropriate with the developer running Release. \n\nThanks!",
						to: '$DEFAULT_RECIPIENTS',
						from: "${env.JENKINS_RELEASE_EMAIL}",
						subject: "Orthoinference graph-qa for v${currentRelease}",
						attachmentsPattern: "**/graph-qa/reports/GraphQA_Summary_v${currentRelease}.csv"
						
					)
				}
			}
		}
		// All databases, logs, reports, and data files generated by this step are compressed before moving them to the Reactome S3 bucket. 
		// All remaining files/folders are then deleted that are not a part of the release-orthoinference repository.
		stage('Post: Archive Outputs'){
			steps{
				script{
					def s3Path = "${env.S3_RELEASE_DIRECTORY_URL}/${currentRelease}/orthoinference"
					sh "mkdir -p databases/ data/ reports/"
					sh "mv --backup=numbered *_${currentRelease}_*.dump.gz databases/"
					sh "mv graph-qa/logs/* logs/"
					sh "mv *.txt data/"
					// Keep this in orthoinference directory for symlink
					sh "mv data/report*sorted.txt ."
					sh "mv graph-qa/reports/* reports/"
					sh "gzip data/* logs/* reports/*"
					sh "aws s3 --no-progress --recursive cp databases/ $s3Path/databases/"
					sh "aws s3 --no-progress --recursive cp logs/ $s3Path/logs/"
					sh "aws s3 --no-progress --recursive cp data/ $s3Path/data/"
					sh "aws s3 --no-progress --recursive cp reports/ $s3Path/reports/"
					sh "rm -r databases logs data reports orthopairs"
					sh "rm -rf graph-importer*"
					sh "rm -rf graph-qa*"
					sh "rm -rf release-jenkins-utils*"
				}
			}
		}
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
// Utility function that checks if a git directory exists. If not, it is cloned.
def cloneOrPullGitRepo(String repoName) {
	// This method is deceptively named -- it can also check if a directory exists
	if(!fileExists(repoName)) {
		sh "git clone ${env.REACTOME_GITHUB_BASE_URL}/${repoName}"
	} else {
		sh "cd ${repoName}; git pull"
	}
}
