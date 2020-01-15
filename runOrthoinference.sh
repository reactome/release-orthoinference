#!/bin/bash
set -e

DIR=$(dirname "$(readlink -f "$0")") # Directory of the script -- allows the script to invoked from anywhere
cd $DIR

## Update repo
git pull
## Create new jar file with orthoinference code
mvn clean compile assembly:single

## Ensures the correct jar file is obtained regardless of orthoinference project version
orthoinference_jar_file=$(ls target/orthoinference-*-jar-with-dependencies.jar)

## Run orthoinference for each species
allSpecies=(mmus rnor cfam btau sscr drer xtro ggal dmel cele ddis spom scer pfal)
for species in "${allSpecies[@]}"
do
	echo "java -jar $orthoinference_jar_file $species > orthoinference_$species.out";
	java -jar $orthoinference_jar_file $species > orthoinference_$species.out;
done

echo "Orthoinference complete"