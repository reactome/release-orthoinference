#!/bin/bash
## This goes through the output report from orthoinference and organizes the lines into the order expected by Curation.
## author: jcook
set -e

## Parse command-line arguments.
releaseNumber=
while (( "$#" )); do
        case "$1" in
                -r|--release)
                        releaseNumber=$2;
                        shift 2;
                        ;;
                -*|--*=)
                        echo "Error: Unsupported flag $1"
                        exit 1
        esac
done

## If missing proper argument, explain usage.
if [ -z "$releaseNumber" ]
then
        echo "Output sorted orthoinference report file"; 
        echo "Usage: bash formatOrthoinferenceReport.sh --release releaseNumber ";
        exit 1
fi

## Create new sorted report file.
sortedReportName=report_ortho_inference_test_reactome_$releaseNumber_sorted.txt
rm -f $sortedReportName
touch $sortedReportName
## Add header
echo "#Inferred reactions report for v$releaseNumber" >> $sortedReportName

## Iterate through ordered species list, grep for species in original report file, and append the result to the new report file.
sortedSpecies=(ddis pfal spom scer cele sscr btau cfam mmus rnor ggal xtro drer dmel)
for species in "${sortedSpecies[@]}"
do
	grep "$species" report_ortho_inference_test_reactome_$releaseNumber.txt >> $sortedReportName
done

echo "Successfully generated $sortedReportName"
