#!/bin/bash

# runs a set of commands I've been repeatedly using to monitor jobs and outputs


# initialize variables
dir=""
declare -i empty=0
declare -i content=0
declare -i subset=0
declare -i sub=0

###################################
function checkstuff {

echo "_____________________________________________________________________"
echo "Contents of /"$1":"; content=$(ls $1 | wc -l)
if [ $content -gt $empty ]; then
	if [ $content -gt 10 ]; then
        	subset=10
	else
		subset=$(echo $content)
	fi	
	echo "Displaying $subset of $content detected file(s):"; ls -U $1 | head -n 10; echo "..."
else
	echo $1 "is currently empty."
fi
echo ""

}

####################################

##~~DRIVER CODE~~##

## check each output directory of interest; directory's name is used as the parameter during the call to `checkstuff{}`
checkstuff log; checkstuff output; checkstuff Results

## for keeping my own directory tidy, also relocatee all submit (*.sub) files to their own directory
sub=$(find . -maxdepth 1 -type f -name "*.sub" | wc -l)
if [ $sub -gt 0 ]
then
	echo "Found $sub submit files in this directory; they have been moved into /Submit_Files"; mv *.sub Submit_Files
fi

## check software/mrblackburn/pymol for csv files
cd /; echo "Output files (*.csv) in software/mrblackburn/pymol: "; cd software/mrblackburn/pymol; dir=$PWD; declare csv=$(find . -maxdepth 1 -type f -name "*.csv" | wc -l)

if [ $csv -gt $empty ]
then
	if [ $csv -gt 10 ]
	then
		echo "Displaying up to 10 of $csv detected file(s)"; find . -maxdepth 1 -type f -name "*.csv" | head -n 10; find . -maxdepth 1 -type f -name "*.cif" | wc -l
		echo "..."
	else
		echo "Displaying $csv detected file(s):"; find . -maxdepth 1 -type f -name "*.csv"; find . -maxdepth 1 -type f -name "*.cif" | wc -l
		echo "..."
	fi
else
	echo $dir "is currently empty."
fi

echo ""; echo "Checking job queue:"

condor_q
