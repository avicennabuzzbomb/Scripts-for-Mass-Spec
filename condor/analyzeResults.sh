#!/bin/bash

###### checks summary csv for patterns that imply a botched job (botched jobs occur
# because of the difficulty in predicting when a cif file has unusual features or a fasta format
# which the code in SASAquatch.py does not anticipate ###########

# checks for high incidence of buried residues (greater than 20%), for invalid fasta characters, and for unusual residue positions (integers < 1)
# then writes the PDB ID to a new list, which I can check manually in a pymol session)
# also gives a bulk job success rate at the top of the file. This script's output is information for troubleshooting edge cases that SASAquatch.py encounters.

## Get datestamp for the summary, alert user the analysis is starting, and begin writing the summary file
d=$(date +%Y-%m-%d)                                                                                                                                     ## LINE 10																				
echo "Scanning output files for signs of edge cases (failure)..."									
echo "Job assessment for job(s) completed after "$d > eval_$d.csv						
echo "Troubleshooting information compiled by individual PDB ID" >> eval_$d.csv										
echo "" >> eval_$d.csv																			
echo "Job Name,Reason for Followup" >> eval_$d.csv																
	
## Create a temporary file to store the jobs and their results as they are run
touch temp.csv
					
## Record the total number of job outputs being analyzed																																
filecount=$(ls log | wc -l); echo "Total number of jobs =" $filecount																	
																			## LINE 20

## declare and initialize threshold ($threshold) value for declaring an error (initialize as int type)
declare -i threshold=10
declare -i threshold2=35           # After verifying flagged ID's manually in pymol, I've seen missing occupancy as high as 33% of residues in a valid structure.
				   # Therefore 35(%) is the new `threshold2` value to check against, to allow those valid files to pass.

for filename in Results/*.csv; do

	# keep track of which file we're looking at currently (retrieves * from SASA_*_ALL, where * is only the PDB query ID; simplifies analysis)
	name=$(basename -s ".csv" " $filename" | sed 's/SASA_//g' | sed 's/_ALL//g')

	## vomit the csv contents into a pipe, and within the pipe count the number of rows under the header matching an unwanted pattern		## LINE 30
	let exposed=$(grep "Exposed" $filename | wc -l)
	
	## detect proportion of absent density ("Not present in structural model", or simply "N/A") in each file
	let NA=$(grep "N/A" $filename | wc -l)

	# get total number rows
	let totRows=$(cat $filename | wc -l)   # first, count all the rows in the file
	let fasta=$(grep "FASTA" $filename | wc -l)    # count up the subheaders in the file ("FASTA" only appears in subheader rows)
																			## LINE 40		
        # exclude instances of "Exposed" that occur in subheaders from residue rows									
        let exposed=$exposed-$fasta															
																			
	let fasta="fasta * 2 + 1"          # combine subheader count with whitespace and the main header						
	let total_rows=$totRows-$fasta     # subtract the header and whitespace rows from the residue rows						
	let bund="total_rows - exposed"	   # number of buried or undetectable residues in the model output																	
				
	## 1ST CHECK: is the current result file populated with data, or is it empty?									
	if [ $total_rows -gt 0 ]; then
							
		# calculate percentage of detected residues (bash is clunky at this, so it has to be this way:)							## LINE 50
		perc=$(echo "$exposed / $total_rows" | bc -l)													
		perc=$(echo "$perc * 100" | bc -l | xargs printf "%.*f\n" "$4")   # last section of pipe round val of $perc to <= 2 sig figs 

		# calculate the percentage of detected residues that have no density (instances of "N/A")
		percNA=$(echo "$NA / $total_rows" | bc -l)
		percNA=$(echo "$percNA * 100" | bc -l | xargs printf "%.*f\n" "$4")
	
		## 2ND CHECK: does the current result file contain too many residues flagged "Buried"? 
		if [ $perc -le $threshold ]; then
			# if 10% or less of the protein is rated as exposed, record that PDB ID for followup
			echo $name", (2) SASA: residues rated exposed vs. total residues: "$exposed"/"$total_rows "("$perc"% of the protein)" >> temp.csv     
		
		## 3RD CHECK: does the file describe an unusual amount of missing electron density?
		elif [ $percNA -gt $threshold2 ]; then																   		
			# if more than 1/3 of residues in the structure flag as no density detected, record the PDB ID for followup
			echo $name", (3) Occupancy: residues lacking modeled density in structure: "$NA"/"$total_rows "("$percNA"% of the protein)" >> temp.csv
		fi
	else 
		echo $name", (1) Null: gave an empty output file; possibly due to a density/occupancy stat error [check if this is an EM structure]" >> temp.csv
	fi
done

## 4TH CHECK: record the PDB IDs which needed too much memory to finish on the queue (they are left behind after *csv collection)
cd /
for filename in software/mrblackburn/pymol/*csv; do
	name=$(basename -s ".csv" " $filename" | sed 's/SASA_//g' | sed 's/_ALL//g')
	echo "$name, (4) Memory: job failed to complete on queue [probable job size could be too large]" >> home/mrblackburn/temp.csv
done

cd home/mrblackburn

## Finally, sort the flagged filenames in the temporary file by their error description (in field 2), and output to the eval file. Also add summary statistics to the bottom.
sort -k2 temp.csv >> eval_$d.csv

linecount=$(cat eval_$d.csv | wc -l)	# count entries in output file
let linecount="linecount - 4"		# ignore headers and whitespace (number of PDB names should be less than the number of Results files if failed jobs occurred)

echo ""; echo "Original number of queries: "$filecount; echo "Number of failed PDB queries: "$linecount
echo "" >> eval_$d.csv; echo "Original number of queries: "$filecount >> eval_$d.csv; echo "Number of failed PDB queries: "$linecount >> eval_$d.csv

est_success=$(echo "100 - $linecount / $filecount * 100" | bc -l | xargs printf "%.*f\n" "$4"); echo "Est. success rate is "$est_success"%"
echo "Estimated success rate is "$est_success"%" >> eval_$d.csv; echo ""; echo "Evaluation file has been stored in the directory 'evaluations' "

## Directory cleanup at the end
mv *eval_* evaluations
rm temp.csv
