#!/bin/bash

##### TODO!!!!!!! (12/05/2022): 1) Update code to correctly flag undetected peptides as 'undetected' instead of defaulting their calculated %change to 0%
##### TODO!!!!!!! (12/05/2022): 2) Add additional code so that the final output file gets re-printed to a new file, with replicates correctly positioned in the same row as each other
##### TODO!!!!!!! (12/05/2022): ...in the form of "Sequence","Position","%change1","%change2","%change3","%change4" and so on.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function 4_CALCULATE!() {

    temp1=tempByFile.csv; temp2=tempByFile_Peptide.csv; temp3=temp3.csv

    echo -e $msg1"Now matching unique sequences and precursor abundances and summing total abundances associated with each labeling event... "
    
    # store an all-uppercased version of file3 in temp; simplifies peptide ID matching further down
    # i is for loop breaks during implementation and testing (see further down)
    i=1

    # make a temporary copy of trimmedFiltered with all chars uppercased (affects filenames and peptide sequences)
    awk -F"," '{print toupper($0)}' $file3 > $temp

    # Begin searching for mods in each file, one file at a time
    for f in ${fnames[*]}; do

        echo "Replicate Filename","Position in Master Protein","Master Peptide Sequence","Labeled","Unlabeled","% labeled" >> $file4
        fupper=$(echo ${f^^})   # stores uppercased version of filename "f". Necessary because this is used to search an uppercased image of file3 (see codeline #15)

        ## Handle the data in subsets by exporting to temporary files and reading from those files (simplifies search process). j-iterator assists with seq-position matching.
        j=1
        for m in ${mseqs[*]}; do

            # Sort first by file ID, and use awk substr() to strip away flanking "[X]." and ".[Y]" characters on the peptide sequence string
            awk -F"," -v f=$fupper '$0 ~ f{
                AnnotatedSequence = substr($1, 5, length($1)-8); Modifications = $2; Charge = $3; SpectFile = $4; PrecAbundance = $5; Label = $6};
                {print AnnotatedSequence "," Modifications "," Charge "," SpectFile "," PrecAbundance "," Label}' $temp > $temp1

            # Extract data associated only with current value of "m" from within the file-sorted data
            awk -F"," -v m=$m '$1 == m{ print $0 }' $temp1 > $temp2

            # Log file for recording PSM pattern-matching behavior (for quality control/debugging purposes)
            echo -e "__________________________________________________________________\nCurrent file f = "$fupper", current sequence m = |"$m"|\nMatched PSMs shown below:\n" >> $log
            cat $temp2 >> $log
            echo -e "__________________________________________________________________\n\n" >> $log 

            ## Extract the position of the current peptide. (j iterates until each sequence 'm' is matched to 'position', then resets in the outer loop for the next file 'f'.
            position=$(awk -F"," -v j=$j 'NR==j {print $2}' Sequence_Position.csv); echo -e "\t\t\tCurrent file = "$fupper", position = "$position" with sequence = "$m

            ## Sum each abundance type separately, then calulate percent abundance of labeled peptides.
            Unlabeled=$(awk -F"," '$0 !~ /GEE/{ sum += $5 } END { print sum }' $temp2)    # changed from temp2 to temp1 during testing
            Labeled=$(awk -F"," '$0 ~ /GEE/{ sum += $5 } END { print sum }' $temp2)       # changed from temp2 to temp1 during testing
            Perclabeled=$(awk -F"," -v Labeled=$Labeled -v Unlabeled=$Unlabeled 'BEGIN{ sum = Labeled + Unlabeled;
                            if ( sum > 0 )
                                perc = 100*Labeled/sum
                            else
                                perc = 0 }END{ print perc }' $file4)

            ## Finally, print the labeling data to final output file.
            echo $f,$position,$m,$Labeled,$Unlabeled,$Perclabeled >> $file4   ### SUCCESS! The calculations are accurate (match manual). Done in ~1 min or less!
                                                                    ### NOTE: May be possible to make this leaner by compressing this interior loop's body
                                                                    ### into a single awk statement and eliminate temporary file handling. Worry about that later.
            let "j++"          # used for awking out the corresponding position
        
        done

        ## a test condition to limit the calculations during debugging. Note that the string comparison here is not working because passing "test" does not trigger a break.
        #echo "labelkey == "$labelkey
        if [[ $i -eq 5 && $# != 0 && $labelkey == "TEST" ]]; then
            echo -e $msg3"Testing has triggered a -break- statement to cancel remaining calculations.\n"; break
        fi

        ## add a newline between each block of data for easier reading
        echo "" >> $file4

        let "i++"              # currently unused

    done
    echo -e $msg2"Done. Results of labeling calculations are stored in '"$file4"'."
    rm temp*    # delete all temporary files at the end to keep working directory clean

    return

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function 3_CleanData() {

    ## Now remove redundant records where applicable (uses a temporary file); NOTE functions as expected (compared to manual work in Excel)
    echo -e $msg1"Now checking '"$file3"' for redundant PSMs and sorting the unique records by associated spectrum file... "
    touch $temp
    awk -F "," 'NR == 1 {print $1",", $2",", $3",", $4",", $5",", $6}' $file3 > $temp
    awk -F "," 'NR >= 2 {print $1",", $2",", $3",", $4",", $5",", $6}' $file3 | sort -k4n | uniq >> $temp
    cat $temp > $file3

    ## Track the data
    trim=$(awk 'NR >= 2' $file3 | wc -l); diff=$((lenfile-trim))

    ## User message
    if [ $trim != 0 ]; then
        echo -e $msg2"Found and removed "$diff" duplicate records ("$trim" records retained in $file3)."
    else
        echo -e $msg3"Found 0 duplicate records; all records retained."
    fi

    ## Use the var 'numPSMfiles' to re-sort PSM records by source filename (which corresponds to original number of input files)
    awk -F"," 'NR == 1' $file3 > $temp
    for ((i = 0; i < $numPSMfiles; i++)); do

        ## in order of appearance, use each unique filename in fnames to gather all data associated with it into a temporary file
        identifier=$(echo ${fnames[i]}); grep $identifier $file3 >> $temp

    done

    # bounce the sorted records back to File 3 and reset the temp file
    awk -F"," '{print $0}' $temp > $file3; rm $temp

    ## Track the data
    trim1=$(awk 'NR >= 2' $file3 | wc -l)
    if [ $trim1 == $trim ]; then
        echo -e $msg2"Records have been sorted by source file (all records retained in $file3)."
    else
        echo -e $msg3"Error during record sorting! Review contents of "$file3" to identify what changed."
    fi

    ## use file3 to populate the mseqs array with unique GEE-labeled peptide sequences (update 6/01/2022: this works as desired (matches manual analysis in Excel))
    mseqs=($( awk -F"." 'NR > 1 && $0 ~ /GEE/{ print toupper($2) }' $file3 | sort | uniq ))
    numseqs=$(echo ${#mseqs[@]}); echo -e $msg2$numseqs" unique peptide sequences identified in "$file3

    ## Generate a reference file containing master sequence / position relationships
    if [ -e Sequence_Position.csv ]; then rm Sequence_Position.csv; fi
    echo -e $msg2"Writing sequence-position list..."
    mpos=""          # stores the current position in the master  sequence                               
    
    # Use the array mseqs to cross reference file3A (GEE_sequencePositions.txt): find the current sequence ("m") in the list file and return the corresponding position values
    for m in ${mseqs[@]}; do
        mpos=$(grep -w -m1 $m $file3A | cut -f2); mpos=$(echo $mpos | awk -F"[" '{print $2}' | awk -F"]" '{print $1}')
        echo $m,$mpos >> Sequence_Position.csv
    done

    # Finally, sort the unique reference peptide IDs by their position in the protein. This file will be used for sequence-matching in the final file!
    sort -t"," -k 2n -o Sequence_Position.csv{,}
    echo -e $msg2"Sequence-position list file completed."

    unset mseqs; mseqs=($( awk -F"," '{ print $1 }' Sequence_Position.csv ))

    return
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function 2_ExtractData() {

    # Use awk to export the desired columns from MergedPSMs to new text file "trimmed"
    echo -e $msg1"Extracting data... "

    # Awk loop: this builds an associative awk array, where each element is an entire record. The columns are then printed by the pattern matched in each element.
    # because of how this works, the RAW output files need to be named simplistically such as *_MB1_* and so on; otherwise, add code to assign arbitrary file IDs at
    # the very beginning to avoid this problem.
    awk -F "\t" '
    NR==1 {
        for (i=1; i<=NF; i++) {
            f[$i] = i
        }
    }
    { print $(f["Annotated Sequence"]) "\t" $(f["Modifications"]) "\t" $(f["Charge"]) "\t" $(f["Spectrum File"]) "\t" $(f["Precursor Abundance"]) "\t" $(f["XCorr"]) "\t" $(f["Rank"])}
    ' $file1 >> $file2

    echo -e $msg2"Copied 'Annotated Sequence' 'Modifications' 'Charge' 'Spectrum File' 'Precursor Abundance' 'XCorr' 'Rank' to "$file2

    # Repeat for to get sequence/position matches from merged peptide groups.
    awk -F "\t" '
    NR==1 {
        for (i=1; i<=NF; i++) {
            f[$i] = i
        }
    }
    { print $(f["Annotated Sequence"]) "\t" $(f["Positions in Master Proteins"]) "\t" $(f["Modifications"]) }
    ' $file1A >> $temp 

    # Now isolate UNIQUE sequence-position pairs which have been GEE labeled; store in file3A
    awk -F"\t" 'NR==1' $temp > $file2A; awk -F"\t" 'NR>=2' $temp | grep "GEE" | sort | uniq >> $file2A
    awk -F"\t" 'NR==1 {print $1"\t",$2}' $file2A > $file3A; awk -F"\t" 'NR > 1 {print $1"\t",$2}' $file2A | sort | uniq >> $file3A
    numUniqPepGrps=$(awk 'NR > 1' $file3A | wc -l); echo -e $msg2"Extracted "$numUniqPepGrps" unique, labeled sequence groups with their positions from among "$actual1A" total peptide groups."

    # First print headers to File 3
    echo -e $msg1"Filtering out low-scoring PSMs and annotating labeled PSMs..."
    Seq=$(awk -F "\t" 'NR==1 {print $1}' $file2)
    Mods=$(awk -F "\t" 'NR==1 {print $2}' $file2)
    Z=$(awk -F "\t" 'NR==1 {print $3}' $file2)
    name=$(awk -F "\t" 'NR==1 {print $4}' $file2) 
    Precursor=$(awk -F "\t" 'NR==1 {print $5}' $file2)
    LabelStatus="Label?"
    Score=$(awk -F "\t" 'NR==1 {print $6}' $file2)
    Rank=$(awk -F "\t" 'NR==1 {print $7}' $file2); echo -e $Seq,$Mods,$Z,$name,$Precursor,$LabelStatus,$Score,$Rank > $file3

    ### Find records meeting required Xcorr and Rank scores, and tag labeled PSMs; print to File 3
    awk -F "\t" 'BEGIN{ seq = ""; mods = ""; charge = 0; name = ""; abundance = 0; label = ""; xcorr = 0; rank = 0 };

            NR >= 2 { seq = $1; mods = $2; charge = $3; name = $4; abundance = $5; xcorr = $6; rank = $7; 

            if ( mods ~ /GEE/ )
                label = "GEE";

            else
                label = "";

            if ( xcorr >= 2 && rank == 1)
                print seq "," mods "," charge "," name "," abundance "," label "," xcorr "," rank}' $file2 >> $file3

    ## Report the number of removed records for failing XCorr threshold   ### NOTE functions as expected (compared to manual work in Excel)
    lenfile=$(awk 'NR >= 2' $file3 | wc -l)
    filtered=$((actual-lenfile))

    ## Copy to backup file
    head -n1 $file3 > backup_$file3
    awk -F"," 'NR >= 2' $file3 | sort -k4n >> backup_$file3

    ## User message
    if [ $actual != 0 ]; then
        echo -e $msg2"Found and removed "$filtered" PSMs with insufficient Xcorr score and/or rank ("$lenfile" records retained)."
        echo -e $msg2"These changes are also copied in 'backup_"$file3"', review as needed." 
    else
        echo -e $msg3"Found 0 failing PSMs; all records retained."
    fi

    return
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function 1_MergeFiles() {

    echo -e $msg1"Beginning file merge..."; echo -e $msg2"Identified "$numPSMfiles" PSM output files and "$numPepGroupfiles" Peptide Group output files to merge."

    # Get the unique headers from PSM files into the new merged PSM file, and repeat for Peptide Groups.
    awk -F"\t" 'NR==1' PDoutputTextFiles/*PSMs* > $file1; awk -F"\t" 'NR==1' PDoutputTextFiles/*PeptideGroups* > $file1A
    
    # Now merge the file contents into their respective mergefiles, and remove all "".
    for filename in PDoutputTextFiles/*PSMs*; do

        # sed's -i flag edits the file 'in-place' without requiring a temporary file. This line eliminates "" from the file.
        awk -F"\t" 'NR > 1' $filename >> $file1; sed -i 's/\"//g' $file1

    done
    for filename in PDoutputTextFiles/*PeptideGroups*; do
        awk -F"\t" 'NR > 1' $filename >> $file1A; sed -i 's/\"//g' $file1A          #NOTE sed's -i flag edits the file 'in-place' without requiring a temporary file
    done

    # count all non-header rows that should have been exported
    numgoalPSMs=$(awk -F"\t" 'NR > 1' PDoutputTextFiles/*PSMs* | wc -l); numgoalPepGrps=$(awk -F"\t" 'NR > 1' PDoutputTextFiles/*PeptideGroups* | wc -l)

    # count actual export numbers
    actual=$(awk 'NR > 1' $file1 | wc -l); actual1A=$(awk 'NR > 1' $file1A | wc -l)

    ## check
    echo -e $msg2"Found "$numgoalPSMs" PSM records and "$numgoalPepGrps" Peptide Groups records to export; exported "$actual" PSM records and "$actual1A" Peptide Groups records (excluded redundant header records)."
    echo -e $msg2"Merged data were saved in '"$file1"' and in '"$file1A"'."  

    return
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function 0_Initialize() {
    # Begin by storing output filenames in variables for ease of use
    ext=.csv

    ## Create variables for storing analysis snapshots in outfiles
    file1=mergedPSMs.txt; file1A=mergedPeptideGroups.txt
    file2=trimmed.txt; file2A=SequencePositions.txt
    file3=trimmedFiltered$ext; file3A=GEE_sequencePositions.txt
    file4=final$ext; temp=temp$ext; log=AnalysisLog.txt

    # user message flag
    msg1=$(echo "\n>> "); msg2=$(echo "\t~ "); msg3=$(echo "\t! ")

    ## Before starting, remove previous output files from the working directory
    if [ -e $file1 ]; then rm $file1; fi
    if [ -e $file1A ]; then rm $file1A; fi
    if [ -e $file2 ]; then rm $file2; fi
    if [ -e $file2A ]; then rm $file2A; fi
    if [ -e $file3 ]; then rm $file3; fi
    if [ -e $file3A ]; then rm $file3A; fi
    if [ -e $file4 ]; then rm $file4; fi
    if [ -e $temp ]; then rm $temp; fi
    if [ -e $log ]; then rm $log; fi

    ## Initialize empty output files, descriptive variables and data structures
    touch $file1; touch $file1A; touch $file2; touch $file3; touch $file4; touch $log   # output files to store an image of the data at each step of analysis
    declare -i numPSMfiles; numPSMfiles=0                                               # var to keep track of data processing steps
    declare -i numPepGroupfiles; numPepGroupfiles=0                                     # ditto
    declare -i numgoalPSMs; numgoalPSMs=0                                               # ditto
    declare -i numgoalPepGrps; numgoalPepGrps=0                                         # ditto
    declare -i actual; actual=0                                                         # ditto
    declare -i actual1A; actual1A=0                                                     # ditto
    declare -a fnames=()                                                                # Array. Store the basenames of the original .raw files here
    declare -a -u mseqs=()                                                              # Array. Store the unique (uppercased) master peptide sequences across all .raw files here

    return
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## BEGIN ANALYSIS HERE
## Check for command line arguments
if [[ $# != 0 ]]; then 
    labelkey=$(echo ${1^^})     # accepts a string, case insensitive and content insensitive. labelkey will be used to search for the label specified by this string.
    graph=$(echo ${2^^})        # accepts true or false, case insensitive. Allows user to request the final data to be graphed for them (will call a separate python script).        
    map=$(echo ${3^^})          # accepts true or false, case insensitive. Allows user to request the final data to be mapped to a protein structure in a pymol session (will call the same python script).
fi

## Begin (Step 0): Initialize variables and data structures; also tally the number of input files present
0_Initialize

# Calculate the number of files to process
numPSMfiles=$(ls PDoutputTextFiles/*PSMs* | wc -l); numPepGroupfiles=$(ls PDoutputTextFiles/*PeptideGroups* | wc -l)

# Populate fnames array in advance, so that successive methods can access its list of names; fnames stores the .raw source filenames associated with the data
for filename in PDoutputTextFiles/*PSMs*; do
    name=$(basename $filename | sed 's/_PSMs.txt//' | awk -F"-" '{print $1}'); name=$(echo $name".raw"); fnames+=( $name )
done

# Display the identified filenames to user
fnamestring=""; echo -e "\nFound the following filename abbreviations: "
for i in ${fnames[*]}; do
    fnamestring=$(echo -e $fnamestring$i" | ")
done
echo -e $fnamestring"\n"

## Step 1: Merge the output files by calling 'mergeFiles'
1_MergeFiles

## Step 2: Extract the desired information by calling ExtractData
2_ExtractData

## Step 3: Clean the data (remove redundant records and sort the unique records)
3_CleanData

## Step 4: Calculate percent labeling of each unique PSM in each file.
4_CALCULATE!


#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| 11/09/2022: Everything to this point has been tested and works correctly! |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######
#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| Script is now complete; accurately calculates abundances and %labeling |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######
#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| NEXT: Debug for ambiguous cleavage in peptide IDs; re-tool to extract |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######
#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| unique peptide sequences from the PSMs in output file3 rather than from peptide groups |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######

### NOTES ###
######################################################################################################################################
## Extract relevant data from output file 3 into separate arrays for cross-referencing by array index. Printing by column number here is safe because
## the output file made by this script are always formatted the same way regardless of how they started.

## To get the average of the %labeled across files, use a counter to track the number of replicates. Sum all matched percent values, then divide at the end by the counter value.
## standard deviation on the other hand will be tricky. Perhaps instead of doing the statistics (including averaging) here, it might be time to call a script using matplotlib or R for this work.
## Because (well done!) it would be the final step!


## NOTE: Excel's 'remove duplicates' function is case insensitive. Therefore for PSMs labeled in different sites, this destroys data.
## NOTE: An approach that preserves this information while collapsing PSM sequences into a master sequence would be to create 2 arrays:
## NOTE: The first stores the unique instances of each string (by bash's uniq() function); this is the expanded list of sequences accounting for
## NOTE: multiple label sites in the same peptide sequence. The second stores only the unique peptide "master" (case-insensitive) sequence - ie.,
## NOTE: the list of sequences converted to uppercase and piped through sort | uniq. The second array also will have additionl elements, containing corresponding
## NOTE: values for unlabeled, labeled and percent all initialized as zero (0). Ie., alternating like this: masterList=(seq1, unlabeled1, labeled1, percent1, seq2, unlabeled2, labeled2, percent2 ... seq_n, unlabeled_n, labeled_n).
## NOTE: Bash does not support above 1D arrays but there is no limit to array size. So, use index with parameter expansion to access and update elements. This will require
## NOTE: retrieving the index of the sequence match, then adding +1 (for unlabeled's index) or +2 (for labeled's index) to get the correct index of the element to sum with new abundance values.
## NOTE: A loop with an awk function will iterate through the first array; awk will pattern match and do summation of peptide abundances (updating
## NOTE: the corresponding variables unlabeled, labeled and percent for this pass). A simple if-matches /*GEE*/ -else statement will control which
## NOTE: variable (labeled or unlabeled) gets updated. Then, the current value of $seq will be converted to uppercase,
## NOTE: matched to its corresponding master sequence in the second array, and then the second array's values will be summed with current values of 
## NOTE: labeled and unlabeled forthe current master sequence. During this looping, need to enforce pattern-matching to File ID value to keep replicates
## NOTE: distinct until final averaging and plotting.



## 11/09/2022 - extensions below still pending
## 6/02/2022 - HEADS UP - several major extensions to this need to be made to this code_________________________________________________________________________________________
## 1) This script is originally written to handle multiple result files from a "batched" search. Unbatched compresses all results into a single file.
##    Make sure that this script knows how to distinguish between the two formats and that it no longer tries to merge files by default. A simple if branch
##    controlling Function 1 should do the trick. Something like "if [[ numPSMfiles -eq 1 ]]; then pass; else; ...run existing code to build fnames array and merge files... "
##
## 2) Add methods for handling hydroxyradical footprinting. Think carefully about how to do this one - may be best to add functionality that lets the user submit a list file
##    as a command-line argument (like a csv, with a list of strings which are the modification keywords). Or, the user may specify a number of mods they want searched, and then
##    manually enter them one by one. Think about it.
##
## 3) This one should actually be completed first: appending positions in master proteins to the final file while simultaneously correcting amino acid index by the length of any
##    added tags. Use the unique sequence ID-position pairs in file3A (unique peptide groups) to do the pattern matching here.
##