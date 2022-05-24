#!/bin/bash

## 1) Iterates through all output text files (exported as-is from Proteome Discoverer) and gathers data into a single file.
## 2) Then, exports from that file only unique records.
## 3) Then, pulls only the following columns for downstream analysis: `Annotated Sequence`, `Modifications`, `Charge`, `Precursor Abundance`, `Spectrum File`, and `File ID`.
## 4) Finally, to the same file, adds a column recording modification status as yes (1) or no (0), done by pattern matching to the substring "GEE"


#### UPDATE: Column extraction algorithm in step 1 has been changed to seeking columns by their actual name, rather than their position in the file (which can vary).
#### This ensures the same data is always pulled regardless of file's contents. Ex) The column named "Precursor Abundance" will always have Precursor Abundance values
#### in it, whereas column #10, #12, or #36 may all have the desired Precursor Abundance values (depending on the file's contents). Contents between files are most 
#### likely to vary if one result was searched with the Precolator node while the other was not.

## function to be called at the very end of data analysis, when the final datafile needs peptide positions to correctly graph % label by peptide and map labels to structure
function 4_CALCULATE!() {

    echo -e $msg1"Now matching unique sequences and precursor abundances and summing total abundances associated with each labeling event... "

    len=$(awk -F"," 'NR > 1' $file3)

    # Use the number of unique peptide IDs as the searching index
    #for (( i = 1; i -le $len; i++ )); do
    for f in ${fnames[*]}; do

        echo "replicate name is "$f


        # This is the nested part of the loop. Here, a master peptide ID is stored in m. Then both replicate and m will be passed to awk as vars for pattern-matching and calculations.
        for m in ${mseqs[*]}; do
            echo -e $msg3"sequence is "$m
        
            ## NEXT: some awk block here that takes these bash variables (see below) and uses them to search trimmedFiltered and sum abundances.

        done
        # awk will then use those vars to pattern-match abundances with the master sequence. Remember that within awk the value of each PSM string must be uppercased
        # and then checked for a match with the current value of m (if ~ /*/ etc.)

        echo ""

        # 1. awk provides a function, "toupper", which capitalizes the string given to it.
        # $ awk '{print toupper($0)}' file1
        # HELLO
        # WELCOME
        # STRING

    done

    return

}

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

        ## in order of appearance, use each unique filename in fnames and use it to gather all values associated with it; then print to a temporary file
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

    return
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function 2_ExtractData() {

    # Use awk to export the desired columns from MergedPSMs to new text file "trimmed"
    echo -e $msg1"Extracting data... "

    awk -F "\t" '
    NR==1 {
        for (i=1; i<=NF; i++) {
            f[$i] = i
        }
    }
    { print $(f["Annotated Sequence"]) "\t" $(f["Modifications"]) "\t" $(f["Charge"]) "\t" $(f["Spectrum File"]) "\t" $(f["Precursor Abundance"]) "\t" $(f["XCorr"]) "\t" $(f["Rank"])}
    ' $file1 >> $file2

    echo -e $msg2"Copied 'Annotated Sequence' 'Modifications' 'Charge' 'Spectrum File' 'Precursor Abundance' 'XCorr' 'Rank' to "$file2

    # Repeat for to get sequence/position matches from merged peptide groups 
    awk -F "\t" '
    NR==1 {
        for (i=1; i<=NF; i++) {
            f[$i] = i
        }
    }
    { print $(f["Annotated Sequence"]) "\t" $(f["Positions in Master Proteins"]) "\t" $(f["Modifications"]) }
    ' $file1A >> $temp 

    # Now isolate unique sequence-position pairs which have been GEE labeled; store in file3A
    awk -F"\t" 'NR==1' $temp > $file2A; awk -F"\t" 'NR>=2' $temp | grep "GEE" | sort | uniq >> $file2A
    awk -F"\t" 'NR==1 {print $1"\t",$2}' $file2A > $file3A; awk -F"\t" 'NR > 1 {print $1"\t",$2}' $file2A | sort | uniq >> $file3A
    numUniqPepGrps=$(awk 'NR > 1' $file3A | wc -l); echo -e $msg2"Extracted "$numUniqPepGrps" unique, labeled sequence groups with their positions from among "$actual1A" total peptide groups."

    # Now store the unique IDs in a bash array for later use
    mseqs=($( awk -F"\t" 'NR > 1 {print $1}' $file3A))

    ## First print headers to File 3
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

    # count all original files to be processed
    numPSMfiles=$(ls PDoutputTextFiles/*PSMs* | wc -l); numPepGroupfiles=$(ls PDoutputTextFiles/*PeptideGroups* | wc -l)

    echo -e $msg1"Beginning file merge..."; echo -e $msg2"Identified "$numPSMfiles" PSM output files and "$numPepGroupfiles" Peptide Group output files to merge."

    # Get the unique headers from PSM files into the new merged PSM file, and repeat for Peptide Groups.
    awk -F"\t" 'NR==1' PDoutputTextFiles/*PSMs* > $file1; awk -F"\t" 'NR==1' PDoutputTextFiles/*PeptideGroups* > $file1A
    
    # Now merge the file contents into their respective mergefiles, and remove all ". Also, build an array fnames() containing unique filenames for later use
    for filename in PDoutputTextFiles/*PSMs*; do

        # sed's -i flag edits the file 'in-place' without requiring a temporary file
        awk -F"\t" 'NR > 1' $filename >> $file1; sed -i 's/\"//g' $file1

        # use _ as a delimiter to print only the core of the filename
        name=$(basename $filename | awk -F"_" '{print $2}' | awk -F"-" '{print $1}'); fnames+=( $name )

    done

    for filename in PDoutputTextFiles/*PeptideGroups*; do

        awk -F"\t" 'NR > 1' $filename >> $file1A; sed -i 's/\"//g' $file1A          # sed's -i flag edits the file 'in-place' without requiring a temporary file

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

    ## Create variables for merged-output files
    file1=mergedPSMs.txt; file1A=mergedPeptideGroups.txt

    ## Create  variables for data analysis snapshot files
    file2=trimmed.txt; file2A=SequencePositions.txt; file3=trimmedFiltered$ext; file3A=GEE_sequencePositions.txt; temp=temp$ext

    # user message flag
    msg1=$(echo "\n>> "); msg2=$(echo "\t~ "); msg3=$(echo "\t! ")

    ## Before starting, remove previous output files from the working directory
    if [ -e $file1 ]; then rm $file1; fi
    if [ -e $file1A ]; then rm $file1A; fi
    if [ -e $file2 ]; then rm $file2; fi
    if [ -e $file2A ]; then rm $file2A; fi
    if [ -e $file3 ]; then rm $file3; fi
    if [ -e $file3A ]; then rm $file3A; fi
    if [ -e $temp ]; then rm $temp; fi

    ## Initialize empty output files, descriptive variables and data structures
    touch $file1; touch $file1A; touch $file2; touch $file3            # output files to store an image of the data at each step of analysis
    declare -i numPSMfiles; numPSMfiles=0                              # var to keep track of data processing steps
    declare -i numPepGroupfiles; numPepGroupfiles=0                    # ditto
    declare -i numgoalPSMs; numgoalPSMs=0                              # ditto
    declare -i numgoalPepGrps; numgoalPepGrps=0                        # ditto
    declare -i actual; actual=0                                        # ditto
    declare -i actual1A; actual1A=0                                    # ditto
    declare -a fnames=()                                               # Array. Store the basenames of the original .raw files here
    declare -a -u mseqs=()                                             # Array. Store the unique master peptide sequences across all .raw files here
    
    return
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function assignPositions() {   # NOTE now that this block is a function, need to ensure that the variables it uses are visible to it (they may need to be passed to it!)
    
    ## 5/23/2022 THIS FUNCITON IS NO LONGER NEEDED HERE. DOING ASSIGNMENTS BY PRINTING UNIQUE IDS WITH THEIR CORRESPONDING POSITIONS IS ANALOGOUS TO CREATING A PERMANENT ARRAY
    ## OF THESE RECORDS, WHICH CAN BE EASILY ACCESSED BY STANDARD FILE READING. USE THIS FOR REFERENCE WITH OTHER CODE BLOCKS, BUT EVENTUALLY DEPRECATE THIS.

    # initialize the variables
    col1=""; refseq=""; col2=""; col3=""; col4=""; col5=""; col6=""; position=""

    # use the rowcount in trimmmed.csv to store the position value, row by row, and print it appended with the other row values to a new file (updated version of trimmed.csv)
    echo -e "Matching sequences with their positions in the master protein..."

    ## REPURPOSE THIS BLOCK FOR USE IN THE SEARCHING-PRECURSOR SUMMING CODE
    ## this can potentially be done in awk (more efficient).
    for (( i = 1; i -le $lenTrim; i++ )); do
        
        # store the current PSM sequence value
        col1=$(awk -F "," -v i=$i 'NR==i {print $1}' $file3)

        # Explanation of awk substr: Syntax is to print from the original input string ($0) starting at char position 4 (non-0 indexed, so starting at '.')
        # and continuing for the *original* string's length, minus 6 chars (ie, ignoring the flanking [char], 6 chars total)
        refseq=$(echo ${col1^^} | awk '{print substr($0, 4, length($0)-6)}')
        col2=$(awk -F "," -v i=$i 'NR==i {print $2}' $file3)
        col3=$(awk -F "," -v i=$i 'NR==i {print $3}' $file3)
        col4=$(awk -F "," -v i=$i 'NR==i {print $4}' $file3)
        col5=$(awk -F "," -v i=$i 'NR==i {print $5}' $file3)
        col6=$(awk -F "," -v i=$i 'NR==i {print $6}' $file3)

        ## if this is the first row, it is a header - store that value as such
        if [[ ! $i -gt 1 ]]; then
            position=$(echo "Position in Master Protein")

        ## repurpose this code for searching... the position match functionality is no longer needed but this has applications in Step 4
        ## if this is not the first row, it contains position values - extract the exact sequence-position match
        else 
            position=$(grep -w -m1 $refseq $file2A | cut -f2)   #-m1 forces grep to stop after the first match. Without it, near-exact matches are also made
                                                                # important: these are subsequences of the same peptide. It will be necessary to check for this situation
                                                                # and then collapse the ID values together as part of the matching process.
                                                                # Furthermore, this should be refined to retrieve *only* the integer characters so that they can be seamlessly
                                                                # used later in selection-expressions when running the pymol mapping script in a subshell to this script.

            echo "Sequence "$col1" was matched to positions "$position
        fi

        echo "$col1,$position,$col2,$col3,$col4,$col5,$col6" >> $temp
    
        ## only while testing this code block
        #if [[ $i -eq 50 ]]; then break; fi

    done; awk -F"," '{print $0}' $temp > $file3
    
    return
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Check for command line arguments
if [[ $# != 0 ]]; then 
    labelkey=$(echo ${1^^})     # accepts a string, case insensitive and content insensitive. labelkey will be used to search for the label specified by this string.
    graph=$(echo ${2^^})        # accepts true or false, case insensitive. Allows user to request the final data to be graphed for them (will call a separate script).        
    map=$(echo ${3^^})          # accepts true or false, case insensitive. Allows user to request the final data to be mapped to a protein structure in a pymol session (will call a separate script).
fi

## Begin (Step 0): Initialize variables and data structures
0_Initialize

## Step 1: Merge the output files by calling 'mergeFiles'
1_MergeFiles

## Step 2: Extract the desired information by calling ExtractData
2_ExtractData

## Step 3: Clean the data (remove redundant records and sort the unique records)
3_CleanData

## Step 4: Calculate percent labeling of each unique PSM in each file.
4_CALCULATE!


#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| 5/20/2022: Everything to this point has been tested and works correctly! |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######
#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| Now break all of this up into separate functions, to be called down here |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######
#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| This will organize the code and make development easier to work on |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######




### FINAL SEGMENT: DO ABUNDANCE MATCHING, SUMMATION AND % LABELING CALCULATIONS. USE *ONLY* UNIQUE, GEE-LABELED IDs AS THE SEARCH KEYS
######################################################################################################################################
## Extract relevant data from output file 3 into separate arrays for cross-referencing by array index. Printing by column number here is safe because
## the output file made by this script are always formatted the same way regardless of how they started.


## Put it all in one awk statement, similar to above (line 183). Think about it, why bother creating multiple arrays when you coud just awk the entire file
## (which itself is structured as an array anyway)
## DEPRECATE EVERYTHING FROM THIS POINT, AND MAKE IT INTO ONE BIG AWK STATEMENT INSTEAD!
## OPTIONAL: BUILD A SMALL TEST FILE STRUCTURED IN THE SAME WAY AS FILE3. Or, just use File 3 itself.
## step 1 ('BEGIN'): populate awk array with unique peptide IDs that were labeled with GEE.
## step 2: populate a second awk array with unique raw file names.
## step 3: awk loop - for element in raw array, iterate through the peptide ID array. For the current value of the peptide ID array, search the trimmedFiltered.csv .
## for IDs which match the current element of both arrays. When a GEE match is found, add its value to a variable storing GEE abundance. Else, store add into a non-GEE variable.
## A third variable, %label, is equal to GEE/nonGEE+GEE * 100 and updates in real time.
## step 4: Continue until the given ID value has no more matches. At the bottom of this loop, print the current raw file, the ID value, the ID's summed GEE abundance,
## its summed nonGEE abundance, and the ID's %labeled value.
## step 5: Repeat step 4 until all ID values have been iterated and printed with their abundances.
## step 6: Repeat the nested loops until all raw files have been sampled for the unique IDs and abundances and everything has been printed. This is the calculated abundances file.
## step 7: New awk block: reading the calculated abundances file, use pattern matching with the same awk array of peptide IDs as above to iterate through the abundances file and get theh average.
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