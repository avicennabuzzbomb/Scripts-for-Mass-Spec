#!/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function 4_CALCULATE!() {

    temp1=tempByFile.csv; temp2=tempByFile_Peptide.csv; temp3=temp3.csv
    # temp3.csv is for re-sorting the PSMs by their sequence positions. At each pass of this loop (for f in...), print all results (except header!) instead to temp3; then
    # use awk -F"," '{print $0}' $temp3 | sort >> $file4                                                                        

    # for loop breaks during implementation and testing (see further down)
    i=0; j=0

    echo -e $msg1"Now matching unique sequences and precursor abundances and summing total abundances associated with each labeling event... "
    
    # store an all-uppercased version of file3 in temp; simplifies peptide ID matching further down
    awk -F"," '{print toupper($0)}' $file3 > $temp

    for f in ${fnames[*]}; do

        echo "Replicate Filename","Master Peptide Sequence","Labeled","Unlabeled","% labeled" >> $file4
        fupper=$(echo ${f^^})   # stores uppercased version of filename "f"

        # Here, a master peptide ID is referenced in m. Then both replicate and m will be passed to awk as vars for pattern-matching and calculations.
        for m in ${mseqs[*]}; do

            # Slice the data into separate streams by file and then by master peptide for awk to handle one a time
            # (eliminates unnecessary pattern matching steps during the search-sum process). Pattern matching is implicitly done by using temporary files to carve up the data
            # into subcategories, and then summing all abundance of each category. Hacky but it will have to do until I can determine why pattern matching is not working here.
            awk -F"," -v f=$fupper '$0 ~ f{ print $0 }' $temp > $temp1                  # filter data by current value of filename (f); WORKS
            awk -F"," -v m=$m '$0 ~ m{ print $0 }' $temp1 > $temp2                      # filter data again by current value of master sequence (m); WORKS

            ## Extract the position of the current peptide
            position=$(grep -w -m1 $m $file3A | cut -f2)
            #echo "Current sequence is "$m "at positions "$position

            ## Sum the abundances separated
            Unlabeled=$(awk -F"," '$0 !~ /GEE/{ sum += $5 } END { print sum }' $temp2)   #; echo "Unlabeled abundance is "$Unlabeled
            Labeled=$(awk -F"," '$0 ~ /GEE/{ sum += $5 } END { print sum }' $temp2)      #; echo -e "Labeled abundance is "$Labeled"\n_____________________"

            # this code block is currently incorrect; as is it expects a specific input file. Deal with this last. 8/03/2022 is this still true?
            Perclabeled=$(awk -F"," -v Labeled=$Labeled -v Unlabeled=$Unlabeled 'BEGIN{ sum = Labeled + Unlabeled;
                            if ( sum > 0 )
                                perc = 100*Labeled/sum
                            else
                                perc = 0 }END{ print perc }' $file4)

            ## finally, print to file
            echo $f,$m,$position,$Labeled,$Unlabeled,$Perclabeled >> $file4   ### SUCCESS! The calculations are accurate (match manual). Done in ~1 min or less!
                                                                    ### NOTE: May be possible to make this leaner by compressing this interior loop's body
                                                                    ### into a single awk statement and eliminate temporary file handling. Worry about that later.
            let "j++"          # currently unused
        
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

    ## use file3 to populate the mseqs array with unique GEE-labeled peptide sequences (update 6/01/2022: this works as desired (matches manual analysis in Excel))
    mseqs=($( awk -F"." 'NR > 1 && $0 ~ /GEE/{ print toupper($2) }' $file3 | sort | uniq ))
    numseqs=$(echo ${#mseqs[@]}); echo -e $msg2$numseqs" unique peptide sequences identified in "$file3


    ## From 8/20/2022:
    ## TODO: Do the positional matching with the unique sequences HERE; either populate the array with additional, corresponding elements
    ## TODO: (ie., values for "Position in Master Protein" from the Peptide Groups file), or just redirect everything to create a new reference temporary folder.
    ## TODO: This will create the opportunity to pre-sort unique sequences by position values, and also eliminates redundantly position-matching during the calculation
    ## TODO: phase (see function 4_CALCULATE!() )
    ### Currently testing. Sending unique, unsorted peptide sequence IDs to a separate temp file, then appending their positions from Peptide Groups and re-sorting by
    ### position in master protein sequence. Doing it here is 2-birds-1-stone

    if [ -e Sequence_Position.csv ]; then rm Sequence_Position.csv; fi

    echo -e $msg2"Writing sequence-position list..."

    mpos=""              # position in master
    for m in ${mseqs[@]}; do

        ## Extract the position of the current peptide
        mpos=$(grep -w -m1 $m $file3A | cut -f2); mpos=$(echo $mpos | awk -F"[" '{print $2}' | awk -F"]" '{print $1}')
        echo $m,$mpos >> Sequence_Position.csv

        # Before sorting, consider trying  awk -F"," 'NR>1{print $2}' posref.csv | awk -F"\[" '{print $2}' | awk -F"\]" '{print $1}'
        # to separate out the numeric values from the position string. This actually works for isolating the integer range from the string, confirmed by testing.
        # Note: Bundling this awk pipe with bash code throws an escape special-character error.

    done

    ## 8/22/2022: Completed- new and efficient way of matching and sorting peptide IDs by their sequence position ranges.
    # Finally, sort the unique reference peptide IDs by their position in the protein. This file will be used for sequence-matching in the final file!
    sort -t"," -k 2n -o Sequence_Position.csv{,}
    echo -e $msg2"Sequence-position list file completed."

    return
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function 2_ExtractData() {

    # Use awk to export the desired columns from MergedPSMs to new text file "trimmed"
    echo -e $msg1"Extracting data... "

    # this builds an associative awk array, where each element is an entire record. The columns are then printed by the pattern matched in each element.
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

    # Repeat for to get sequence/position matches from merged peptide groups 
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

    # Now store the unique IDs in a bash array for later use, using substr() to eliminate the []. characters (regex characters cannot be string literals unless escaped)
    #mseqs=($( awk -F"\t" 'NR > 1 {print $1}' $file3A | awk '{print substr($0, 4, length($0)-6)}' ))  #TODO DEPRECATE THIS


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

    echo -e $msg1"Beginning file merge..."; echo -e $msg2"Identified "$numPSMfiles" PSM output files and "$numPepGroupfiles" Peptide Group output files to merge."

    # Get the unique headers from PSM files into the new merged PSM file, and repeat for Peptide Groups.
    awk -F"\t" 'NR==1' PDoutputTextFiles/*PSMs* > $file1; awk -F"\t" 'NR==1' PDoutputTextFiles/*PeptideGroups* > $file1A
    
    # Now merge the file contents into their respective mergefiles, and remove all ". Also, build an array fnames() containing unique filenames for later use
    for filename in PDoutputTextFiles/*PSMs*; do

        # sed's -i flag edits the file 'in-place' without requiring a temporary file. This line eliminates "" from the file.
        awk -F"\t" 'NR > 1' $filename >> $file1; sed -i 's/\"//g' $file1

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
    file2=trimmed.txt; file2A=SequencePositions.txt; file3=trimmedFiltered$ext; file3A=GEE_sequencePositions.txt; file4=final$ext; temp=temp$ext

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

    ## Initialize empty output files, descriptive variables and data structures
    touch $file1; touch $file1A; touch $file2; touch $file3; touch $file4            # output files to store an image of the data at each step of analysis
    declare -i numPSMfiles; numPSMfiles=0                                            # var to keep track of data processing steps
    declare -i numPepGroupfiles; numPepGroupfiles=0                                  # ditto
    declare -i numgoalPSMs; numgoalPSMs=0                                            # ditto
    declare -i numgoalPepGrps; numgoalPepGrps=0                                      # ditto
    declare -i actual; actual=0                                                      # ditto
    declare -i actual1A; actual1A=0                                                  # ditto
    declare -a fnames=()                                                             # Array. Store the basenames of the original .raw files here
    declare -a -u mseqs=()                                                           # Array. Store the unique (uppercased) master peptide sequences across all .raw files here

    return
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
function assignPositions() {   # NOTE now that this block is a function, need to ensure that the variables it uses are visible to it (they may need to be passed to it!)
    
    ## 5/23/2022 THIS FUNCTION IS NO LONGER NEEDED HERE. DOING ASSIGNMENTS BY PRINTING UNIQUE IDS WITH THEIR CORRESPONDING POSITIONS IS ANALOGOUS TO CREATING A PERMANENT ARRAY
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
## BEGIN ANALYSIS HERE
## Check for command line arguments
if [[ $# != 0 ]]; then 
    labelkey=$(echo ${1^^})     # accepts a string, case insensitive and content insensitive. labelkey will be used to search for the label specified by this string.
    graph=$(echo ${2^^})        # accepts true or false, case insensitive. Allows user to request the final data to be graphed for them (will call a separate python script).        
    map=$(echo ${3^^})          # accepts true or false, case insensitive. Allows user to request the final data to be mapped to a protein structure in a pymol session (will call the same python script).
fi



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

## Begin (Step 0): Initialize variables and data structures; also tally the number of input files present
0_Initialize

# Calculate the number of files to process
numPSMfiles=$(ls PDoutputTextFiles/*PSMs* | wc -l); numPepGroupfiles=$(ls PDoutputTextFiles/*PeptideGroups* | wc -l)

# Populate fnames array in advance, so that successive methods can access its list of names
### WARNING this name extraction code is buggy and leads to inconsistencies in source file matching. Re-work this so it extract PD's assigned File ID instead.
for filename in PDoutputTextFiles/*PSMs*; do
    #name=$(basename $filename | awk -F"_" '{print $2}' | awk -F"-" '{print $1}'); fnames+=( $name )    ## original
    name=$(basename $filename | awk -F"_" '{print $2}' | awk -F"-" '{print $1}'); echo $name "added to fnames."; fnames+=( $name )     ## testing
done

## Step 1: Merge the output files by calling 'mergeFiles'
1_MergeFiles

## Step 2: Extract the desired information by calling ExtractData
2_ExtractData

## Step 3: Clean the data (remove redundant records and sort the unique records)
3_CleanData

## Step 4: Calculate percent labeling of each unique PSM in each file.
4_CALCULATE!


#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| 6/01/2022: Everything to this point has been tested and works correctly! |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######
#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| Script is now complete; accurately calculates abundances and %labeling |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######
#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| NEXT: Debug for ambiguous cleavage in peptide IDs; re-tool to extract |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######
#######~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~| unique peptide sequences from the PSMs in output file3 rather than from peptide groups |~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#######

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