#!/bin/bash

## 1) Iterates through all output text files (exported as-is from Proteome Discoverer) and gathers data into a single file.
## 2) Then, exports from that file only unique records.
## 3) Then, pulls only the following columns for downstream analysis: `Annotated Sequence`, `Modifications`, `Charge`, `Precursor Abundance`, `Spectrum File`, and `File ID`.
## 4) Finally, to the same file, adds a column recording modification status as yes (1) or no (0), done by pattern matching to the substring "GEE"

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Begin by storing output filenames in variables for ease of use
# specify extension of output files
ext=.csv

# output file creation
file1=mergedPSMs.txt
file1A=mergedPeptideGroups.txt
file2=trimmed$ext
file3=final$ext
temp=temp$ext

# user message flag
msg1=$(echo "\n>> "); msg2=$(echo "\t! ")

## clean up any files left over previously
if [ -e $file1 ]; then
    rm $file1
fi

if [ -e $file1A ]; then
    rm $file1A
fi

if [ -e $file2 ]; then
    rm $file2
fi

if [ -e $file3 ]; then
    rm $file3
fi

## prepare empty output files, descriptive variables and data structures
touch $file1; touch $file1A; touch $file2; touch $file3            # output files to store an image of the data at each step of analysis
declare -i numPSMfiles; numPSMfiles=0                              # var to keep track of data
declare -i numPepGroupfiles; numPepGroupfiles=0
declare -i numgoal; numgoal=0                                      # var to keep track of data
declare -i actual; actual=0                                        # var to keep track of data
declare -a fnames=()                                               # associative array to store the source file names for referencing

## First, get only one instance of the headers row (record #1)  [works for my purposes, wastes compute power though -  in future modify this so a loop isn't used]
## Also, store base filenames in an array for later sorting of the data

# count all original files to be processed
numPSMfiles=$(ls PDoutputTextFiles/*PSMs* | wc -l); numPepGroupfiles=$(ls PDoutputTextFiles/*PeptideGroups* | wc -l)
echo -e $msg1"Beginning file merge ..."; echo -e $msg2"Identified "$numPSMfiles" PSM output files and "$numPepGroupfiles" Peptide Group output files to merge."

# Merge PSMs into one file
for filename in PDoutputTextFiles/*PSMs*; do
    awk -F "\t" 'NR == 1' $filename > $file1
    # get unique portion of each filename and store it in the array fnames() for later
    name=$(tail -n1 $filename | cut -f29 | sed 's/"//g' | sed 's/.raw//g' | sed 's/[0-9][0-9][0-9][0-9][0-9][0-9]//g' | sed 's/_//g')
    fnames+=( $name )
done
for filename in PDoutputTextFiles/*PSMs*; do
    awk -F "\t" 'NR >= 2' $filename >> $file1
done

# Merge Peptide Groups into one file
for filename in PDoutputTextFiles/*PeptideGroups*; do
    awk -F "\t" 'NR == 1' $filename > $file1A
done
for filename in PDoutputTextFiles/*PeptideGroups*; do
    awk -F "\t" 'NR >= 2' $filename >> $file1A
done

## Data cleaning
# Eliminate all instances of "" in both merged files
sed -i 's/\"//g' $file1; sed -i 's/\"//g' $file1A

# count all non-header rows that should have been exported
numgoal=$(cat PDoutputTextFiles/*txt | wc -l)
numgoal="$((numgoal-numfiles))"

# count actual export numbers
actual=$(awk 'NR >= 2' $file1 | wc -l)

## check
echo -e $msg2"Found "$numgoal" records to export; exported "$actual" records."; echo -e $msg2"Merged data was saved in '"$file1"'."

## Begin Step 2

# Store column headers
Seq=$(awk -F "\t" 'NR == 1 {print $7}' $file1)
Mods=$(awk -F "\t" 'NR == 1 {print $8}' $file1)
Charge=$(awk -F "\t" 'NR == 1 {print $13}' $file1)
Filename=$(awk -F "\t" 'NR == 1 {print $29}' $file1)
Precursor=$(awk -F "\t" 'NR == 1 {print $36}' $file1)
GEE=$(echo "GEE?")

echo -e $msg1"Extracting the following values to "$file2": '"$Seq"' '"$Mods"' '"$Charge"' '"$Filename"' '"$Precursor"'\n"
echo $Seq,$Mods,$Charge,$Filename,$Precursor,$GEE > $file2

## Remove records where the Xcorr score is less than 2.0
echo -e $msg1"Now filtering for PSMs with a minimum XCorr of 2 and annotating GEE label events..."

## AWK desired columns to new CSV file "trimmed" and also add a column marking whether a GEE or hydroGEE mod is present -- SUCCESS, everything through this block works
awk -F "\t" 'BEGIN{ seq = ""; mods = ""; charge = 0; name = ""; ID = 0; xcorr = 0; abundance = 0; GEEstatus = "" };

        NR >= 2 { seq = $7; mods = $8; charge = $13; name = $29; xcorr = $32; abundance = $36;

        if ( mods ~ /GEE/ )
            GEEstatus = "GEE";
        else
            GEEstatus = "";

        if ( xcorr >= 2 )
            print seq "," mods "," charge "," name "," abundance "," GEEstatus;

        };' $file1 >> $file2

## Report the number of removed records for failing XCorr threshold   ### NOTE functions as expected (compared to manual work in Excel)
lenfile=$(awk 'NR >= 2' $file2 | wc -l)
filtered=$((actual-lenfile))

## Copy to backup file
head -n1 $file2 > backup_$file2
awk -F"," 'NR >= 2' $file2 | sort -k4n >> backup_$file2

## User message
if [ $actual != 0 ]; then
    echo -e $msg2"Found and removed "$filtered" PSMs with Xcorr score < 2 ("$lenfile" records retained)."; echo -e $msg2"These changes are also copied in 'backup_"$file2"', review as needed.\n" 
else
    echo -e $msg2"Found 0 failing PSMs; all records retained.\n"
fi

## Now remove redundant records where applicable (uses a temporary file); NOTE functions as expected (compared to manual work in Excel)
echo -e $msg1"Now checking '"$file2"' for redundant PSMs..."
touch $temp
awk -F "," 'NR == 1' $file2 > $temp
awk -F "," 'NR >= 2' $file2 | sort -k4n | uniq >> $temp
cat temp$ext > $file2

## Track the data
trim=$(awk 'NR >= 2' $file2 | wc -l)
diff=$((lenfile-trim))

## User message
if [ $trim != 0 ]; then
    echo -e $msg2"Found and removed "$diff" duplicate records ("$trim" records retained).\n"
else
    echo -e $msg2"Found 0 duplicate records; all records retained.\n"
fi

## Now calculate percent labeling of each unique PSM in each file.
echo -e $msg1"Sorting records by source .raw file ... matching unique sequences and precursor abundances..."

## Use the var 'numPSMfiles' to re-sort PSM records by source filename (which corresponds to original number of input files)
awk 'NR == 1' $file2 > $temp

for ((i = 0; i < $numPSMfiles; i++)); do     ## this requires file IDs to start with 1; otherwise it cannot match some records and ends up excluding valuable data

    identifier=$(echo ${fnames[i]})
    grep $identifier $file2 >> $temp

done

cat $temp > $file2


## Update trimmed.csv with positions in the master protein sequence ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# First, extract only unique records (exclude headers) in file1A into an array
declare -a MasterPositions=(); MasterPositions=($(awk -F "\t" 'NR > 1 {print $4, $12}' $file1A | sort -k1 | uniq))  # use for later.............


## Extract relevant data from output file 2 into separate arrays for cross-referencing by array index 
declare -a PSMs=(); PSMs=($(awk -F "," 'NR > 1 {print $1}' $file2))
declare -a Precabundances=(); Precabundances=($(awk -F "," 'NR > 1 {print $5}' $file2))
declare -a Modstatus=(); Modstatus=($(awk -F "," 'NR > 1 {print $6}' $file2))
declare -a source=(); source=($(awk -F "," 'NR > 1 {print $4}' $file2))
declare -a SourceList=(); SourceList=($(awk -F "," 'NR > 1 {print $4}' $file2 | uniq))

## Set up data structures for searching the above arrays and calculating abundances
declare -a MasterSeqs=()                    # unique, case-insensitive representations of each PSM are stored here
declare -i Mastercount=0                    # incremented during the for-loop (below) to count the number of unique PSMs that have been found
pre=0                                       # current value of precursor abundance until it gets added to the rolling value of labeled or unlabeled
modstring=""                                # stores current value of Modstatus; used to check where the value of pre should be stored
labeled=0                                   # precursor abundance stored here when a GEE mod is detected
unlabeled=0                                 # precursor abundance stored here when a GEE mod is not found
percLabeled=0                               # calculated from the final summed values stored in the labeled and unlabeled vars
declare -i index=0                          # stores the current index value of the PSMs array (incremented during loop). Used to extract corresponding values from the other arrays.
master=""                                   # stores the uppercased (case-insensitive) equivalent of the current element in the PSMs() array.
raw=""                                      # stores the source filename value corresponding to the current 

## Identify and store the unique "master" peptide sequences from the available PSMs in a separate array
for element in "${PSMs[@]}"; do
    master=$(echo ${element^^})
    if [[ ! " ${MasterSeqs[@]} " =~ " ${master} " ]]; then
        # if a unique "master" peptide sequence is not yet recorded, store it in MasterSeqs
        MasterSeqs+=( $master )
        let "++Mastercount"
    fi
done

echo -e $msg2"Identified "$Mastercount" unique master peptide sequences from among "${#PSMs[@]}" unique PSMs across "${#SourceList[@]}" .raw files."

## NOTE I believe it is necessary to begin with a while loop in here; the while loop would use the source file's name [ex., `JB1`] as a condition for doing the searches.
## ie, 'while source == JB1; do ... search commands; when source filename changes, restart the search using the new filename'
for src in "${SourceList[@]}"; do    # the outermost for loop controls how the data is sorted; it forces results to be associated with the file they originate from.
    
    break     # only until work continues on these blocks
    
    for element in "${PSMs[@]}"; do
    
        pre=$(echo ${Precabundances[index]}); modstring=$(echo ${Modstatus[index]}); raw=$(echo ${source[index]})  # <--for each source .raw filename, get all data for each PSM associated with that filename
        master=$(echo ${element^^}) # convert current PSM to its master sequence. TODO need to implement a way to detect when this master sequence changes. This needs to happen
                                    # inside of the test condition [if raw == src]

        # check that the current data's source file matches the outer loop's current selection of unique source filename
        if [[ " $raw " == " $src " ]]; then   ## TODO this should really be an until loop, using case statements. Or it should only be case statements.

            # detect the GEE modification and store the precursor value in the appropriate variable
            if [[ " $modstring " == "GEE" ]]; then
                labeled=$(echo $pre)
            else
                unlabeled=$(echo $pre)
            fi

            : # TODO: extend code here to lasso the PSM master sequence with its corresponding source file affiliation, and labeled/unlabeled abundance values

        else
            # whenever the current data source does not match the outer loop selection, skip to the next iteration of the nested for loop
            # (may want to put this condition first, rather than second to save computing power)
            #### or possibly better, can this branchpoint be redone as an until loop?
            continue
        fi
    done
done







####### UP NEXT ####### 
## 1) Create new file, "final.csv"; add headers "Unique Sequence", and then (using a loop), "FileIDxx-unlabeled", "FileIDxx-labeled", "FileIDxx-%labeled" for each file ID
## 2) Extract unique sequences to this new file;
## 3) Use each *current* value of the first column after record 1 to sum all associated abundances in a given file ID (ie, timepoint-type) (for fileID#... etc) and print them to unlabeled and labeled columns
## 4) Use those values to calculate a percent labeling, and print it to %labeled.
## 5) Underneath these calculations, print the unique sequences again, this time with columns printed containing the average and standard deviations of each percentage labeling event.
## 6) Warning: as is, code is agnostic to how the file IDs relate to each other. Need to explicitly group file IDs into distinct treatment groups (ie., J1-J4 and SP1-SP4 are always control timepoints)

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
