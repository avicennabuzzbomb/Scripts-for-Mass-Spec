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
file2=trimmed$ext
file3=final$ext
temp=temp$ext

# user message flag
msg1=$(echo "\n>> "); msg2=$(echo "\t! ")

## clean up any files left over previously
if [ -e $file1 ]; then
    rm $file1
fi

if [ -e $file2 ]; then
    rm $file2
fi

if [ -e $file3 ]; then
    rm $file3
fi

## prepare empty output files, descriptive variables and data structures
touch $file1; touch $file2; touch $file3    # output files to store an image of the data at each step of analysis
declare -i numfiles; numfiles=0             # var to keep track of data
declare -i numgoal; numgoal=0               # var to keep track of data
declare -i actual; actual=0                 # var to keep track of data
declare -a fnames=()                        # associative array to store the source file names for referencing

# First, get only one instance of the headers row (record #1)  [works for my purposes, wastes compute power though -  in future modify this so a loop isn't used]
# Also, store base filenames in an array for later sorting of the data
for filename in PDoutputTextFiles/*txt; do

    awk -F "\t" 'NR == 1' $filename > $file1

    # get unique portion of each filename and store it in the array fnames() for later
    name=$(tail -n1 $filename | cut -f29 | sed 's/"//g' | sed 's/.raw//g' | sed 's/[0-9][0-9][0-9][0-9][0-9][0-9]//g' | sed 's/_//g')
    fnames+=( $name )

done


# Then get all other rows in all files (print all rows, starting with row #2)
for filename in PDoutputTextFiles/*txt; do

    # get the data to the new file
    awk -F "\t" 'NR >= 2' $filename >> $file1

done

## Data cleaning: Eliminate all instances of "" in the merged file
sed -i 's/\"//g' $file1

# count all original files that were processed
numfiles=$(ls PDoutputTextFiles/*txt | wc -l)
echo -e $msg1"Beginning file merge ..."; echo -e $msg2"Identified "$numfiles" output files to merge."

# count all non-header rows that should have been exported
numgoal=$(cat PDoutputTextFiles/*txt | wc -l)
numgoal="$((numgoal-numfiles))"

# count actual export numbers
actual=$(awk 'NR >= 2' $file1 | wc -l)

## check
echo -e $msg2"Found "$numgoal" records to export; exported "$actual" records."; echo -e $msg2"Merged data was saved in '"$file1"'."

# Begin Step 2
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
awk -F"," 'NR == 1' $file2 > $temp
awk -F"," 'NR >= 2' $file2 | sort -k4n | uniq >> $temp
cat temp$ext > $file2; rm $temp

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

## Use the var 'numfiles' to re-sort records by source filename (which corresponds to original number of input files)
awk 'NR == 1' $file2 > $temp

for ((i = 0; i < $numfiles; i++)); do     ## this requires file IDs to start with 1; otherwise it cannot match some records and ends up excluding valuable data

    identifier=$(echo ${fnames[i]})
    grep $identifier $file2 >> $temp
    #echo $identifier

done

cat $temp > $file2; rm $temp
#trim=$(awk 'NR >= 2' $file2 | wc -l)
#echo -e $msg2"$file2 has "$trim" records total."

## extract all unique "master" sequences into array to begin calculating labeling percentages (case-sensitive, so converting to uppercase is required)
## this array `MasterSeqs()` is not used for searching; instead it will be cross-referenced to the contents of `trimmed` and used to print in `final.csv`
declare -a targets=()
declare -a MasterSeqs=()
targets=($(awk 'NR > 1' $file2 | cut -f1 -d"," | sort | uniq))
for element in "${targets[@]}"; do
   element=$(echo ${element^^})
   if [[ ! " ${MasterSeqs[@]} " =~ " ${element} " ]]; then
      MasterSeqs+=( $element )   
   fi
done
echo -e $msg2"Identified "${#MasterSeqs[@]}" unique master peptide sequences."



####### UP NEXT ####### 
## 1) Create new file, "final.csv"; add headers "Unique Sequence", and then (using a loop), "FileIDxx-unlabeled", "FileIDxx-labeled", "FileIDxx-%labeled" for each file ID
## 2) Extract unique sequences to this new file;
## 3) Use each *current* value of the first column after record 1 to sum all associated abundances in a given file ID (ie, timepoint-type) (for fileID#... etc) and print them to unlabeled and labeled columns
## 4) Use those values to calculate a percent labeling, and print it to %labeled.
## 5) Underneath these calculations, print the unique sequences again, this time with columns printed containing the average and standard deviations of each percentage labeling event.
## 6) Warning: as is, code is agnostic to how the file IDs relate to each other. Need to explicitly group file IDs into distinct treatment groups (ie., J1-J4 and SP1-SP4 are always control timepoints)

