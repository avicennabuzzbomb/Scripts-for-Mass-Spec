## testing arrays in bash
declare -a seqList=()   # initialize empty array

testfile=test1.csv
testfile2=test2.csv

## copy "Annotated Sequences" into 
touch $testfile; awk -F"," 'NR >= 2{print $1}' trimmed.csv | sort | uniq > $testfile

## appends "Annotated Sequences" one by one into array
while read p; do
    seqList+=( $p ) 
done < $testfile

#echo "Array length is "${#seqList[*]}   ## prints array length
#echo -e "Contents: \n"${seqList[*]}     ## prints array contents

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

## 2 arrays example test code:

#intialize array of master sequences and variables to store associated abundance values:
declare -a masterList=()
unlabeled=0
labeled=0

## use 1st array to populate 2nd one.
## NOTE this for-loop functions as desired:
for seq in "${seqList[@]}"; do
    masterseq=$(echo ${seq^^}) ## convert to uppercase (master sequence)

    ## now, ensure case-insensitive sequence is only stored once
    if [[ ! " ${masterList[*]} " =~ " ${masterseq} " ]]; then

        ## NOTE to use awk to update the variables for masterList(), those variables need to directly capture awk's output via echo.
        ## ie., something like...... labeled=$(echo $seq | awk '... if ( GEE ~ /GEE/ ); labeled+=abundance; print labeled')
        ## and...... unlabeled=$(echo $seq | awk ' ... if ( ! GEE ~/GEE/ ); unlabeled+=abundance; print unlabeled')
        masterList+=( $seq $labeled $unlabeled )  ## append to array with reserved spots for summed-unlabeled and summed-labeled abundances

    fi

done

## Check contents of the associative array -masterList- which contains master sequence, the abundance of labeled peptides with that sequence,
## and the abundance of unlabeled peptides matching that sequence
arrSize=$(echo ${#masterList[@]})
index=$((arrSize - 1))
echo -e "Contents of master list in indices [0,"$index"]: \n"
for ((i=0; i<=$index; ++i)); do
    echo ${masterList[i]}
done