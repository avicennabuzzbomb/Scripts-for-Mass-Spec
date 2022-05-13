#!/bin/bash

## test code for while'ing through trimmed.csv and storing its values

# Get positions of unique peptide IDs
file1A=mergedPeptideGroups.txt

if [[ -e updatedTrimmed.csv ]]; then
    rm updatedTrimmed.csv
fi

touch temp1.txt 

## EXPERIMENTAL BLOCK STARTS HERE

# First, copy only the unique sequence-position pairs from mergedPeptideGroups.txt into the temp file. 
awk -F "\t" 'NR >= 2 {print $0}' $file1A | cut -f4,12 | sort -k1 | uniq > temp1.txt

# initialize the variables
col1=""; refseq=""; col2=""; col3=""; col4=""; col5=""; col6=""; position=""  # WORKS!

# store the number of rows in trimmed.csv that need a matching position value
lenTrim=$(awk '{print 0}' trimmed.csv | wc -l); echo -e "\nlength of trimmed.csv is "$lenTrim"."

echo -e "Matching sequences with their positions in the master protein..."

# use the rowcount in trimmmed.csv to store the position value, row by row, and print it appended with the other row values to a new file (updated version of trimmed.csv)
for ((i = 1; i <= $lenTrim; i++ )); do

    col1=$(awk -F "," -v i=$i 'NR==i {print $1}' trimmed.csv)

    # Explanation of awk substr: Syntax is to print from the original input string ($0) starting at char position 4 (non-0 indexed, so starting at '.')
    # and continuing for the *original* string's length, minus 6 chars (ie, the flanking [char] combined is 6 chars)
    refseq=$(echo ${col1^^} | cut -f1 | awk '{print substr($0, 4, length($0)-6)}'); echo -e "\nrefseq currently equals: "$refseq
    col2=$(awk -F "," -v i=$i 'NR==i {print $2}' trimmed.csv)
    col3=$(awk -F "," -v i=$i 'NR==i {print $3}' trimmed.csv)
    col4=$(awk -F "," -v i=$i 'NR==i {print $4}' trimmed.csv)
    col5=$(awk -F "," -v i=$i 'NR==i {print $5}' trimmed.csv)
    col6=$(awk -F "," -v i=$i 'NR==i {print $6}' trimmed.csv)

    if [[ ! $i -gt 1 ]]; then
        position=$(echo "Position in Master Protein")
    else
        position=$(grep -w -m1 $refseq temp1.txt | cut -f2)     #-m1 forces grep to stop after the first match
        echo "Sequence "$col1" was matched to positions "$position
    fi

    echo "$col1,$position,$col2,$col3,$col4,$col5,$col6" >> updatedTrimmed.csv      ### everything here works exactly as desired. Incorporate this block into the main script!!

done