#!/bin/bash

name=""         # filename
value=0         # current value of precursor abundance
sumGEE=0        # working sum of GEE-labeled precursor abundances
sumOther=0      # working sum of precursor abundances without GEE label
sumAll=0        # all the PSMs detected
rowcount=0      # the number of rows in the current filename; needed to iterate through each file agnostic to its length
perc=0          # percent of GEE label abundance in the file
label=""        # string that stores a label value to search for a specific pattern. 

#echo "Filename","GEE-PSM abundance","Other-PSM abundance","Percent GEE abundance"
echo "Filename","GEE-PSM abundance","Other-PSM abundance","Percent GEE abundance" > outputSums.csv

for file in PDoutputTextFiles/*PSMs*; do
    
    # store current file's name and number of data rows (excluding headers)
    name=$(tail -n1 $file | cut -f29 | sed 's/"//g' | sed 's/.raw//g' | sed 's/[0-9][0-9][0-9][0-9][0-9][0-9]//g' | sed 's/_//g'); echo "Name is "$name
    rowcount=$(awk -F "\t" 'NR > 1 {print $1}' $file | wc -l); echo "File "$name" has "$rowcount "rows."

    # use current value of i to set where awk begins reading values
    for ((i=2; i<=$rowcount; i++)); do

        # get current precursor abundance and label
        value=$(awk -F"\t" -v i=$i 'NR==i {print $34}' $file); echo "Found precursor value "$value" at i="$i"."
        label=$(awk -F"\t" -v i=$i 'NR==i {print $8}' $file); echo "Found label description "$label" ..."

        if [[ $label == "*GEE*" ]]; then
            #sumGEE=$(awk -v sumGEE=$sumGEE -v value=$value 'BEGIN {print (sumGEE + value)}'); echo -e "\nValue of sumGEE is now "$sumGEE
            sumGEE=$(echo "scale=2; $sumGEE + $value" | bc -l); echo -e "\tThis has GEE label(s)."
        else
            #sumOther=$(awk -v sumOther=$sumOther -v value=$value 'BEGIN {print (sumOther + value)}'); echo "Value of sumOther is now "$sumOther
            sumOther=$(echo "scale=2; $sumOther + $value" | bc -l); echo -e "\tNo label of interest found."
        fi

        #sumAll=$(awk -v sumAll=$sumAll -v value=$value 'BEGIN {print (sumAll + value)}'); echo "Value of sumAll is now "$sumAll
        sumAll=$(echo "scale=2; $sumAll + $value" | bc -l)

        #break ## temporary until this code works

    done
    
    #perc=$(awk -v sumGEE=$sumGEE -v sumAll=$sumAll 'BEGIN {print (100 * (sumGEE / sumAll))}')
	perc=$(echo "scale=2; $sumGEE / $sumAll" | bc -l)													
	perc=$(echo "$perc * 100" | bc -l | xargs printf "%.*f\n" "$4")

    #echo $name,$sumGEE,$sumOther,$perc
    echo $name,$sumGEE,$sumOther,$perc >> outputSums.csv
    

    #break ## temporary until this code works

done