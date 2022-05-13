## testing arrays in bash
declare -a -u seqList=()   # initialize empty indexed array; the additional -u forces its elements to become uppercased as they are appended

## stream text directly into array
seqList=($( awk 'NR>=2' trimmed.csv | cut -f1 -d"," ))

## appends "Annotated Sequences" one by one into array

#while read p; do
#    seqList+=( $p ) 
#done < trimmed.csv

declare -i index
index=${#seqList[@]}

## check array contents
echo -e "Contents of master list in indices [0,"$index"]: \n"
for ((i=0; i<=$index; ++i)); do
    echo ${seqList[i]}
done