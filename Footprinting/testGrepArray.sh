# for testing how to populate an array using grep's output

file=trimmed.csv
declare -a targets=()
declare -a uniqSeq=()

targets=($(awk 'NR > 1' $file | cut -f1 -d"," | sort | uniq))

for element in "${targets[@]}"; do

   element=$(echo ${element^^})
   
   if [[ ! " ${uniqSeq[@]} " =~ " ${element} " ]]; then
   #if [[ ! " ${element} " =~ " ${uniqSeq[@]} " ]]; then
      
      uniqSeq+=( $element )   

   fi

done

echo -e "targets() contains "${#targets[@]}" elements;\n"; echo -e "uniqSeq() contains "${#uniqSeq[@]}" elements."