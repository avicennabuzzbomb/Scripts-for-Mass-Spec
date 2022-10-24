### for pulling out D/E positions from any SASA file

if [ -e *D_E.csv ]; then rm *D_E.csv*; fi

for file in SASA*ALL.csv; do
    name=$(echo $file | awk -F"_" '{print $1$2}'); name=$(echo $name"_D_E.csv"); touch $name
    awk 'NR <= 3{print $0}' $file > $name
    awk -F"," 'NR > 3{
        if ($1 ~ /D/ || $1 ~ /E/)
            print $0
    }' $file >> $name
    
done

cat $name