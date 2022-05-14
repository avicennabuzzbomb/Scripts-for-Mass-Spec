### develop a block to handle column extraction by name alone using awk
file1=test.txt
untrimmed=untrimmed.txt

if [[ -e untrimmed.txt ]]; then
    rm untrimmed.txt
fi

awk -F"\t" '
NR==1 {
    for (i=1; i<=NF; i++) {
        f[$i] = i
    }
}
{ print $(f["Annotated Sequence"]) "\t" $(f["Modifications"]) "\t" $(f["Charge"]) "\t" $(f["Spectrum File"]) "\t" $(f["Precursor Abundance"]) "\t" $(f["XCorr"]) "\t" $(f["Rank"])}
' $file1 >> $untrimmed
