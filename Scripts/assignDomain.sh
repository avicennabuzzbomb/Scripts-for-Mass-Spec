####Fetch a list of crosslink spectral matches {peptide sequences} and their linked residues, and assign the Xlink to a domain####
####For a single protein, enter a name, its domain names, and the range (indices) of those domains####
####For now, test with AHA-defined domains as a preset, and match crosslinks down a list to their domain assignment####
####Input: a .csv list of spectral matches and residue matches, ex) imaginary crosslink between 650 ARRNG[K]FAIRYNG and 345 [K]NNSPECNKR####

####Domain names need to be stored in variables that are associated with numerical ranges####

DATE=$(date)

domain1="Actuator"                                                                                                                                                                                                                                      #Actuator: residues 12–57 and 129–233
domain2="Phosphorylation"                                                                                                                                                                                                                               #Phosphorylation: residues 308–337 and 489–625
domain3="Nucleotide-binding"                                                                                                                                                                                                                            #Nucleotide-binding: residues 338–488
domain4="C-terminal"                                                                                                                                                                                                                                    #C-terminal: residues 821-924
domain5="cytoplasmic loop, TM6/TM7"                                                                                                                                                                                                                     #cytoplasmic loop, TM6/TM7: residue 696



for filename in test/*
do
    mkdir testResult #creates an output folder to store the new files
    analysis=$(basename -s ".csv" "test/$filename" | sed s/".csv"//g) 
    echo $analysis
    echo "$DATE"
    echo "$DATE" > testResult/$analysis"_domainAssignments".csv
    echo ""
    echo "The domains in AHA2 are: "
    echo "$domain1, $domain2, $domain3, $domain4, $domain5"
    echo ""
    echo "Domain assignments of crosslinked residues in AHA2:" >> testResult/$analysis"_domainAssignments".csv
    echo "Available input files are: "
    ls test

    #use awk to look for resA and resB values (columns $1 and $2) and loop through them to compare to a numerical range [domain]
    #awk can print individual fields at a time - in one awk statement.

    #loop through each line record, cut irrelevant columns, and use if-then statements to compare lysine numbers [if x, then domain]
    #initialize variables resA and resB, then begin reading at line 2.

    # 'resA=0' stores integer value of the crosslinked lysine [K] in peptide A
    # 'resB=0' stores integer value of the crosslinked lysine [K] in peptide B

    awk 'BEGIN{ resA = 0; resB = 0; domainName = "" }; NE >= 2 { resA = $1; resB = $2 }; if ( resA >= 12 && =< 57 || resA >= 129 && 233 ) domainName = "Actuator"; else if ( resA >= 308 && resA <= 337 || resA <= 489 && resA <= 625 ) domainName = "Phosphorylation"; else if ( resA >= 338 && resA <= 488 ) domainName = Nucleotide-binding; else if (resA >= 821 && resA <=924) domainName = "C-terminal"; else if ( resA == 696 ) domainName = "Cytoplasmic loop, TM6/TM7"} END{ print ResA domainName };' testDomainResult

   # cat $filename | cut -d "," -f 1,2,3,4,5,6,7 >> testResult/$analysis"_domainAssignments".csv
   # grep $filename ""

done
    
#next: read file contents for Xlink patterns (column 2)

