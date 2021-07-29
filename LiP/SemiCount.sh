## bash shell script for analyzing peptide sequence data and annotating a peptide ID as 'tryptic' or 'semitryptic'

out=$(echo "LiP_OUTPUT")
input=$(echo "LiP_INPUT")
declare -i detectedFiles
detectedFiles=0

if [ ! -e $out ]; then
    mkdir $out
    printf "Output directory was not found; created a folder named $out. Your output files will be stored there.\n"
fi

let detectedFiles=$( ls $input/*txt | wc -l )

if [ $detectedFiles -gt 0 ]; then
    printf "\nFound $detectedFiles text files in the input folder $input. Beginning annotation:\n_______________________________________________________________________________\n"
    for filename in $input/*.txt
    do
        # get the basename for copying the original filenames, and set headers to use in the output file
	    name=$(basename -s ".txt" " ./$input/$filename"); name=$(echo $name"_annotated"); echo "Tallying peptides in filename $name..."
        header=$(echo "Annotated Sequence")
        header2=$(echo "Modifications")
        header3=$(echo "NtermCut")
        header4=$(echo "CtermCut")
        header5=$(echo "Cleavage Type")
        header6=$(echo "PSMs")

        # create the new output file
        echo $header,$header2,$header3,$header4,$header5,$header6 > $out/$name.csv

        # Prints all desired fields, including N- and C-terminal cutsites. Classifies peptides as "semi" or "tryptic", then prints all with updated cleavage "type";
        # accounts for peptides that can only get cut once (ie, peptides with one end being the Protein's N- or C-terminus, designated as [-])
        awk -F "\t" 'BEGIN{ seq = ""; mods = ""; type = ""; PSMs = ""; Nterm = ""; Cterm = ""; C_noCut = "" };
    
                NR >= 2 { seq = $4; mods = $5; type = "BUKKAW"; PSMs = $8; Nterm = substr(seq, 3, 1); Cterm = substr(seq, (length(seq)-5), 1); C_noCut = substr(seq, (length(seq)-2), 1);
            
                if (( Nterm == "R" || Nterm == "K" || Nterm == "-" ) && ( Cterm == "R" || Cterm == "K" || C_noCut == "-" ))
                    type = "Tryptic";

                else
                    type = "Semitryptic";

                print seq "," mods "," Nterm "," Cterm "," type "," PSMs };
            
                END{ print "" };' $filename >> $out/$name.csv

    # additional steps pending in END action (ie, calculate the proportions of each sequence out of the total PSMs to get a weighted tally of peptide sequence in the result; print at bottom of file)!
    # use grep (patterns matching the substring "Semi") and pipe to "wc -l" to count the number of semitryptics. Also do this for all lines minus the headers. Store outputs
    # in separate variables and use to calculate % semitryptic. Print these outputs at the bottom of the current output file
    
        let semiCount=$( grep "Semitryptic" $out/$name.csv | wc -l ); let tryptCount=$( grep "Tryptic" $out/$name.csv | wc -l )
    
        printf "\nUnique Semitryptics","Unique Tryptics","\n"$semiCount,$tryptCount >> $out/$name.csv

        # grep "Semitryptic" $out/$name.csv | awk 'BEGIN{ sumPSMs += $6 } END{ print sumPSMs };

    done

    # TODO Once output files are done, generate a summary output file that shows a table comparing/describing all of the other new output files
    # Statement to user confirming outputs
    printf "\n\n...done!\n\n\nYour output file(s) are located in an output folder named $out, shown here:\n______________________________________________________________________________________________\n"; ls; printf "\n\nAnd your output file(s) in $out/ are listed below:\n______________________________________________________________________________________________\n";ls $out

else

echo "Warning - no input folders were found in $input. Make sure your input files are 1) .txt files and 2) copied into $input before starting."

fi