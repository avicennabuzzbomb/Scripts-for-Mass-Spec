file3=trimmedFiltered.csv
file4=fileArray.csv
file5=masterIDs.csv

## THIS IS WHAT WORKS (before END, that is just for printing) - note that the array is filled in a random order in this method.
awk -F"," 'NR > 1{ fileArray[$4] = 1 } END{ for (i in fileArray) print i}' $file3 > $file4
awk -F"," 'NR > 1{ fileArray[$1] = 1 } END{ for (i in fileArray) print i}' $file3 > $file5

# note: this method works; however arrays in awk are far more complicated than expected. It doesn't seem trivial to create them
# with specific conditions. Accessing specific arrays by their index is also complicated because they are associative (ie., index is
# an arbitrary string that can be a string of any number(s)). So accessing them is not trivial. Also, filling them with specific elements
# are stored first in a variable is also not trivial.
# .
# .
# .
# On reflecting upon these problems, I think it will be best to create reference files - each a list of unique filenames and a list of unique peptide IDs.
# Bash can rapidly create them as specified with great ease, and because they will be small, it can loop them fairly quickly.
# Each element can be stored via a nested for loop (ie, for i in filename; do; for j in masterIDs; do; ... then pass those var names into awk for
# pattern matching and cumulative addition. No need to ever touch an awk array. They are far too wacky. Additionally, unique GEE master IDs can be matched to their
# positions in the same ref file: col A = seq, col B = position; sending only current col A value to awk for matching and summing by current filename value.)

## AS OF 5/22/2022 - INSTEAD OF AWK ARRAYS, CREATE AUXILIARY FILES CONTAINING THE DESIRED UNIQUE VALUES.
## USE BASH AND/OR AWK WHERE FEASIBLE TO MAKE THOSE FILES. SPECIFICALLY WHEN MAKING THE IDs FILE: USE AWK TO GEE-PATTERN-MATCH TO ONLY THOSE IDs
## WHICH ARE LABELED (AFTER ALL, THOSE ARE THE ONLY ONES I CARE ABOUT ANYWAY; IF ALL OTHERS ARE COMPLETELY UNLABLED, THAT'S ALL I NEED TO KNOW ABOUT THOSE ONES.
## THEN, FOR MASTER IDs SPECIFICALLY, CALL THE POSITIONS FUNCTION TO APPEND THEIR POSITION-IN-MASTER-PROTEIN VALUES INTO THAT FILE. THEN, USE BASH TO FOR-LOOP
## THROUGH 1) FILENAME; 2) NESTED UNDER 1, FOR-LOOP THROUGH COL-A OF MASTER IDs FILE. THOSE CURRENT VALUES OF EACH WILL THEN BE PASSED TO AWK
## FOR *ALL* CALCULATIONS TO BE PERFORMED ON THE TRIMMED FILE.

## NOW, for the current values of the filename/masterpeptide ID , search the trimmedFiltered.csv for IDs which match the current elements
## of both files. When a GEE match is found, add its value to a variable storing GEE abundance. Else, store it into the non-GEE abundance variable.
## A third variable, %label, is equal to GEE/nonGEE+GEE * 100 and updates in real time.
## step 4: Continue until the given ID value has no more matches. At the bottom of this loop, print the current raw file, the ID value, the ID's summed GEE abundance,
## its summed nonGEE abundance, and the ID's %labeled value.
## step 5: Repeat step 4 until all ID values have been iterated and printed with their abundances.
## step 6: Repeat the nested loops until all raw files have been sampled for the unique IDs and abundances and everything has been printed. This is the calculated abundances file.
## step 7: New awk block: reading the calculated abundances file, use pattern matching with the same awk array of peptide IDs as above to iterate through the abundances file and get theh average.
## To get the average of the %labeled across files, use a counter to track the number of replicates. Sum all matched percent values, then divide at the end by the counter value.
## standard deviation on the other hand will be tricky. Perhaps instead of doing the statistics (including averaging) here, it might be time to call a script using matplotlib or R for this work.
## Because (well done!) it would be the final step!









