### DUMP FAULTY CODE HERE FOR ISOLATED TROUBLESHOOTING! ###




## otherwise, next block to add can be developed here. From main script, do this:
### FINAL SEGMENT: DO ABUNDANCE MATCHING, SUMMATION AND % LABELING CALCULATIONS. USE *ONLY* UNIQUE, GEE-LABELED IDs AS THE SEARCH KEYS
######################################################################################################################################
## Extract relevant data from output file 3 into separate arrays for cross-referencing by array index. Printing by column number here is safe because
## the output file made by this script are always formatted the same way regardless of how they started.


## Put it all in one awk statement, similar to above (line 183). Think about it, why bother creating multiple arrays when you coud just awk the entire file
## (which itself is structured as an array anyway)
## DEPRECATE EVERYTHING FROM THIS POINT, AND MAKE IT INTO ONE BIG AWK STATEMENT INSTEAD!
## OPTIONAL: BUILD A SMALL TEST FILE STRUCTURED IN THE SAME WAY AS FILE3. Or, just use File 3 itself.
## step 1 ('BEGIN'): populate awk array with unique peptide IDs that were labeled with GEE.
## step 2: populate a second awk array with unique raw file names.
## step 3: awk loop - for element in raw array, iterate through the peptide ID array. For the current value of the peptide ID array, search the trimmedFiltered.csv .
## for IDs which match the current element of both arrays. When a GEE match is found, add its value to a variable storing GEE abundance. Else, store add into a non-GEE variable.
## A third variable, %label, is equal to GEE/nonGEE+GEE * 100 and updates in real time.
## step 4: Continue until the given ID value has no more matches. At the bottom of this loop, print the current raw file, the ID value, the ID's summed GEE abundance,
## its summed nonGEE abundance, and the ID's %labeled value.
## step 5: Repeat step 4 until all ID values have been iterated and printed with their abundances.
## step 6: Repeat the nested loops until all raw files have been sampled for the unique IDs and abundances and everything has been printed. This is the calculated abundances file.
## step 7: New awk block: reading the calculated abundances file, use pattern matching with the same awk array of peptide IDs as above to iterate through the abundances file and get theh average.
## To get the average of the %labeled across files, use a counter to track the number of replicates. Sum all matched percent values, then divide at the end by the counter value.
## standard deviation on the other hand will be tricky. Perhaps instead of doing the statistics (including averaging) here, it might be time to call a script using matplotlib or R for this work.
## Because (well done!) it would be the final step!









