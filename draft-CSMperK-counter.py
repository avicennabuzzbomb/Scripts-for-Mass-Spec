##TODO This is inspired by a problem I encountered while analyzing the number of CSMs of each crosslinkable lysine in AHA2

##TODO My "human algorithm" required:
##          1. Reviewing several Xlink datasets on AHA2 (each from a different sample prep method and/or instrument method)
##          2. Then, making a list of unique crosslink IDs, where Peptide A x Peptide B and 
##             Peptide B x Peptide A (if A = A and B = B) are equivalent
##          3. Then, summing all CSMs per Xlink ID per occurrence of A x B and B x A (if A = A and B = B) to have a 
##             CSM tally per unique ID
##          4. Finally, counting the occurrence of each unique crosslinking lysine in the list of unique crosslink IDs,
##             and further summing all those CSMs to get an approximation of overall crosslinking frequency of each lysine
##             (useful for comparing to SASA values)
##
##TODO Implement code that does this automatically for any number of lists in txt, csv, delimited or xlsx format, and outputs
# the list of those lysines with their associated CSM tally; then finally merge with SASA values?

print("Finish manual data, then write this when you have time.")

# TODO Implement proposed algorithm:
# Once file is opened and being read, go through the list from top to bottom, and set a variable
# 'unique_ID' equal to the first new value A x B for that cell. Then, create the reverse peptide (B x A)
# from the contents of that cell or variable using parameter expansion, and search the list and all files for that; use a counter
# to verify uniqueness. Also, something like "for unique_ID in file; CSMtotal = 0; if unique_ID == "A x B"; CSMcount = {int #CSMs in associated cell}; CSMtotal += CSMcount
#  #second if statement, not nested; if unique_ID == "B x A"; CSMcount = {int #CSMs in associated cell}; CSMtotal += CSMcount; ##end the for loop, but nest inside of another for loop
# to repeat this action on each file that is present in some directory or named in some batch.