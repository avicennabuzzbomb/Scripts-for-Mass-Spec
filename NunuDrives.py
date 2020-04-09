# driver script to pass arguments to SASAquatch.py, including a list of PDB IDs

# TODO Implement a list reader?
# TODO Implement a for loop iterating through the list of PDB IDs, calling SASAquatch.py with each PDB ID as an argument
# TODO Implement an additional argument, taken from the header of the list of PDB IDs, specifying SASA parameters such as threshold, dot_density and dot_solvent;
# TODO Most importantly, have one of these arguments be 1) whether the SASA is for a specific amino acid / set of amino acids, or 2) if ALL amino acids are requested per protein.

import re
import os

# collect the arguments into a list
batch_file = open("batch.txt", "r")
query_list = batch_file.readlines()
batch_file.close()

# organize the list elements so that depth is considered separately from the query list
query_string = query_list[0]
query_list = query_string.split(',')
depth = query_list[len(query_list) - 1]
query_list = query_list[:-1]

list_index = 0
count = 0
for i in query_list:
    if count < len(query_list) - 1:        
        query = query_list[count]

        os.system("SASAquatch.py " + query + depth)

        count += 1
        print("current query: " + query)

print("depth = " + depth)
print("len(query_list) = " + str(len(query_list)))