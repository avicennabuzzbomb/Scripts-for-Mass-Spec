# This is a "driver" script to pass arguments to other scripts, such as SASAquatch.py. It reads the contents of a batch text file,
# (for example, batch.txt) such as a list of PDB IDs. It then parses the list into a pair of variables, 'query' and 'depth', and loops through the list
# for each possible value of 'query'. At each iteration in the loop, it passes those variables to the script it is calling (currently, SASAquatch.py).
# NOTE that this command line call inside the loop is programmable, so this script could in principle be used to drive ANY script written for PyMOL.

import re
import os

script_name = "SASAquatch.py"
print("PyMOL script: " + script_name)

# collect the arguments into a list
batch_file = open("batch.txt", "r")
query_list = batch_file.readlines()
batch_file.close()

# organize the list elements so that depth is considered separately from the query list
query_string = query_list[0]
query_list = query_string.split(',')
depth = query_list[len(query_list) - 1]
del query_list[-1]

# loops through the query list, calling the named script with the current values of 'query' and 'depth'
count = 0
for i in query_list:         
    os.system("pymol -c " + script_name + " " + i + " " + depth)
    print("current query: " + i)
    count += 1

print("depth = " + depth)
print("len(query_list) = " + str(len(query_list)))