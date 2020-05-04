### Use the binomial coefficient n choose k to make a list of all unique possible lysine-lysine crosslinks, ignoring spatial constraints.
### Use a fasta input file, downloadable from RCSB PDB.

import csv      # for writing the table output
import re       # for string operations
import time     # for demonstrating speed of this process

start = time.time()

# first get all lines from file by storing each line in a new List of strings called 'fasta'. Use fasta sequences deposited on uniprot whenever possible
# (tair has an additional chunk of amino acids which are not in the uniprot sequence and are not in any structure models of AHA2)
with open('AHA2_Uniprot.txt', 'r') as f:
    fasta = f.readlines()

# concatenate list contents to a new string, then remove all non-letter characters
fastastr = ""
for i in fasta:
    fastastr += i

fastastr = re.sub('\t', '', fastastr)
fastastr = re.sub(' ', '', fastastr)
fastastr = re.sub('\n', '', fastastr)
fastastr = re.sub('[0-9]', '', fastastr)

# a list to store the lysine positions
lysines = []

for i in range(len(fastastr)):
    count = i + 1
    if fastastr[i] == 'K':
        current = fastastr[i] + str(count) 
        lysines.append(current)

## open a new csv file to store the unique crosslinks table
header = ["Position A", "Position B"]

with open('UniqueCrosslinks.csv', 'w', newline = '') as file:
    writer = csv.writer(file, delimiter = ',')
    writer.writerow(header)

    ## make a new list of possible lysine combinations
    crosslinks = []
    spacer = " x "

    for A in range(len(lysines)):
        for B in range(len(lysines)):
            # treat crosslinked positions [A][B] and [B][A] as identical
            pair = lysines[A] + spacer + lysines[B]
            altpair = lysines[B] + spacer + lysines[A]
            # exclude intermolecular crosslinks
            if pair == altpair:
                continue
            # include only unique crosslinks
            if pair not in crosslinks:
                if altpair not in crosslinks:
                    crosslinks.append(pair)
                    print(pair)
                    current_row = [lysines[A], lysines[B]]
                    writer.writerow(current_row)
                    current_row.clear()

    print("Number of lysines in sequence: " + str(len(lysines)))
    print("Number of possible unique crosslinks: " + str(len(crosslinks)))

stop = time.time()
time = stop - start
print("\nTime spent assembling table of unique crosslinks: " + str(time) + " seconds.")