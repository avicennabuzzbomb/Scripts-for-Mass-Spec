import os
import csv

# implement code here that checks the file extension of each file that it attempts to open' also,
# make this be filename-agnostic so that it can look into a directory (where the user has placed their file(s) to
# be analyzed), take each file, check its extension, and then use the appropriate file reader before importing its
# contents into the dictionary.


# ext = open(f'name','wb')        ## f strings in python 3?
#TODO check file extension by getting file's whole name as a string then distinguishing what should be done based on extension

ext = "DSSO upperphase XLs data set.txt"

if ext.endswith('.txt'):      # if statements work; once storing a filename in var 'ext' is implemented, files can be handled more appropriately
    print("TXT")
    pass

if ext.endswith('.csv'):
    print("CSV")
    pass

#TODO implement code block incorrectly formatted data from being used (that way, no data loss due to duplicate keys in the first column), 
#and return an error message to the user saying how the data should be formatted before running the script.
 
####################################################################################################
# This code stores only two columns (specified by index in line 8) from the file in the dictionary, row by row.
# In each row, the fields are paired as key:value pairs. 
"""
with open('dictionary-tests.csv', mode='r') as infile:
    reader = csv.reader(infile)
    fruitdict = {rows[0]:rows[2] for rows in reader}

print(fruitdict) ##OUTPUT: {'Fruit': 'Clime', 'mango': 'tropical', 'apple': 'temperate', 'pear': 'temperate', 'peach': 'temperate/tropical'}
"""
####################################################################################################
# This code sets the column header as the key and the entire list of entries as the associated value
# (as a List). This is closer to what I want.
"""
reader = csv.DictReader(open('DSSO upperphase XLs data set.csv'))
result = {}
for row in reader:
    for column, value in row.items():
        result.setdefault(column, []).append(value)
    
# print(result) ##OUTPUT: Entire file contents, with Keys as column headers and entire column (row by row) as the associated
              ## values in a List
positionA = result.get("Position A", "")
positionB = result.get("Position B", "")

for i in positionA: #warning: endless print loop. Need to find a way to retrieve the index of a value, and use that
    print(positionA) #to find a value for another key (example: Position A is a key, but an entry corresponding to a value)
                     #associated with it will have the same index as anything else in the same row. The index then should be able to    
                    #find any other value in the row (common index for every key) for any key.
for i in positionB:
    print(positionB)
"""
####################################################################################################
# This code sets the row header as the key and all other entries in the row as the associated value
# (as a List). In this structure, the key would be each crosslink ID if applied to the dataset.
# Problem: if the key values end up being the same, those keys and associated values are overwritten by
# any other hashes having the same key!
reader = csv.reader(open('DSSO upperphase XLs data set.csv'))
result1 = {}
for row in reader:
    key = row[0]
    if key in result1:
        # implement duplicate row handling here; data that is used as a key overwrites itself (include associated values) if
        # the key is a duplicate.
        result1[key] == row
        pass
    result1[key] = row[1:]

print(result1) ##OUTPUT: {'Fruit': ['Seed type', 'Clime'], 'mango': ['stone', 'tropical'], 'apple': ['pit', 'temperate'], 'pear': ['pit', 'temperate'], 'peach': ['stone', 'temperate/tropical']}
"""
####################################################################################################
reader = csv.reader(open('dictionary-tests.csv'))
result2 = {}
for row in reader:
    key = row[0]

    if key in result2:
        # implement duplicate row handling here
        pass
    
    result2[key] = row[1:]

print(result2)

####################################################################################################
reader3 = csv.reader(open('DSSO-ArabidopsisProteome.csv'))
result3 = {}
for row in reader3:
    key = row[0]
    result3[key] = row[1:]

print(result3)
"""