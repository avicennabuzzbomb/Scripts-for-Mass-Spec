import os
import csv

# implement code here that checks the file extension of each file that it attempts to open' also,
# make this be filename-agnostic so that it can look into a directory (where the user has placed their file(s) to
# be analyzed), take each file, check its extension, and then use the appropriate file reader before importing its
# contents into the dictionary.

####################################################################################################
# This code stores only two columns (specified by index in line 8) from the file in the dictionary, row by row.
# In each row, the fields are paired as key:value pairs. 
with open('dictionary-tests.csv', mode='r') as infile:
    reader = csv.reader(infile)
    fruitdict = {rows[0]:rows[2] for rows in reader}

print(fruitdict) ##OUTPUT: {'Fruit': 'Clime', 'mango': 'tropical', 'apple': 'temperate', 'pear': 'temperate', 'peach': 'temperate/tropical'}

####################################################################################################
# This code sets the column header as the key and the entire list of entries as the associated value
# (as a List). This is closer to what I want.
reader = csv.DictReader(open('DSSO upperphase XLs data set.csv'))
result = {}
for row in reader:
    for column, value in row.items():
        result.setdefault(column, []).append(value)
    
print(result) ##OUTPUT: {'Fruit': ['mango', 'apple', 'pear', 'peach'], 'Seed type': ['stone', 'pit', 'pit', 'stone'], 'Clime': ['tropical', 'temperate', 'temperate', 'temperate/tropical']}


####################################################################################################
# This code sets the row header as the key and all other entries in the row as the associated value
# (as a List). In this structure, the key would be each crosslink ID if applied to the dataset.
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



