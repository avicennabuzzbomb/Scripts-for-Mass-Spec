import fileinput  # Class of file reader-writer methods
import re         # Class of string-replacement operations

if __name__ == "__main__":
    pass

# Reads all the contents of a file into a List;
# prints the new string into a new file, but scrubs numbers and whitespace.
# This is meant for "cleaning" DNA and amino acid sequences for easy processing by
# 3rd party applications (like TOPCONS or ProtParam) that require it but do not do it
# for you.

# REQUIRED: your input file must be a .txt file named 'fasta'

# initialize empty string
fasta = ""

# first get all lines from file by storing each line in a new List of strings called 'lines'
with open('fasta.txt', 'r') as f:
    lines = f.readlines()

# remove all digits and whitespace from the List 'lines'
pattern = '[0-9]'
spaces = ' '
lines = [re.sub(pattern, '', i) for i in lines] 
lines = [re.sub(spaces, '', i) for i in lines]

# combine the List elements from 'lines' into a single string; line breaks are stripped during concatenation.
for line in lines:
    fasta += line.strip()

with open('cleanedFasta.txt', 'w') as f:
    f.writelines(fasta)

# close the input stream
f.close()
