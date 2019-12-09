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

# first get all lines from file by storing each line in a new List of strings called 'lines'
with open('fasta.txt', 'r') as f:
    lines = f.readlines()
    #print(type(lines))

# remove spaces from the List 'lines'
lines = [line.replace(' ', '') for line in lines]

# remove all digits from the List 'lines'
pattern = '[0-9]'
lines = [re.sub(pattern, '', i) for i in lines] 

# finally, write the List 'lines' into a new file
with open('cleanedFasta.txt', 'w') as f:
    f.writelines(lines)

# close the input stream
f.close()
