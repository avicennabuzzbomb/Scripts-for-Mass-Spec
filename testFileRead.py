## Create file read object to read a .txt file line by line and field by field to perform operations on field values

if __name__=='__main__':  # Run this script when invoked, instead of the modules imported into it
    import os             # for interacting with directories
    import glob           # for pattern-matching filenames

    for file in glob.glob('*.txt'):           # text files only
        readInputFile = open(file, "r")      # open file stream
        if readInputFile.mode == "r":
            print(file)
            contents = readInputFile.readlines()
            for x in contents:
                print(x)

    print("DONE")