## Create file read object to read a .txt file line by line and field by field to perform operations on field values

if __name__=='__main__':  # Run this script when invoked, instead of the modules imported into it
    import os             # for interacting with directories
    import glob           # for pattern-matching filenames

#    for file in os.listdir('test/'):             # all files in the test/ subdirectory
    for file in glob.glob('*.txt'):           # text files only
        readInputFile = open(file, "r")      # open file stream
        if readInputFile.mode == "r":
            print(file)
            contents = readInputFile.readlines()
            for x in contents:               #TODO: currently this works (can print contents of files in the working directory)
                print(x)                     #TODO: need to specify columns (for .csv files) to print to new output file; good for grabbing data of interest from the raw output.
                                             #TODO: but also need to read in individual fields (one row at a time) and perform heavyLabel-mzShift.py operations on them,
                                             #TODO: to be printed to a final output file. Maybe also call assignDomains.sh?
    print("DONE, needs additional code")