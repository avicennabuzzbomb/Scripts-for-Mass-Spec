### Python script for automatically calculating the expected m/z of 15N/15N, 15N/14N, and 14N/15N labeled dipeptides identified ######
### by crosslinking mass spectrometry (XL-MS) experiments ######

if __name__=='__main__':  # Run this script when invoked, instead of the modules imported into it

    import numpy   # may need for accurate math (double check this)
    import decimal # need for correct rounding TODO: build in decimal rounding, halfway cases always round up
    import csv     # tools for reading and writing into csv files

    """ProteomeDiscoverer's reported mass of a neutron; needed for calculating mass
       shifts in heavy-nitrogen labeled peptides"""
    neutron = 0.997035                          

    R = 3.98814     # nitrogen mass shift in Arginine, R after heavy N labeling (4 nitrogens)
    H = 2.991105    # nitrogen mass shift in Histidine, H after heavy N labeling (3 nitrogens)
    K = 1.99407     # nitrogen mass shift in Lysine, K after heavy N labeling (2 nitrogens)
    N = 1.99407     # nitrogen mass shift in Asparagine, N after heavy N labeling (2 nitrogens)
    Q = 1.99407     # nitrogen mass shift in Glutamine, Q after heavy N labeling (2 nitrogens)
    W = 1.99407     # nitrogen mass shift in Tryptophan, W after heavy N labeling (2 nitrogens)

    seq = ""   # sequence of peptide A of the 14N/14N dipeptide (string)
    peplength = 0    # length of peptide A; needed to loop through the string during calculation
    z = 0       # z, charge of the peptide
    mZ = 0      # m/z value of the dipeptide ID

    # FUNCTION: calculate the m/z value of all possible 15N-labeled dipeptides by adding 15N peptide shifts to original m/z value (mZ).
    # Returns a List of the labeled peptide expected m/z values.
    def calc_Labeled(pepA, pepB, mZ, charge):
        """
        pepA and pepB must be strings, and receive the seqences of peptides A and B respectively.
        mZ must be a float, and receives the mass:charge values (m/z) of the unlabeled dipeptide AB.
        charge must be an integer, and receives the charge (z) of the unlabeled dipeptide AB.
        """
        assert type(pepA)==str, "Error, pepA must be a string representing the peptide"
        assert type(pepB)==str, "Error, pepB must be a string representing the peptide"
        assert type(mZ)==float, "Error, mZ must be a float representing the m/z of the dipeptide"
        assert type(charge)==int, "Error, z must be an integer representing the charge of the dipeptide"

        heavyshift = calc_15Nshift(seq, charge)   # These statements calculate peptide

        H = round(heavyshift + mZ, 3)       # fully labeled (Heavy-Heavy), rounded to 3 decimal places
        L = round(mZ, 3)                                   # unlabeled (Light-Light), rounded to 3 decimal places
        #TODO: rounding may not be desireable every time; may want to make round conditional on an input - a boolean and an integer to say, "yes, round to this many places"
        #TODO: ...and that should be a different method, which is contigent on the boolean 'if response == y, boolean == true and call with (labeledPeps,x) where x is the integer
        #TODO: specifying the decimal places.

        labeledPeps = [L, H]      # store predicted m/z in a List and output the List       
        return labeledPeps


    # FUNCTION: peptide sequence-agnostic calculation of m/z shift in a sequence due to 15N labeling; other methods that differentiate peptides A and B should call it
    def calc_15Nshift(seq, charge):
        """
        seq is the first argument and must be a string.
        len is the second argument and must be an integer.
        charge is the third argument and must be an integer.
        """
        # Assert input argument types
        assert type(seq)==str, "Error, seq must be a string representing the peptide"
        assert type(charge)==int, "Error, charge must be an integer represent string (peptide) charge"

        length = len(seq)
        countR = seq.count('R')
        countH = seq.count('H')
        countK = seq.count('K')
        countN = seq.count('N') 
        countQ = seq.count('Q')
        countW = seq.count('W')

        countMods = 0
        cleanedStr = ""

        ## TODO COME BACK TO THIS; this needs to be its own module (or method) for peak predictor

        # remove modification symbols from the string so they don't wrongly affect the amino acid calculation
        for char in seq:
            if char.islower() == True:
                continue
            else:
                cleanedStr += char
                print(cleanedStr)

        countOther = seq.count('[') + seq.count(']') + seq.count('.')

        countRemaining = length - sum([countR, countH, countK, countN, countQ, countW, countOther])   # count the backbone nitrogens of amino acids without N-containing R-groups; ignores brackets automatically (no need for string cleanup)
        mZshift = (R*countR + H*countH + K*countK + N*countN + Q*countQ + W*countW + neutron*countRemaining) / charge    # add backbone neutrons of above to the summed shifts of amino acids with N-containing R-groups
        # print("Peptide sequence is",seq," and its m/z shift is ",mZshift)
        return mZshift

    # FUNCTION: opens file stream and reads in CSM data (currently, exported .txt files from ProteomeDiscoverer)
    def read_CSMs(): #TODO finish this code block once you learn how to read individual fields.
        """
        Reads .csv or .txt files and iteratively stores rows in a variable.
        Fields in row are pulled into separate variables, to be used as params
        for calling calc_Labeled to get mixed peptide specs. These are stored in
        list as in the test method below, and the list is used to print a new output
        file containing the analysis results
        """
        #TODO CODE BLOCK FOR OPENING AN INPUT STREAM OBJECT HERE (read a list of CSMs into this script and call functions with that imported data)

        #FIXME A full .txt table of the required data can be output from a PDresultview file. To do this, make sure no filters are applied in any of the associated tables. Example of headers below:
        #FIXME "Checked"	"Confidence"	"Max. XlinkX Score"	"Sequence A"	"Accession A"	"Position A"	"Sequence B"	"Accession B"	"Position B"	"Crosslinker"	"Crosslink Type"	"# CSMs"	"Protein Descriptions A"	"Protein Descriptions B"

        return
    
"""
    with open('DSSO-ArabidopsisProteome.csv', mode='r') as infile:
        reader = csv.reader(infile)
    #    with open('DSSO-ArabidopsisProteome.csv', mode='w') as outfile:
    #        writer = csv.writer(outfile)
        datadict = {rows[0]:rows[1] for rows in reader}

    for key, value in datadict.items():
        print(key, ' : ', value)      
"""

# collect the arguments into a list
batch_file = open("batch.txt", "r")
query_list = batch_file.readlines()
batch_file.close()

    # create a new file object Output, open an empty file with it to 'w'rite into. 'w+' means create the file to write into if it does not exist.
    #TODO csv.writer and writer.row(<val>) is superior for formatting an output .csv file; change this when there is time.
Output = open('Predicted_mixedPeaks.csv', 'w+')
print("HEADER:: Predicted m/z for light-light(LL), light-heavy(LH), heavy-light(HL), heavy-heavy(HH) dipeptides.", file = Output)  
print("", file = Output)
COUNT = 1
