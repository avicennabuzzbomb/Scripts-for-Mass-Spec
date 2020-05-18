### Python script for automatically calculating the expected m/z of 15N/15N, 15N/14N, and 14N/15N labeled dipeptides identified ######
### by crosslinking mass spectrometry (XL-MS) experiments ######

if __name__=='__main__':  # Run this script when invoked, instead of the modules imported into it

    import re      # for string cleanup
    import string  # substring manipulation
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
    def calc_Labeled(seq, mZ, charge):
        """
        pepA and pepB must be strings, and receive the seqences of peptides A and B respectively.
        mZ must be a float, and receives the mass:charge values (m/z) of the unlabeled dipeptide AB.
        charge must be an integer, and receives the charge (z) of the unlabeled dipeptide AB.
        """
        assert type(seq)==str, "Error, pepA must be a string representing the peptide"
        assert type(mZ)==float, "Error, mZ must be a float representing the m/z of the dipeptide"
        assert type(charge)==int, "Error, z must be an integer representing the charge of the dipeptide"

        heavyshift = calc_15Nshift(seq, charge)   # These statements calculate peptide

        H = round(heavyshift + mZ, 3)       # fully labeled (Heavy-Heavy), rounded to 3 decimal places
        L = round(mZ, 3)                                   # unlabeled (Light-Light), rounded to 3 decimal places
        #TODO: rounding may not be desireable every time; may want to make round conditional on an input - a boolean and an integer to say, "yes, round to this many places"
        #TODO: ...and that should be a different method, which is contigent on the boolean 'if response == y, boolean == true and call with (labeledPeps,x) where x is the integer
        #TODO: specifying the decimal places.

        labeledPeps = [str(L) + " (light)", str(H) + " (heavy)"]      # store predicted m/z in a List and output the List       
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

        # remove modification symbols from the string so they don't wrongly affect the amino acid calculation
        seq = re.sub('[a-z]', '', seq)  # remove lowercase (aa modifications due to MS)
        
        if "[" in seq:
            seq = seq[4:(len(seq) - 4)]

        length = len(seq)
        countR = seq.count('R')
        countH = seq.count('H')
        countK = seq.count('K')
        countN = seq.count('N') 
        countQ = seq.count('Q')
        countW = seq.count('W')

        countRemaining = length - sum([countR, countH, countK, countN, countQ, countW])   # count the backbone nitrogens of amino acids without N-containing R-groups; ignores brackets automatically (no need for string cleanup)
        mZshift = (R*countR + H*countH + K*countK + N*countN + Q*countQ + W*countW + neutron*countRemaining) / charge    # add backbone neutrons of above to the summed shifts of amino acids with N-containing R-groups
        return mZshift


###########__DRIVER_CODE__###########
# collect the rows into a List object. Each row of the input .csv is an element in this List.
#batch_file = open("PSMyeastSearch.csv", "r")
batch_file = open("PSMyeastSearch.csv", "r")
query_list = batch_file.readlines()
batch_file.close()

# create a new file that contains the heavy-labeled peptide m/z
with open('LH_PSMs.csv', 'w', newline = '') as file:
    # create a .csv file writer object                                                                         
    writer = csv.writer(file, delimiter = ',')

    # parse the whole list into individual row-lists
    rownum = 0
    Xcorr = ""
    header = []
    current_row = []
    light_heavy = []
    for element in query_list:
        if rownum == 0:   
            # this branch causes the first row (headers) to be stored separately
            # (peptide sequence [0], protein ID [2], m/z [4], charge [3], RT retention time [7], and Xcorr [9])
            current_row = element.split(",")   # slice the string-element by the "," delimiter into separate elements for the new list
            Xcorr = current_row[9]
            Xcorr = Xcorr.strip("\n")
            header = [current_row[0], current_row[2], current_row[4], current_row[3], current_row[7], Xcorr]
            print(header)
            writer.writerow(header)  
            rownum += 1
            continue    
        else:             
            # this branch gets the shifted m/z of the current peptide and stores it
            current_row = element.split(",")   # slice the string-element by the "," delimiter into separate elements for the new list
            seq = str(current_row[0])
            z = int(current_row[3])
            mZ = float(current_row[4])

            # to console
            print(seq + " | " + str(z) + " | " + str(mZ))
        
            # get m/z of fully labeled peptides
            # labeledPeps = [str(L) + " (light)", str(H) + " (heavy)"]      # store predicted m/z in a List and output the List       
            labeledPeps = calc_Labeled(seq, mZ, z)

            # to console
            print(labeledPeps)

            #This block prints the contents of the array of lists: label names, and predicted m/z to the new file.
            Xcorr = current_row[9]
            Xcorr = Xcorr.strip("\n")
            next_row = [current_row[0], current_row[2], labeledPeps[0], current_row[3], current_row[7], Xcorr]
            heavy_row = [" "," ", labeledPeps[1], " ", " ", " "]
            writer.writerow(next_row)
            writer.writerow(heavy_row)

print("\nDONE!")