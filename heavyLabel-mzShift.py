### Python script for automatically calculating the expected m/z of 15N/15N, 15N/14N, and 14N/15N labeled dipeptides identified ######
### by crosslinking mass spectrometry (XL-MS) experiments ######

if __name__=='__main__':  # Run this script when invoked, instead of the modules imported into it

    import numpy # may need for accurate math (double check this)

    """ProteomeDiscoverer's reported mass of a neutron; needed for calculating mass
    #shifts in heavy-nitrogen labeled peptides"""
    neutron = 0.997035                          

    R = 3.98814     # nitrogen mass shift in Arginine, R after heavy N labeling (4 nitrogens)
    H = 2.991105    # nitrogen mass shift in Histidine, H after heavy N labeling (3 nitrogens)
    K = 1.99407     # nitrogen mass shift in Lysine, K after heavy N labeling (2 nitrogens)
    N = 1.99407     # nitrogen mass shift in Asparagine, N after heavy N labeling (2 nitrogens)
    Q = 1.99407     # nitrogen mass shift in Glutamine, Q after heavy N labeling (2 nitrogens)
    W = 1.99407     # nitrogen mass shift in Tryptophan, W after heavy N labeling (2 nitrogens)

    seqA = ""   # sequence of peptide A of the 14N/14N dipeptide (string)
    lenA = 0    # length of peptide A; needed to loop through the string during calculation
    seqB = ""   # sequence of peptide B of the 14N/14N dipeptide (string)
    lenB = 0    # length of peptide B; needed to loop through the string during calculation
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
        assert type(z)==int, "Error, z must be an integer representing the charge of the dipeptide"

        heavyAshift = calc_15Nshift(pepA, charge)   
        heavyBshift = calc_15Nshift(pepB, charge)

        HH = heavyAshift + heavyBshift + mZ                  # fully labeled (Heavy-Heavy)
        HL = heavyAshift + mZ     # peptide A is labeled, peptide B is light (Heavy-Light)
        LH = heavyBshift + mZ     # peptide B is labeled, peptide A is light (Light-Heavy)
        LL = mZ                   # unlabeled (Light-Light)

        labeledPeps = [LL, LH, HL, HH] # store predicted m/z in a List and output the List       
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
        countRemaining = length - sum([countR, countH, countK, countN, countQ, countW])
        mZshift = (R*countR + H*countH + K*countK + N*countN + Q*countQ + W*countW + neutron*countRemaining) / charge
        print("Peptide sequence is",seq," and its m/z shift is ",mZshift)
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
        #TODO CODE BLOCK FOR OPENING AN  INPUT STREAM OBJECT HERE (read a list of CSMs into this script, call functions, write new output file)
        return

    # FUNCTION: Calls read_CSMs(); calls calc_Labeled() with parameters read from file contents; outputs pepNames and labeledPeps to a new result file.
    def ANALYZE_CSMs(): #TODO this may be redundant to read_CSMs; either bundle everything under read and make it main, or make ANALYZE the main function and have it call read.
        return

    ##____TEST PARAMS AND CODE FOR THE FUNCTION calcMzShift____##
    def testCalcs(testStrA,testStrB,test_mZ,testZ):
        dipeptide = testStrA + "-" + testStrB
        print("The dipeptide is ",dipeptide,"and its charge state is +",testZ)
        print("")

        labeledPeps = calc_Labeled(testStrA,testStrB,test_mZ,testZ)  # list of expected m/z
        pepNames = ['LL','LH','HL','HH']                             # list of label names

        for x in range(len(pepNames)):
            print("Predicted m/z of labeled peptides are",pepNames[x],": ",labeledPeps[x])
    
        return
    """
    #List of parameters for the testing function testCalcs()
    testStrA = "NETVDLEKIPIEEVFQQLK"
    testStrB = "VSKGAPEQILELAK"
    test_mZ = 978.774                                                                    
    testZ = 4
    testStrC = "IPIEEVFQQLKCSR"
    testStrD = "ELSEIAEQAKR" 
    test_mZ2 = 885.199

    #TEST 1
    print("")
    print("___TEST 1___")
    testCalcs(testStrA,testStrB,test_mZ,testZ) #TEST 1 calls testCalcs() with TEST 1 params

    #TEST 2
    print("")
    print("___TEST 2___")
    testCalcs(testStrC,testStrD,test_mZ2,testZ) #TEST 2 calls testCalcs() with TEST 2 params
    """

    response = 'y'   #default in order to enter the loop

    while response == 'Y' or response == 'y':
        seqA = input("Enter peptide sequence A: ")
        print("seqA is ",seqA)
        seqB = input("Enter peptide sequence B: ")
        print("seqB is ",seqB)
        mZ = input("Enter mZ of seqA-B: ")
        mZ = float(mZ)          #have to typecast because input defaults to str
        print("mZ is ",mZ)      
        charge = input("Enter charge, z: ")
        charge = int(charge)    #have to typecast because input defaults to str
        print("Charge: z = ",charge)

        testCalcs(seqA,seqB,mZ,charge)

        response = input("Analysis done. Enter another dipeptide?")

            













