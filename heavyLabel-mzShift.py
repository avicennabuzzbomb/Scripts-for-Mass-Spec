### Python script for automatically calculating the expected m/z of 15N/15N, 15N/14N, and 14N/15N labeled dipeptides identified ######
### by crosslinking mass spectrometry (XL-MS) experiments ######

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
    assert type(pepA)==str, "Error, seq must be a string representing the peptide"
    assert type(pepB)==str, "Error, seq must be a string representing the peptide"
    assert type(mZ)==float, "Error, seq must be a float representing the m/z of the dipeptide"
    assert type(z)==int, "Error, seq must be an integer representing the charge of the dipeptide"

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
    #mZshift = mshift / charge

    # print("Peptide mass shift is ",mshift)
    print("Peptide sequence is",seq," and its m/z shift is ",mZshift)

    return mZshift



##____TEST PARAMS FOR calcMzShift____##
#FIXME combine TEST 1 and TEST 2 into one method
def testCalcs(testStrA,testStrB,test_mZ,testZ):
    dipeptide = testStrA + "-" + testStrB
    print("The dipeptide is ",dipeptide,"and its charge state is +",testZ)
    print("")

    labeledPeps = calc_Labeled(testStrA,testStrB,test_mZ,testZ)  # list of expected m/z
    pepNames = ['LL','LH','HL','HH']                             # list of label names

    for x in range(len(pepNames)):
        print("Predicted m/z of labeled peptides are",pepNames[x],": ",labeledPeps[x])
    
    return

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


#TODO CODE BLOCK FOR OPENING AN  INPUT STREAM OBJECT HERE (read a list of CSMs into this script, call functions, write new output file)














