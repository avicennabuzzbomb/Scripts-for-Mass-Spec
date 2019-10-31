### Python script for automatically calculating the expected m/z of 15N/15N, 15N/14N, and 14N/15N labeled dipeptides identified ######
### by crosslinking mass spectrometry (XL-MS) experiments ######

if __name__=='__main__':  # Run this script when invoked, instead of the modules imported into it

    import numpy # may need for accurate math (double check this)
    import csv   # tools for reading and writing into csv files

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

        heavyAshift = calc_15Nshift(pepA, charge)   # These statements calculate peptide
        heavyBshift = calc_15Nshift(pepB, charge)   # the expected N15 shift of each. 

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

        countBrackets = seq.count('[')
        countBrackets1 = seq.count(']')

        countRemaining = length - sum([countR, countH, countK, countN, countQ, countW, countBrackets, countBrackets1])
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
    
    # create a new file object Output, open an empty file with it to 'w'rite into.
    Output = open('Predicted_mixedPeaks.csv', 'w+')
    print("HEADER:: Predicted m/z for light-light(LL), light-heavy(LH), heavy-light(HL), heavy-heavy(HH) dipeptides.", file = Output)  
    print("", file = Output)
    COUNT = 1

    ## Block below specifically for my purposes (copy-paste directly from Thao's data, then run script over this list)
    pepA_list = ["ELSEIAEQA[K]R",
                "ELSEIAEQA[K]R",
                "TLHGLQP[K]EAVNIFPEK",
                "ELHTL[K]GHVESVVK",    
                "[K]HIVGMTGDGVNDAPALK",    # need to check this one in PD (CSM 5)
                "IPIEEVFQQL[K]CSR",        # This is wrong, m/z does not match RT (CSM 6), possible error by Thao
                "IPIEEVFQQL[K]CSR",
                "VS[K]GAPEQILELAK",
                "ELSEIAEQA[K]R",
                "LSVD[K]NLVEVFCK",
                "LSVD[K]NLVEVFCK",
                "TLHGLQP[K]EAVNIFPEK",
                "GVE[K]DQVLLFAAMASR",
                "EVHFLPFNPVD[K]R",
                "EVHFLPFNPVD[K]R",
                "EVHFLPFNPVD[K]R",
                "VLSIID[K]YAER",
                "TGTLTLN[K]LSVDK",
                "VENQDAIDAAMVGMLADP[K]EAR",
                "VENQDAIDAAMVGMLADP[K]EAR",
                "VS[K]GAPEQILELAK",
                "VS[K]GAPEQILELAK",
                "[K]ADIGIAVADATDAAR",
                "VS[K]GAPEQILELAK",
                "VS[K]GAPEQILELAK",
                "EVHFLPFNPVD[K]R",
                "ADGFAGVFPEH[K]YEIVK",
                "EVHFLPFNPVD[K]R",
                "GAPEQILELA[K]ASNDLSK",
                "EEEEEEEEEEEEEEK",         # no sequence! check PD; RT is 62.07.
                "VLSIID[K]YAER",
                "NLVEVFC[K]GVEK",
                "V[K]PSPTPDSWK",
                "MITGDQLAIG[K]ETGR",
                "IPIEEVFQQL[K]CSR",
                "NETVDLE[K]IPIEEVFQQLK",   # "positive control" test for this script, picked by Thao
                "IQIFGPN[K]LEEK",          # this is an example where the monoisotopic peak is not reliable and have to pick something in the spectra as the starting point
                "ELHTL[K]GHVESVVK",
                "VENQDAIDAAMVGMLADP[K]EAR",
                "[K]ADIGIAVADATDAAR",
                "VDQSALTGESLPVT[K]HPGQEVFSGSTCK",
                "[K]HIVGMTGDGVNDAPALK",
                "NETVDLE[K]IPIEEVFQQLK",
                "V[K]PSPTPDSWK",
                "TAFTM[K]K",
                "TLHGLQP[K]EAVNIFPEK",
                "EVHFLPFNPVD[K]R",
                "IPIEEVFQQL[K]CSR",
                "MTAIEEMAGMDVLCSD[K]TGTLTLNK",
                "NETVDLE[K]IPIEEVFQQLK",
                "EVHFLPFNPVD[K]R",
                "VDQSALTGESLPVT[K]HPGQEVFSGSTCK",
                "VENQDAIDAAMVGMLADP[K]EAR",
                "EAVNIFPE[K]GSYR",
                "MITGDQLAIG[K]ETGR",
                "LLEGDPL[K]VDQSALTGESLPVTK",
                "ELHTL[K]GHVESVVK",
                "TLHGLQP[K]EAVNIFPEK",
                "IPIEEVFQQL[K]CSR",
                "VS[K]GAPEQILELAK",
                "GAPEQILELA[K]ASNDLSK",
                "T[K]ESPGAPWEFVGLLPLFDPPR",
                "VS[K]GAPEQILELAK",
                "VENQDAIDAAMVGMLADP[K]EAR",
                "TLHGLQP[K]EAVNIFPEK",
                "MITGDQLAIG[K]ETGR",
                "VS[K]GAPEQILELAK",
                "MTAIEEMAGMDVLCSD[K]TGTLTLNK",
                "EVHFLPFNPVD[K]R",
                "VSKGAPEQILELA[K]ASNDLSKK",
                "EAVNIFPE[K]GSYR",
                "EVHFLPFNPVD[K]R",
                "EVHFLPFNPVD[K]R",
                "VS[K]GAPEQILELAK",
                "LLEGDPL[K]VDQSALTGESLPVTK",
                "EAVNIFPE[K]GSYR",
                "EAVNIFPE[K]GSYR",
                "VENQDAIDAAMVGMLADP[K]EAR",
                "MITGDQLAIG[K]ETGR",
                "MTAIEEMAGMDVLCSD[K]TGTLTLNK",    #slide 97: data is too low quality to be useful.
                "IPIEEVFQQL[K]CSR",
                "VLSIID[K]YAER",
                "EVHFLPFNPVD[K]R",
                "LLEGDPL[K]VDQSALTGESLPVTK",
                "WSEQEAAILVPGDIVSI[K]LGDIIPADAR",
                "IQIFGPN[K]LEEK",
                "ADGFAGVFPEH[K]YEIVK",
                "V[K]PSPTPDSWK",
                "V[K]PSPTPDSWK",
                "IQIFGPN[K]LEEK",
                "VENQDAIDAAMVGMLADP[K]EAR",        #slide 108, data poor quality; test different peaks.
                "LLEGDPL[K]VDQSALTGESLPVTK",
                "GAPEQILELA[K]ASNDLSK",
                "VENQDAIDAAMVGMLADP[K]EAR",
                "LSQQGAIT[K]R",
                "EVHFLPFNPVD[K]R",
                "ADGFAGVFPEH[K]YEIVK",
                "LLEGDPL[K]VDQSALTGESLPVTK",
                "VDQSALTGESLPVT[K]HPGQEVFSGSTCK",
                "V[K]PSPTPDSWK",
                "MTAIEEMAGMDVLCSD[K]TGTLTLNK",
                "MITGDQLAIG[K]ETGR",
                "VENQDAIDAAMVGMLADP[K]EAR",
                "LLEGDPL[K]VDQSALTGESLPVTK",
                "MITGDQLAIG[K]ETGR",
                "IPIEEVFQQL[K]CSR"]          #slide 123, data poor quality; test different peaks.

    pepB_list = ["LSQQGAIT[K]R",
                "LSQQGAIT[K]R",
                "EVHFLPFNPVD[K]R",
                "LSQQGAIT[K]R",        
                "[K]ADIGIAVADATDAAR",      # need to check this one in PD (CSM 5)
                "ELSEIAEQA[K]R",           # This is wrong, m/z does not match RT (CSM 6), possible error by Thao
                "IQIFGPN[K]LEEK",
                "ASNDLS[K]K",
                "DYG[K]EER",
                "LSQQGAIT[K]R",
                "LSQQGAIT[K]R",
                "ELSEIAEQA[K]R ",
                "LSVD[K]NLVEVFCK",
                "ELSEIAEQA[K]R",
                "ELSEIAEQA[K]R",
                "LSQQGAIT[K]R",
                "ASNDLS[K]K",
                "LSQQGAIT[K]R",
                "EVHFLPFNPVD[K]R",
                "GVE[K]DQVLLFAAMASR",
                "LSQQGAIT[K]R",
                "VLSIID[K]YAER",
                "ELSEIAEQA[K]R",
                "EVHFLPFNPVD[K]R",
                "ELSEIAEQA[K]R",
                "TGTLTLN[K]LSVDK",
                "[K]ADIGIAVADATDAAR",
                "ASNDLS[K]K",
                "EVHFLPFNPVD[K]R",
                "EEEEEEEEEEEEEEK",         # no sequence! check PD; RT is 62.07.
                "LSQQGAIT[K]R",
                "LSQQGAIT[K]R",
                "LSQQGAIT[K]R",
                "ELSEIAEQA[K]R",
                "EVHFLPFNPVD[K]R",
                "VS[K]GAPEQILELAK",        # "positive control" test for this script, picked by Thao
                "LSQQGAIT[K]R",            # this is an example where the monoisotopic peak is not reliable and you have to use the next (or the one after) as the starting point
                "EAVNIFPE[K]GSYR",
                "LSQQGAIT[K]R",
                "LSQQGAIT[K]R",
                "EVHFLPFNPVD[K]R",
                "LSQQGAIT[K]R",
                "ELSEIAEQA[K]R",
                "EVHFLPFNPVD[K]R",
                "DYG[K]EER",
                "EAVNIFPE[K]GSYR",
                "VLSIID[K]YAER",
                "VS[K]GAPEQILELAK",
                "VS[K]GAPEQILELAK",
                "EVHFLPFNPVD[K]R",
                "IQIFGPN[K]LEEK",          # Picked a different peak to start (vs. 634.934)
                "VS[K]GAPEQILELAK",
                "VS[K]GAPEQILELAK",
                "ELSEIAEQA[K]R",
                "VLSIID[K]YAER",
                "VS[K]GAPEQILELAK",
                "ELSEIAEQA[K]R",
                "DYG[K]EER",
                "VLSIID[K]YAER",
                "IQIFGPN[K]LEEK",
                "[K]VLSIIDK",
                "VS[K]GAPEQILELAK",
                "TGTLTLN[K]LSVDK",
                "IQIFGPN[K]LEEK",
                "LSQQGAIT[K]R",
                "VS[K]GAPEQILELAK",
                "LSVD[K]NLVEVFCK",
                "LSQQGAIT[K]R",
                "LSVD[K]NLVEVFCK",
                "VLSIID[K]YAER",
                "LSQQGAIT[K]R",
                "[K]VLSIIDK",
                "NLVEVFC[K]GVEK",
                "ASNDLSK[K]VLSIIDKYAER",
                "EVHFLPFNPVD[K]R",
                "EVHFLPFNPVD[K]R",
                "DYG[K]EER",
                "VLSIID[K]YAER",
                "LSQQGAIT[K]R",
                "LSVD[K]NLVEVFCK",     #slide 97: data is too low quality to be useful.
                "LSQQGAIT[K]R",
                "QVVPE[K]TK",
                "QVVPE[K]TK",
                "VENQDAIDAAMVGMLADP[K]EAR",
                "LSQQGAIT[K]R",
                "ES[K]LLK",
                "EVHFLPFNPVD[K]R",
                "VS[K]GAPEQILELAK",
                "IQIFGPN[K]LEEK",
                "ELSEIAEQA[K]R",
                "LSVD[K]NLVEVFCK",       #slide 108, data poor quality; test different peaks.
                "[K]ADIGIAVADATDAAR",
                "QVVPE[K]TK",
                "[K]VLSIIDK",
                "[K]VLSIIDK",
                "YEIV[K]K",
                "LSQQGAIT[K]R",
                "IQIFGPN[K]LEEK",
                "LSQQGAIT[K]R",
                "LLEGDPL[K]VDQSALTGESLPVTK",
                "EVHFLPFNPVD[K]R",
                "QVVPE[K]TK",
                "GAPEQILELA[K]ASNDLSK",
                "LSQQGAIT[K]R",
                "IQIFGPN[K]LEEK",
                "DYG[K]EER"]                #slide 123, data poor quality; test different peaks.

    mZ_list = [844.777,
               633.834,
               735.989,
               709.388,
               885.199,                    # need to check this one in PD (CSM 5)
               885.199,                    # This is wrong, m/z does not match RT (CSM 6), possible error by Thao
               830.685,
               626.329,
               776.367,
               703.122,
               703.122,
               838.690,
               861.448,
               757.887,
               757.888,
               714.878,
               776.072,
               662.864,
               1025.503,
               1059.766,
               686.130,
               737.405,
               747.881,
               648.348,
               729.138,
               786.919,
               906.205,
               524.268,
               728.578,
               1077.765,
               642.100,
               894.133,
               625.828,
               755.884,
               701.160,
               978.774,                     # "positive control" test for this script, picked by Thao
               669.363,                     # this is an example where the monoisotopic peak is not reliable and you have to use the next (or the one after) as the starting point
               811.421,
               901.448,
               704.872,
               943.666,
               771.153,
               926.480,
               749.882,
               627.291,
               718.372,
               766.152,
               847.453,
               1093.043,
               806.219,
               635.135,                     # Picked a different peak to start (vs. 634.934)
               920.667,
               996.753,
               980.821,
               764.150,
               1013.297,
               752.396,
               744.371,
               803.424,
               764.667,
               739.906,
               1024.044,
               758.168,
               980.236,
               795.679,
               808.178,
               792.428,
               997.740,
               661.942,
               987.766,
               692.860,
               668.365,
               694.904,
               987.768,
               1042.047,
               816.930,
               641.550,
               952.724,
               713.368,
               1110.537,      #slide 97: data is too low quality to be useful.
               752.147,
               598.827,
               671.605,
               1228.615,
               1059.572,
               764.086,
               733.170,
               721.133,
               704.367,
               712.374,
               1013.740,       #slide 108, data poor quality; test different peaks.
               1032.46,
               743.146,
               854.937,
               544.310,
               634.336,
               792.161,
               996.534,
               844.242,
               953.000,
               1121.795,
               669.600,
               1097.045,
               917.994,
               791.411,
               700.840]             #slide 123, data poor quality; test different peaks.

    charge_list = [3,
                   4,
                   5,
                   4,
                   4,                      # need to check this one in PD (CSM 5)
                   4,                      # This is wrong, m/z does not match RT (CSM 6), possible error by Thao
                   4,
                   4,
                   3,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   3,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   5,
                   4,
                   4,
                   4,
                   5,
                   5,
                   4,
                   4,
                   3,
                   4,
                   4,
                   5,
                   4,                   # "positive control" test for this script, picked by Thao
                   4,                   # this is an example where the monoisotopic peak is not reliable and you have to use the next (or the one after) as the starting point
                   4,
                   4,
                   4,
                   5,
                   4,
                   4,
                   4,
                   3,
                   5,
                   4,
                   4,
                   4,
                   5,                   # Picked a different peak to start (vs. 634.934)
                   5,
                   4,
                   3,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   5,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,   #slide 97: data is too low quality to be useful.
                   4,
                   4,
                   4,
                   4,
                   4,
                   3,
                   5,
                   4,
                   4,
                   4,
                   4,   #slide 108, data poor quality; test different peaks.
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   5,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4,
                   4]    #slide 123, data poor quality; test different peaks. 
                   

## MAIN (script starts)
    for x in range(len(pepA_list)):

        seqA = pepA_list[x]
        seqB = pepB_list[x]
        mZ = mZ_list[x]
        charge = charge_list[x]

        dipeptide = seqA + "-" + seqB
        print("The dipeptide is ",dipeptide,"and its charge state is +",charge)
        print("")

        labeledPeps = calc_Labeled(seqA,seqB,mZ,charge) #get list of predicted sequences from calc_Labeled
        pepNames = ["LL","LH","HL","HH"]

        #This block prints the contents of the array of lists: label names, and predicted m/z to the new file.
        print(COUNT,":","Dipeptide sequence is ",dipeptide,"and its charge state is +",charge, file = Output)
        for x in range(len(pepNames)):
            print(pepNames[x], labeledPeps[x], file = Output)
        print("", file=Output)
        COUNT += 1

    #When user stops analysis, close the file.
    Output.close()       
    print("\n\nATTENTION: New file 'Predicted_mixedPeaks.csv' containing the results has been saved in the working directory.")













