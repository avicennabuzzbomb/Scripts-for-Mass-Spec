### Python script for automatically calculating the expected m/z of 15N/15N, 15N/14N, and 14N/15N labeled dipeptides identified ######
### by crosslinking mass spectrometry (XL-MS) experiments ######

if __name__=='__main__':  # Run this script when invoked, instead of the modules imported into it

    import numpy   # may need for accurate math (double check this)
    import decimal # need for correct rounding TODO: build in decimal rounding, halfway cases always round up
    import csv     # tools for reading and writing into csv files

    # ProteomeDiscoverer's reported mass of a neutron; needed for calculating mass shifts in heavy-nitrogen labeled peptides
    neutron = 0.997035                          

    # Amino acid symbol and additional mass due to sidechains containing extra nitrogen (beyond the 1 default nitrogen, present in the backbone)
    R = 3.98814     # nitrogen mass shift in Arginine, R after heavy N labeling (extra mass from 4 additional nitrogens)
    H = 2.991105    # nitrogen mass shift in Histidine, H after heavy N labeling (extra mass from 3 additional nitrogens)
    K = 1.99407     # nitrogen mass shift in Lysine, K after heavy N labeling (extra mass from 2 additional nitrogens)
    N = 1.99407     # nitrogen mass shift in Asparagine, N after heavy N labeling (extra mass from 2 additional nitrogens)
    Q = 1.99407     # nitrogen mass shift in Glutamine, Q after heavy N labeling (extra mass from 2 additional nitrogens)
    W = 1.99407     # nitrogen mass shift in Tryptophan, W after heavy N labeling (extra mass from 2 additional nitrogens)

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
        assert type(charge)==int, "Error, z must be an integer representing the charge of the dipeptide"

        heavyAshift = calc_15Nshift(pepA, charge)   # These statements calculate peptide
        heavyBshift = calc_15Nshift(pepB, charge)   # the expected N15 shift of each. 

        HH = round(heavyAshift + heavyBshift + mZ, 3)       # fully labeled (Heavy-Heavy), rounded to 3 decimal places
        HL = round(heavyAshift + mZ, 3)                     # peptide A is labeled, peptide B is light (Heavy-Light), rounded to 3 decimal places
        LH = round(heavyBshift + mZ, 3)                     # peptide B is labeled, peptide A is light (Light-Heavy), rounded to 3 decimal places
        LL = round(mZ, 3)                                   # unlabeled (Light-Light), rounded to 3 decimal places
        #TODO: rounding may not be desireable every time; may want to make round conditional on an input - a boolean and an integer to say, "yes, round to this many places"
        #TODO: ...and that should be a different method, which is contigent on the boolean 'if response == y, boolean == true and call with (labeledPeps,x) where x is the integer
        #TODO: specifying the decimal places.

        labeledPeps = [LL, LH, HL, HH]      # store predicted m/z in a List and output the List       
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

        countRemaining = length - sum([countR, countH, countK, countN, countQ, countW, countBrackets, countBrackets1])   # count the backbone nitrogens of amino acids without N-containing R-groups; ignores brackets automatically (no need for string cleanup)
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
    # create a new file object Output, open an empty file with it to 'w'rite into. 'w+' means create the file to write into if it does not exist.
    #TODO csv.writer and writer.row(<val>) is superior for formatting an output .csv file; change this when there is time.
Output = open('Predicted_mixedPeaks.csv', 'w+')
print("HEADER:: Predicted m/z for light-light(LL), light-heavy(LH), heavy-light(HL), heavy-heavy(HH) dipeptides.", file = Output)  
print("", file = Output)
COUNT = 1

    ## Block below specifically for my purposes (copy-paste directly from Thao's data, then run script over this list)
pepA_list = ["ELSEIAEQA[K]R",                      #2       # corrected; 1 and 2 were duplicate IDs.
            "TLHGLQP[K]EAVNIFPEK",                 #3
            "ELHTL[K]GHVESVVK",                    #4
            "[K]HIVGMTGDGVNDAPALK",                #5
                "IPIEEVFQQL[K]CSR",                #6
                "IPIEEVFQQL[K]CSR",                #7
                "VS[K]GAPEQILELAK",                #8       # corrected
                "ELSEIAEQA[K]R",                   #9       # <-- left off after the 9th CSM (10/31/2019)
                "LSVD[K]NLVEVFCK",                 #10      #10 and 11 were actually identical; correct on 11/02/2019
                "TLHGLQP[K]EAVNIFPEK",             #11
                "GVE[K]DQVLLFAAMASR",              #12
                "EVHFLPFNPVD[K]R",                 #13
                "EVHFLPFNPVD[K]R",                 #15
                "VLSIID[K]YAER",                   #16
                "TGTLTLN[K]LSVDK",                 #17
                "VENQDAIDAAMVGMLADP[K]EAR",        #18
                "VENQDAIDAAMVGMLADP[K]EAR",        #19
                "VS[K]GAPEQILELAK",                #20
                "VS[K]GAPEQILELAK",                #21
                "[K]ADIGIAVADATDAAR",              #22    
                "VS[K]GAPEQILELAK",                #23
                "VS[K]GAPEQILELAK",                #24
                "EVHFLPFNPVD[K]R",                 #25
                "ADGFAGVFPEH[K]YEIVK",             #26
                "EVHFLPFNPVD[K]R",                 #27
                "GAPEQILELA[K]ASNDLSK",            #28
                "MTAIEEMAGMDVLCSD[K]TGTLTLNK",     #29    # no sequence! check PD; RT is 62.07. Corrected 11/01. <-- left off after the 30th CSM (11/01/2019)
                "VLSIID[K]YAER",                   #30
                "NLVEVFC[K]GVEK",                  #31
                "V[K]PSPTPDSWK",                   #32
                "MITGDQLAIG[K]ETGR",               #33 
                "IPIEEVFQQL[K]CSR",                #34
                "NETVDLE[K]IPIEEVFQQLK",           #35    # "positive control" test for this script, picked by Thao
                "IQIFGPN[K]LEEK",                  #36    # this is an example where the monoisotopic peak is not reliable and have to pick something in the spectra as the starting point
                "ELHTL[K]GHVESVVK",                #37
                "VENQDAIDAAMVGMLADP[K]EAR",        #38
                "[K]ADIGIAVADATDAAR",              #39  
                "VDQSALTGESLPVT[K]HPGQEVFSGSTCK",  #40
                "[K]HIVGMTGDGVNDAPALK",            #41
                "NETVDLE[K]IPIEEVFQQLK",           #42
                "V[K]PSPTPDSWK",                   #43
                "TAFTM[K]K",                       #44
                "TLHGLQP[K]EAVNIFPEK",             #45
                "EVHFLPFNPVD[K]R",                 #46
                "IPIEEVFQQL[K]CSR",                #47
                "MTAIEEMAGMDVLCSD[K]TGTLTLNK",     #48
                "NETVDLE[K]IPIEEVFQQLK",           #49
                "EVHFLPFNPVD[K]R",                 #50
                "VDQSALTGESLPVT[K]HPGQEVFSGSTCK",  #51
                "VENQDAIDAAMVGMLADP[K]EAR",        #52
                "EAVNIFPE[K]GSYR",                 #53
                "MITGDQLAIG[K]ETGR",               #54
                "LLEGDPL[K]VDQSALTGESLPVTK",       #55
                "ELHTL[K]GHVESVVK",                #56
                "TLHGLQP[K]EAVNIFPEK",             #57
                "IPIEEVFQQL[K]CSR",                #58
                "VS[K]GAPEQILELAK",                #59
                "GAPEQILELA[K]ASNDLSK",            #60
                "T[K]ESPGAPWEFVGLLPLFDPPR",        #61
                "VS[K]GAPEQILELAK",                #62
                "VENQDAIDAAMVGMLADP[K]EAR",        #63
                "TLHGLQP[K]EAVNIFPEK",             #64
                "MITGDQLAIG[K]ETGR",               #65
                "VS[K]GAPEQILELAK",                #66
                "MTAIEEMAGMDVLCSD[K]TGTLTLNK",     #67
                "EVHFLPFNPVD[K]R",                 #68
                "VSKGAPEQILELA[K]ASNDLSKK",        #69
                "EAVNIFPE[K]GSYR",                 #70
                "EVHFLPFNPVD[K]R",                 #71
                "EVHFLPFNPVD[K]R",                 #72
                "VS[K]GAPEQILELAK",                #73
                "LLEGDPL[K]VDQSALTGESLPVTK",       #74
                "EAVNIFPE[K]GSYR",                 #75
                "EAVNIFPE[K]GSYR",                 #76
                "VENQDAIDAAMVGMLADP[K]EAR",        #77
                "MITGDQLAIG[K]ETGR",               #78
                "MTAIEEMAGMDVLCSD[K]TGTLTLNK",     #79
                "IPIEEVFQQL[K]CSR",                #80
                "VLSIID[K]YAER",                   #81
                "EVHFLPFNPVD[K]R",                 #82
                "LLEGDPL[K]VDQSALTGESLPVTK",       #83
                "WSEQEAAILVPGDIVSI[K]LGDIIPADAR",  #84
                "IQIFGPN[K]LEEK",                  #85
                "ADGFAGVFPEH[K]YEIVK",             #86
                "V[K]PSPTPDSWK",                   #87
                "V[K]PSPTPDSWK",                   #88
                "IQIFGPN[K]LEEK",                  #89
                "VENQDAIDAAMVGMLADP[K]EAR",        #90
                "LLEGDPL[K]VDQSALTGESLPVTK",       #91
                "GAPEQILELA[K]ASNDLSK",            #92
                "VENQDAIDAAMVGMLADP[K]EAR",        #93
                "LSQQGAIT[K]R",                    #94
                "EVHFLPFNPVD[K]R",                 #95
                "ADGFAGVFPEH[K]YEIVK",             #96
                "LLEGDPL[K]VDQSALTGESLPVTK",       #97
                "VDQSALTGESLPVT[K]HPGQEVFSGSTCK",  #98
                "V[K]PSPTPDSWK",                   #99
                "MTAIEEMAGMDVLCSD[K]TGTLTLNK",     #100
                "MITGDQLAIG[K]ETGR",               #101
                "VENQDAIDAAMVGMLADP[K]EAR",        #102
                "LLEGDPL[K]VDQSALTGESLPVTK",       #103
                "MITGDQLAIG[K]ETGR",               #104
                "IPIEEVFQQL[K]CSR"]                #105

pepB_list = ["LSQQGAIT[K]R",                        #2
                "EVHFLPFNPVD[K]R",                  #3
                "LSQQGAIT[K]R",                     #4
                "[K]ADIGIAVADATDAAR",               #5   
                "ELSEIAEQA[K]R",                    #6   
                "IQIFGPN[K]LEEK",                   #7
                "ASNDLS[K]K",                       #8
                "DYG[K]EER",                        #9
                "LSQQGAIT[K]R",                     #10    
                "ELSEIAEQA[K]R ",                   #11
                "LSVD[K]NLVEVFCK",                  #12
                "ELSEIAEQA[K]R",                    #13
                "LSQQGAIT[K]R",                     #15
                "ASNDLS[K]K",                       #16
                "LSQQGAIT[K]R",                     #17
                "EVHFLPFNPVD[K]R",                  #18
                "GVE[K]DQVLLFAAMASR",               #19
                "LSQQGAIT[K]R",                     #20
                "VLSIID[K]YAER",                    #21
                "ELSEIAEQA[K]R",                    #22
                "EVHFLPFNPVD[K]R",                  #23
                "ELSEIAEQA[K]R",                    #24
                "TGTLTLN[K]LSVDK",                  #25
                "[K]ADIGIAVADATDAAR",               #26
                "ASNDLS[K]K",                       #27
                "EVHFLPFNPVD[K]R",                  #28
                "NLVEVFC[K]GVEK",                   #29  # no sequence! check PD; RT is 62.07. Corrected 11/01. <-- left off after the 30th CSM (11/01/2019).
                "LSQQGAIT[K]R",                     #30
                "LSQQGAIT[K]R",                     #31
                "LSQQGAIT[K]R",                     #32
                "ELSEIAEQA[K]R",                    #33
                "EVHFLPFNPVD[K]R",                  #34
                "VS[K]GAPEQILELAK",                 #35   # "positive control" test for this script, picked by Thao
                "LSQQGAIT[K]R",                     #36   # this is an example where the monoisotopic peak is not reliable and you have to use the next (or the one after) as the starting point
                "EAVNIFPE[K]GSYR",                  #37
                "LSQQGAIT[K]R",                     #38
                "LSQQGAIT[K]R",                     #39
                "EVHFLPFNPVD[K]R",                  #40
                "LSQQGAIT[K]R",                     #41
                "ELSEIAEQA[K]R",                    #42
                "EVHFLPFNPVD[K]R",                  #43
                "DYG[K]EER",                        #44
                "EAVNIFPE[K]GSYR",                  #45
                "VLSIID[K]YAER",                    #46
                "VS[K]GAPEQILELAK",                 #47
                "VS[K]GAPEQILELAK",                 #48
                "EVHFLPFNPVD[K]R",                  #49
                "IQIFGPN[K]LEEK",                   #50   # Picked a different peak to start (vs. 634.934)
                "VS[K]GAPEQILELAK",                 #51
                "VS[K]GAPEQILELAK",                 #52
                "ELSEIAEQA[K]R",                    #53
                "VLSIID[K]YAER",                    #54
                "VS[K]GAPEQILELAK",                 #55
                "ELSEIAEQA[K]R",                    #56
                "DYG[K]EER",                        #57
                "VLSIID[K]YAER",                    #58
                "IQIFGPN[K]LEEK",                   #59
                "[K]VLSIIDK",                       #60
                "VS[K]GAPEQILELAK",                 #61
                "TGTLTLN[K]LSVDK",                  #62
                "IQIFGPN[K]LEEK",                   #63
                "LSQQGAIT[K]R",                     #64
                "VS[K]GAPEQILELAK",                 #65
                "LSVD[K]NLVEVFCK",                  #66
                "LSQQGAIT[K]R",                     #67
                "LSVD[K]NLVEVFCK",                  #68
                "VLSIID[K]YAER",                    #69
                "LSQQGAIT[K]R",                     #70
                "[K]VLSIIDK",                       #71
                "NLVEVFC[K]GVEK",                   #72
                "ASNDLSK[K]VLSIIDKYAER",            #73
                "EVHFLPFNPVD[K]R",                  #74
                "EVHFLPFNPVD[K]R",                  #75
                "DYG[K]EER",                        #76
                "VLSIID[K]YAER",                    #77
                "LSQQGAIT[K]R",                     #78
                "LSVD[K]NLVEVFCK",                  #79
                "LSQQGAIT[K]R",                     #80
                "QVVPE[K]TK",                       #81
                "QVVPE[K]TK",                       #82
                "VENQDAIDAAMVGMLADP[K]EAR",         #83
                "LSQQGAIT[K]R",                     #84
                "ES[K]LLK",                         #85
                "EVHFLPFNPVD[K]R",                  #86
                "VS[K]GAPEQILELAK",                 #87
                "IQIFGPN[K]LEEK",                   #88
                "ELSEIAEQA[K]R",                    #89
                "LSVD[K]NLVEVFCK",                  #90
                "[K]ADIGIAVADATDAAR",               #91
                "QVVPE[K]TK",                       #92
                "[K]VLSIIDK",                       #93
                "[K]VLSIIDK",                       #94
                "YEIV[K]K",                         #95
                "LSQQGAIT[K]R",                     #96
                "IQIFGPN[K]LEEK",                   #97
                "LSQQGAIT[K]R",                     #98
                "LLEGDPL[K]VDQSALTGESLPVTK",        #99
                "EVHFLPFNPVD[K]R",                  #100
                "QVVPE[K]TK",                       #101
                "GAPEQILELA[K]ASNDLSK",             #102
                "LSQQGAIT[K]R",                     #103
                "IQIFGPN[K]LEEK",                   #104
                "DYG[K]EER"]                        #105

mZ_list = [633.834,                                 #2
               735.989,                             #3
               709.388,                             #4
               885.199,                             #5        
               795.154,                             #6        
               830.685,                             #7         
               626.335,                             #8        
               776.367,                             #9
               703.373,                             #10 
               838.690,                             #11
               861.448,                             #12
               758.389,                             #13
               714.878,                             #15
               776.072,                             #16
               662.864,                             #17
               1025.753,                            #18
               1059.766,                            #19
               686.130,                             #20
               737.654,                             #21
               747.881,                             #22
               648.348,                             #23
               729.138,                             #24
               786.919,                             #25
               906.205,                             #26
               524.268,                             #27
               728.980,                             #28
               1077.765,                            #29
               642.100,                             #30
               894.133,                             #31
               625.828,                             #32
               755.884,                             #33
               701.160,                             #34
               979.027,                             #35    # "positive control" test for this script, picked by Thao
               669.363,                             #36    # this is an example where the monoisotopic peak is not reliable and you have to use the next (or the one after) as the starting point
               811.421,                             #37
               901.448,                             #38
               704.872,                             #39
               943.666,                             #40
               771.153,                             #41
               926.728,                             #42
               750.133,                             #43
               627.291,                             #44
               718.372,                             #45
               766.404,                             #46
               847.453,                             #47
               1093.293,                            #48
               806.219,                             #49
               635.536,                             #50     # Picked a different peak to start (vs. 634.934)
               920.868,                             #51
               997.255,                             #52
               981.155,                             #53
               764.150,                             #54
               1013.547,                            #55
               752.396,                             #56
               744.371,                             #57
               803.424,                             #58
               764.667,                             #59
               739.906,                             #60
               1024.546,                            #61
               758.670,                             #62
               980.236,                             #63
               795.930,                             #64
               808.428,                             #65
               798.428,                             #66
               997.740,                             #67       
               662.143,                             #68
               987.766,                             #69
               693.360,                             #70
               668.365,                             #71
               795.410,                             #72
               987.768,                             #73
               1042.553,                            #74
               816.930,                             #75
               641.550,                             #76
               952.724,                             #77
               712.866,                             #78
               1110.047,                            #79
               752.147,                             #80
               598.827,                             #81
               671.605,                             #82
               1228.615,                            #83
               1059.572,                            #84
               764.086,                             #85
               733.170,                             #86
               721.383,                             #87
               704.870,                             #88
               712.625,                             #89
               1013.989,                            #90     #slide 108, data poor quality; test different peaks.
               1033.043,                            #91
               743.146,                             #92
               855.188,                             #93
               544.310,                             #94
               634.586,                             #95
               792.408,                             #96
               996.785,                             #97
               844.625,                             #98
               953.500,                             #99
               1121.795,                            #100
               669.600,                             #101
               1097.045,                            #102
               918.245,                             #103
               791.411,                             #104
               700.840]                             #105    #slide 123, data poor quality; test different peaks.

charge_list = [4,                                   #2
                   5,                               #3
                   4,                               #4
                   4,                               #5                
                   4,                               #6                  
                   4,                               #7
                   4,                               #8            # corrected
                   3,                               #9
                   4,                               #10
                   4,                               #11
                   4,                               #12
                   4,                               #13
                   4,                               #15
                   3,                               #16
                   4,                               #17
                   4,                               #18
                   4,                               #19
                   4,                               #20
                   4,                               #21
                   4,                               #22
                   5,                               #23
                   4,                               #24
                   4,                               #25
                   4,                               #26
                   5,                               #27
                   5,                               #28
                   4,                               #29
                   4,                               #30
                   3,                               #31
                   4,                               #32
                   4,                               #33
                   5,                               #34
                   4,                               #35              # "positive control" test for this script, picked by Thao
                   4,                               #36              # this is an example where the monoisotopic peak is not reliable and you have to use the next (or the one after) as the starting point
                   4,                               #37
                   4,                               #38
                   4,                               #39
                   5,                               #40
                   4,                               #41
                   4,                               #42
                   4,                               #43
                   3,                               #44
                   5,                               #45
                   4,                               #46
                   4,                               #47
                   4,                               #48
                   5,                               #49            
                   5,                               #50            # Picked a different peak to start (vs. 634.934)
                   5,                               #51
                   4,                               #52
                   3,                               #53
                   4,                               #54
                   4,                               #55
                   4,                               #56
                   4,                               #57
                   4,                               #58
                   4,                               #59
                   4,                               #60
                   4,                               #61
                   4,                               #62
                   4,                               #63
                   4,                               #64
                   4,                               #65
                   4,                               #66
                   4,                               #67
                   5,                               #68
                   4,                               #69
                   4,                               #70
                   4,                               #71
                   4,                               #72
                   4,                               #73
                   4,                               #74
                   4,                               #75
                   4,                               #76
                   4,                               #77
                   4,                               #78
                   4,                               #79
                   4,                               #80   
                   4,                               #81
                   4,                               #82
                   4,                               #83
                   4,                               #84
                   3,                               #85
                   5,                               #86
                   4,                               #87
                   4,                               #88
                   4,                               #89
                   4,                               #90   #slide 108, data poor quality; test different peaks.
                   4,                               #91
                   4,                               #92
                   4,                               #93
                   4,                               #94
                   4,                               #95
                   4,                               #96
                   4,                               #97
                   5,                               #98
                   4,                               #99
                   4,                               #100
                   4,                               #101
                   4,                               #102
                   4,                               #103
                   4,                               #104
                   4]                               #105      #slide 123, data poor quality; test different peaks. 
                   

## MAIN (script starts)
for x in range(len(pepA_list)):

    seqA = pepA_list[x]
    seqB = pepB_list[x]
    mZ = mZ_list[x]
    charge = charge_list[x]
    dipeptide = seqA + "-" + seqB
    '''
    print("The dipeptide is ",dipeptide,"and its charge state is +",charge)
    print("")
    '''
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