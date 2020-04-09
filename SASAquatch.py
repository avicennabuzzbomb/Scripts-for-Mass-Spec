## Calculate the Solvent-Accessible Surface Area (SASA) of individual residues at a time or as a group in PyMOL,
## then print to an output file.
## TODO the driver script (NunuDrives.py) will loop through THIS script
## NOTE This script works by starting a fresh PyMOL session and then loading each requested pymol session (or fetching each requested PDB ID). 
## To access individual residues one at a time, a string is made from PyMOL's fasta sequence command, and this string is copied and "cleaned" 
## so that a version containing only the amino acids is left. Then, the string is used to get the position (index) of each residue from python's enumerator.
## Whenever a residue is found in the string, its index is recorded where the first index of the string is set to '1'. Then, PyMOL selection and SASA calculation
## functions are called and written to a new csv output file. 

## NOTE this script is written to be driven by an external script which will loop it, passing a PDB ID to it with each iteration. Header and subheader will be printed on each pass.
"""
The recommended way to run PyMOL-Python scripts is by using PyMOL as the interpreter. This is supported by all versions of PyMOL, 
including the pre-compiled bundles provided by Schrödinger.

Example from a shell:

shell> pymol -cq script.py
With arguments (sys.argv becomes ["script.py", "foo", "bar"]):

shell> pymol -cq script.py -- foo bar
Example from a running PyMOL instance:

PyMOL> run script.py

## NOTE - To run a sript in PyMOL, it needs to be saved in Pymol's working directory.
## NOTE - to find a PyMOL session's working directory, simply use 'pwd' in PyMOL's command line; use cd to change directories.
## NOTE - 'reinitialize' removes all objects from this pymol session's memory (it's a reset as if the program was restarted) 
"""

## TODO TODO TODO Make this MORE agnostic by accounting for aa sequences with nonsense characters (? or -)

from pymol import cmd     # PyMOL's methods and commands; this is required for the script to use PyMOL's functions.
import pymol              # this one might be unnecessary since we specifically need pymol.cmd()
import re                 # methods for string editing
import decimal            # methods for correct rounding
import csv                # methods for handling csv file i/o
import time               # methods for tracking efficiency of the code (CPU time)
import os                 # methods for directory handling
import sys                # methods for taking command line arguments. Script's name is sys.arg[0] by default

query = str(sys.argv[0].upper()) # assign PDB ID to the first argument in the list of possible arguments, uppercase it (arg is case insensitive), and typecast to string
depth = str(sys.argv[1].upper()) # assign query type ('ALL' or aa single letter code) to the second argument. NOTE On its own, this script must take arguments
                                 # at [3] and [4] because the "pymol" and "-cq" flag used to set pymol as the interpreter are being interpreted
                                 # by this script as its arguments. If run by an external script as a subprocess, the ranges are lower [0] and [1] because the 
                                 # interpreter argument is stated when invoking the external script instead of when invoking this one.

## NOTE This Dictionary contains maximum sidechain SASA values (unrounded) for each biological amino acid; calculated with Pymol's get_area command, with dot_solvent = 1 (SASA = true)
## and with dot_density = 4 (maximum accuracy to the discrete surface area calculation); variables declared above methods are accessible to those methods! Key = residue name, value = SASA
Max_SASA = {'R' : 249.4673309,   # arg	
            'H' : 194.6491394,   # his
            'K' : 214.5114136,   # lys
            'D' : 142.6577301,   # asp
            'E' : 173.6901855,   # glu
            'S' : 110.3298187,   # ser
            'T' : 141.3123779,   # thr
            'N' : 155.2926178,   # asn
            'Q' : 184.1813354,   # gln
            'C' : 131.9754333,   # cys
            'G' : 34.49541855,   # gly
            'P' : 153.1377106,   # pro
            'A' : 91.39417267,   # ala
            'V' : 157.4189758,   # val
            'I' : 183.812027,    # ile
            'L' : 181.1264496,   # leu
            'M' : 192.694458,    # met
            'F' : 218.0111237,   # phe
            'Y' : 231.1389923,   # tyr
            'W' : 262.1514587}   # trp

## NOTE This Dictionary of 3-letter codes (required for pyMOL's selection algebra)
aa_codes = {'R' : 'lys',   # arg	
            'H' : 'his',   # his
            'K' : 'lys',   # lys
            'D' : 'asp',   # asp
            'E' : 'glu',   # glu
            'S' : 'ser',   # ser
            'T' : 'thr',   # thr
            'N' : 'asn',   # asn                # TODO TODO TODO several rhomboid protease files throw a key mismatch error because
            'Q' : 'gln',   # gln                # the fasta enumeration returns a '?' at one or more aa positions. Need to troubleshoot
            'C' : 'cys',   # cys                # how the retrieved fasta is pulling a '?' or if it's even present there. Unusual cases.
            'G' : 'gly',   # gly
            'P' : 'pro',   # pro
            'A' : 'ala',   # ala
            'V' : 'val',   # val
            'I' : 'ile',   # ile
            'L' : 'leu',   # leu
            'M' : 'met',   # met
            'F' : 'phe',   # phe
            'Y' : 'tyr',   # tyr
            'W' : 'trp'}   # trp

# NOTE / TODO When running the commands in this script individually in pymol's command line for the protease structures,
# this is returned:
'''
PyMOL>print(cmd.get_fastastr('all'))
>2xov
ERAGPVTWVMMIACVVVFIAMQILGDQEVMLWLAWPFDPTLKFEFWRYFTHALMHFSLMHILFNLLWWWY
LGGAVEKRLGSGKLIVITLISALLSGYVQQKFSGPWFGGLSGVVYALMGYVWLRGERDPQSGIYLQRGLI
IFALIWIVAGWFDLFGMSMANGAHIAGLAVGLAMAFVDSLN

PyMOL>reinitialize
PyMOL>fetch 3zmi
please wait ...
 ExecutiveLoad-Detail: Detected mmCIF
 CmdLoad: loaded as "3zmi".
PyMOL>print(cmd.get_fastastr('all'))
>3zmi
RAGPVTWVMMIACVVVFIAMQILGDQEVMLWLAWPFDPTLKFEFWRYFTHALMHFSLMHILFNLLWWWYL
GGAVEKRLGSGKLIVITLISALLSGYVQQKFSGPWFGGLSGVVYALMGYVWLRGERDPQSGIYLQRGLII
FALIWIVAGWFDLFGMSMANGAHIAGLAVGLAMAFVDSL
'''
# NOTE / TODO: this implies that there is nothing wrong with the string itself, instead that there is an issue retrieving it.

# The cutoff value for relative SASA. Residues with a value greater than 0.25 are considered solvent-exposed, otherwise are considered buried.
threshold = 0.25

# Initialize empty data stuctures for handling python's iterator stream
stored_residues = set()
residues = []
occu_vals = set()

###########################################################################################################################################################
## Helper Method: evaluate if current residue position is actually present in the structure model, and annotate accordingly                              ##
###########################################################################################################################################################
def occupancy(count):
    cmd.iterate("resi " + count, 'occu_vals.add(q)')
    if 0.0 in occu_vals: # if any atoms in the selection fail occupancy check, the selection must not have electron density
        return False
    
    else:
        return True

###########################################################################################################################################################
## Method for getting each amino acid position and using it to make a separate residue selection, and calculate and write SASA to output. (ALL RESIDUES) ##
###########################################################################################################################################################
def find_Allchain(seq, start):
    for count, ltr in enumerate(seq, start): 
        ## Typecasting count (index or aa position) to a string allows the script to use "count" in a PyMOL selection-expression.
        count = str(count)
        residue = "" + ltr + count
            
        # Method call: Qualify this selection is actually present in the structure model
        presence = occupancy(count)
        occu_vals.clear()

        # Calculations are not performed if the electron density is missing from the model; instead, their values are
        # defaulted to "Not present in structure model", "N/A", and "N/A"
        if presence == True:
            # using the reference residue side chain, calculate relative SASA by accessing from the dictionary of max SASA values.
            # Chain A must be specified every time, otherwise the value will be multiplied by the number of identical chains.

            sasa = cmd.get_area("resi " + count + " and sidechain" + " and chain A", 1, 0)
            rel_sasa = sasa / Max_SASA[ltr]
            burial = "Exposed"  # intialized arbitrarily to "Exposed"; if threshold minimum is met, does not update.
        
            if rel_sasa <= threshold:
                burial = "Buried"
        
            sasa = str(sasa)
            rel_sasa = str(rel_sasa)
        
        else:
            sasa = "Not present in structure model"
            rel_sasa = "N/A"
            burial = "N/A"

        ## Each SASA printed to console and then to output
        print(residue + " | " + sasa + " | " + rel_sasa + " | " + burial)
        current_row = [residue, sasa, rel_sasa, burial]
        writer.writerow(current_row)            
    
    return

################################################################################################################################################################
## Method for getting each lysine position and using it to make a separate residue selection, and calculate and write SASA to output. (ONLY SPECIFIED RESIDUES)#
################################################################################################################################################################
def find_res(seq, ch, start):
    for count, ltr in enumerate(seq, start):    
        if ltr == ch:
            ## Typecasting count (index or aa position) to a string allows the script to use count in a PyMOL selection-expression.
            count = str(count)
            residue = "" + ltr + count
            
            # Method call: Qualify this selection is actually present in the structure model
            presence = occupancy(count)
            occu_vals.clear()

            # Calculations are not performed if the electron density is missing from the model; instead, their values are
            # defaulted to "Not present in structure model", "N/A", and "N/A"
            if presence == True:
                # using the reference residue side chain, calculate relative SASA by accessing from the dictionary of max SASA values.
                # Chain A must be specified every time, otherwise the value will be multiplied by the number of identical chains.

                sasa = cmd.get_area("resi " + count + " and sidechain" + " and chain A", 1, 0)
                rel_sasa = sasa / Max_SASA[ltr]
                burial = "Exposed"  # intialized arbitrarily to "Exposed"; if threshold minimum is met, does not update.
        
                if rel_sasa <= threshold:
                    burial = "Buried"
        
                sasa = str(sasa)
                rel_sasa = str(rel_sasa)
        
            else:
                sasa = "Not present in structure model"
                rel_sasa = "N/A"
                burial = "N/A"
    
            ## Each SASA printed to console and then to output
            print(residue + " | " + sasa + " | " + rel_sasa + " | " + burial)
            current_row = [residue, sasa, rel_sasa, burial]
            writer.writerow(current_row)

    return
    
################################################################################################################################################
##~~~~~~~~~~~~~~~~ ___DRIVER CODE___~~~~~~~~~~~~~~~~##  NOTE Writes to the csv file with each method call.

## "Stopwatch" starts now
start_time = time.time()

header = ["SOLVENT ACCESSIBLE SURFACE AREAS OF TARGET PROTEOME"]

with open('SASA_' + query + '_' + depth + '.csv', 'w', newline = '') as file: # write output values into csv row by row, vals in separate columns
                                                                              # NOTE: the last argument, "newline = ''" prevents the file being written 
                                                                              # at every other line  
    writer = csv.writer(file, delimiter = ',')
    
    # Begin writing into the csv output file with a master header describing the job
    writer.writerow(header)

    ## start with a fresh pymol session
    cmd.reinitialize

    ## SASA settings
    cmd.set('dot_solvent', 1)  ## 1 is for solvent surface area. 0 is for total molecular surface area [default]
    cmd.set('dot_density', 4)  ## 1-4; defines quality (accuracy) of the calculation, better=more CPU

    ## When fetching:
    cmd.fetch(query, "Chain A")
    cmd.remove("het")
    unwantedHeader = ">Chain_A_A"

    ## Get the full aa sequence and count all lysines in the selection using this string
    fasta = cmd.get_fastastr('Chain A')
    clean_fasta = re.sub(unwantedHeader, "", fasta)
    clean_fasta = re.sub("\n", "", clean_fasta)     # remove Chain A header and all whitespace, leaving pyMOL's fasta string with capitalized aa single-letter codes.                                               

    # Begin writing into the csv output file. NOTE For a multiprotein job, this provides spacer and subheaders for each individual protein output
    subheader = ["Residue", "Absolute SASA", "Relative SASA", "Exposed or Buried?", "PDB ID: " + query, "FASTA:", clean_fasta]
    writer.writerow("")
    writer.writerow(subheader)

    # Iterate the residue positions into a Set named stored_residues, move the set into a List (so elements can be accessed by index) named residues, 
    # and then use List[0] to set enumerator's startpoint
    cmd.iterate("chain A", 'stored_residues.add(resv)')
    for i in stored_residues:
        residues.append(i)

    print("The first amino acid is at position " + str(residues[0]))  
    start = residues[0]

    # Run the job with the query (requested PDB ID) based on its type. NOTE / TODO Currently, this check works by looking for a single letter only (not robust, needs improving)
    if len(depth) == 1:
        ch = depth
        find_res(clean_fasta, ch, start)

    elif depth == "ALL":
        find_Allchain(clean_fasta, start)
    
## "Stopwatch" stops now; print runtime
stop_time = time.time()
print("\nTime (seconds) taken for SASA calculations: " + str(stop_time - start_time) + "\n\nOutput file: " + "SASA_" + query + "_" + depth + ".csv was saved in the working directory: " + os.getcwd())

'''
NOTE: WARNING. The following exceptions are thrown once the script finishes (it actually works! but throws this pair of exceptions):
Error: unsupported file type: 6sl6
Error: Argument processing aborted due to exception (above).
'''