## Calculate the Solvent-Accessible Surface Area (SASA) of individual residues at a time or as a group in PyMOL,
## then print to an output file.

## This script works by starting a fresh PyMOL session and then loading each requested pymol session (or fetching each requested PDB ID). 
## To access individual residues one at a time, a string is made from PyMOL's fasta sequence command, and this string is copied and "cleaned" 
## so that a version containing only the amino acids is left. Then, the string is used to get the position (index) of each residue from python's enumerator.
## Whenever a residue is found in the string, its index is recorded where the first index of the string is set to '1'. Then, PyMOL selection and SASA calculation
## functions are called and written to a new csv output file.
"""
The recommended way to run PyMOL-Python scripts is by using PyMOL as the interpreter. This is supported by all versions of PyMOL, including the pre-compiled bundles provided by Schrödinger.

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

## TODO TODO TODO Add arguments and a docstring.

from pymol import cmd     # PyMOL's methods and commands; this is required for the script to use PyMOL's functions.
import pymol              # this one might be unnecessary since we specifically need pymol.cmd()
import re                 # methods for string editing
import decimal            # methods for correct rounding
import csv                # methods for handling csv file i/o
import time               # methods for tracking efficiency of the code (CPU time)

## Dictionary containing maximum sidechain SASA values (unrounded) for each biological amino acid; calculated with Pymol's get_area command, with dot_solvent = 1 (SASA = true)
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

## Dictionary of 3-letter codes (required for pyMOL's selection algebra)
aa_codes = {'R' : 'lys',   # arg	
            'H' : 'his',   # his
            'K' : 'lys',   # lys
            'D' : 'asp',   # asp
            'E' : 'glu',   # glu
            'S' : 'ser',   # ser
            'T' : 'thr',   # thr
            'N' : 'asn',   # asn
            'Q' : 'gln',   # gln
            'C' : 'cys',   # cys
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

header = ["Residue", "Absolute SASA", "Relative SASA", "Exposed or Buried?"]
threshold = 0.25   # the cutoff value for relative SASA. Residues with a value greater than 0.25 are considered solvent-exposed, otherwise are considered buried.

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
    print("Length of fasta string (characters) is " + str(num_aa))

    for count, ltr in enumerate(seq, start):  # 6sl6 is truncated but also has no electron density until res 89; need to implement code to detect problems like this
                                               # and account for them at runtime.
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
            sasa = cmd.get_area("resi " + count + " and chain A", 1, 0)
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
        print("\n" + residue + " | " + sasa + " | " + rel_sasa + " | " + burial)
        current_row = [residue, sasa, rel_sasa, burial]
        writer.writerow(current_row)            
    
    return
################################################################################################################################################################
## Method for getting each lysine position and using it to make a separate residue selection, and calculate and write SASA to output. (ONLY SPECIFIED RESIDUES)#
################################################################################################################################################################

## seq: string input (amino acid sequence); ch: the amino acid letter of choice; thr: relative SASA cutoff value for determining exposure
def find_res(seq, ch, start):
    print("Length of fasta string (characters) is " + str(num_aa))

    for count, ltr in enumerate(seq, start):    
        if ltr == ch:
            ## Typecasting count (index or aa position) to a string allows the script to use count in a PyMOL selection-expression.
            count = str(count)
            residue = "" + ltr + count
            ## selection algebra for picking one specific residue at a time and performing operation(s) on it.
            ## note that this specific line of code can be universally useful for making discrete selections.
            cmd.select("sele", "resn " + ltr + " and resi " + count + " and chain A")

            ## SASA and relative SASA are calculated for each residue <ch> in the string fasta and then typecast to string so it can be printed easily.
            # using the reference residue side chain, calculate relative SASA by accessing from the dictionary of max SASA values.
            sasa = cmd.get_area("sele", 1, 0)
            rel_sasa = sasa / Max_SASA[ch]
            # threshold check: exposure is "buried" or "exposed".
            burial = "Exposed"  # intialized arbitrarily to "Exposed"; if threshold minimum is met, does not update.
            
            if rel_sasa <= threshold:
                burial = "Buried"
            
            sasa = str(sasa)
            rel_sasa = str(rel_sasa)
            ## Each SASA printed to console and then to output
            print("\n" + residue + " | " + sasa + " | " + rel_sasa + " | " + burial)
            current_row = [residue, sasa, rel_sasa, burial]
            writer.writerow(current_row)

    return
################################################################################################################################################

##~~~~~~~~~~~~~~~~ ___DRIVER CODE___~~~~~~~~~~~~~~~~##
## Query the user for job parameters:
print("Welcome to SASAfrass. \nThis script automates the calculation the solvent-accessible surface area (SASA) of a target list of amino acid residue sidechains from a list of protein models.\nThe calculation uses PyMOL's built-in SASA algorithm, which requires that the protein has a structural model file that can be read by PyMOL.")
#TODO print("Does your job contain more than one protein?") #if yes, provide the name of the list file or manually enter a list of PDB IDs or filenames
print("")

jobquery = str(input("Enter the PDB ID of the protein model you would like to analyze, then press ENTER:\n"))
jobquery = jobquery.upper() # fully capitalize the input to make it case-insensitive
print("You entered " + jobquery + ".\n")

jobquery1 = str(input("Currently this script can calculate each residue's SASA, or it can calculate the SASA for each occurrence of an amino acid.\nTo get all residues, type 'all' and press ENTER;\nTo get only one kind of residue, type the single letter amino acid code (for example, 'K' for lysines) and press ENTER.\n"))
jobquery1 = jobquery1.upper() # fully capitalize the input to make it case-insensitive
print("You entered " + jobquery1 + ".\n")

## "Stopwatch" starts now
start_time = time.time()

with open('SASA_' + jobquery + '.csv', 'w', newline = '') as file: # write output values into csv row by row, vals in separate columns
    writer = csv.writer(file, delimiter = ',')
    writer.writerow(header)

    ## start with a fresh pymol session
    cmd.reinitialize

    ## SASA settings
    cmd.set('dot_solvent', 1)  ## 1 is for solvent surface area. 0 is for total molecular surface area [default]
    cmd.set('dot_density', 4)  ## 1-4; defines quality (accuracy) of the calculation, better=more CPU

    ## When fetching:
    cmd.fetch(jobquery, "Chain A")
    unwantedHeader = ">Chain_A_A"

    ## Get the full aa sequence and count all lysines in the selection using this string
    fasta = cmd.get_fastastr('Chain A')
    clean_fasta = re.sub(unwantedHeader, "", fasta)
    clean_fasta = re.sub("\n", "", clean_fasta)     # remove Chain A header and all whitespace, leaving pyMOL's fasta string with capitalized aa letters.                                               
    num_aa = len(clean_fasta)

    print("FASTA sequence is: \n" + fasta)
    print("'Clean' FASTA sequence is: \n" + clean_fasta)

    # Iterate the residue positions into a Set, move the set into a List, and then use List[0] to set enumerator's startpoint
    cmd.iterate("chain A", 'stored_residues.add(resv)')

    for i in stored_residues:
        residues.append(i)

    start = residues[0]

    if len(jobquery1) == 1:
        ch = jobquery1
        find_res(clean_fasta, ch, start)

    elif jobquery1 == "ALL":
        find_Allchain(clean_fasta, start)

## "Stopwatch" stops now; print runtime
stop_time = time.time()
print("Time (seconds) taken for SASA calculations: " + str(stop_time - start_time))