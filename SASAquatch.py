## Calculate the Solvent-Accessible Surface Area (SASA) of individual residues at a time or as a group in PyMOL,
## then print to an output file.
## NOTE the driver script (NunuDrives.py) will loop through this script
## NOTE This script works by starting a fresh PyMOL session and then loading each requested pymol session (or fetching each requested PDB ID). 
## To access individual residues one at a time, a string is made from PyMOL's fasta sequence command, and this string is copied and "cleaned" 
## so that a version containing only the amino acids is left. Then, the string is used to get the position (index) of each residue from python's enumerator.
## Whenever a residue is found in the string, its index is recorded where the first index of the string is set to '1'. Then, PyMOL selection and SASA calculation
## functions are called and written to a new csv output file. 

"""
The recommended way to run PyMOL-Python scripts is by using PyMOL as the interpreter. This is supported by all versions of PyMOL, 
including the pre-compiled bundles provided by Schrödinger.

Example from a shell:

shell> pymol -c script.py
With arguments (sys.argv becomes ["script.py", "foo", "bar"]):

shell> pymol -c script.py foo bar

Example from a running PyMOL instance:
PyMOL> run script.py

## NOTE - To run a script from the PyMOL command line, it needs to be saved in Pymol's working directory.
## NOTE - PyMOL's command line commands are Unix-like (same as Bash and Windows command line)
## NOTE - 'reinitialize' removes all objects from this pymol session's memory (it's a reset as if the program was restarted) 
"""

## NOTE: "Bad" residue selections are detected if cmd.count_atom("sele") == 0. The moment a residue is detected like this,
## that entire protein is logged in an error file as failed due to bad selection algebra, and the loop skips to the next query. 

from pymol import cmd     # PyMOL's methods and commands; this is required for the script to use PyMOL's functions.
import pymol              # this one might be unnecessary since we specifically need pymol.cmd()
import re                 # methods for string editing
import decimal            # methods for correct rounding
import csv                # methods for handling csv file i/o
import time               # methods for tracking efficiency of the code (CPU time)
import os                 # methods for directory handling
import sys                # methods for taking command line arguments. Script's name is sys.arg[0] by default when -c flag is used

query = str(sys.argv[3].upper()) # query: assign PDB ID to the first argument in the list of possible arguments, uppercase it (arg is case insensitive), and typecast to string
depth = str(sys.argv[4].upper()) # depth: assign query type ('ALL' or aa single letter code) to the second argument. 
                                 # NOTE This script must take arguments at [3] and [4] because the intrepreter argument "pymol", the "-c" flag, and the script name are being interpreted
                                 # by this script as initial arguments [0], [1], and [2]. Unclear why this is the case. TODO eventually fix that, but for now everything works!

print("query: " + query)
print("depth: " + depth)

## NOTE This Dictionary contains the maximum calculated SASA (dot_solvent = 1, dot_density = 4) value and 3-letter residue code as a tuple value
## associated with each single-character amino acid letter code (the key) currently being analyzed.
AA_attributes = {'R' : ('arg', 249.4673309),   # arginine
                 'H' : ('his', 194.6491394),   # histidine
                 'K' : ('lys', 214.5114136),   # lysine
                 'D' : ('asp', 142.6577301),   # aspartic acid
                 'E' : ('glu', 173.6901855),   # glutamic acid
                 'S' : ('ser', 110.3298187),   # serine
                 'T' : ('thr', 141.3123779),   # threonine
                 'N' : ('asn', 155.2926178),   # aspargine
                 'Q' : ('gln', 184.1813354),   # glutamine
                 'C' : ('cys', 131.9754333),   # cysteine
                 'G' : ('gly', 34.49541855),   # glycine
                 'P' : ('pro', 153.1377106),   # proline
                 'A' : ('ala', 91.39417267),   # alanine
                 'V' : ('val', 157.4189758),   # valine
                 'I' : ('ile', 183.812027),    # isoleucine
                 'L' : ('leu', 181.1264496),   # leucine
                 'M' : ('met', 192.694458),    # methionine
                 'F' : ('phe', 218.0111237),   # phenylalanine
                 'Y' : ('tyr', 231.1389923),   # tyrosine
                 'W' : ('trp', 262.1514587)}   # tryptophan

# The cutoff value for relative SASA. Residues with a value greater than 0.25 are considered solvent-exposed, otherwise are considered buried.
threshold = 0.25

# Initialize empty data stuctures for handling python's iterator stream and generally crappy API
stored_residues = set()
residues = []
occu_vals = set()

#############################################################################################################################################################
## Helper Method: evaluate if electron density exists at the current residue position is actually present in the structure model, and annotate accordingly ##
#############################################################################################################################################################
def occupancy(count):
    ##TODO implement an if branch to check for a valid selection; if invalid, print to error file and skip SASA calculation.
    ##TODO if this script writes the file, each bad query will have its own file; if Nunu writes the file, it will be one file that is a
    ##TODO list of the bad queries.

    
    '''
    TODO
    if cmd.count_atoms("sele") == 0:
        <do not call this method>; else, call this method>
        <write this part carefully; each method forks off before this occupancy check occurs, so each calculation method will need to use this code block>
        <which also means this code block needs to contain the call to occupancy() in each SASA method>
    TODO
    '''


    cmd.iterate("resi " + count, 'occu_vals.add(q)') # 'q' is occupancy, or the presence of electron density in the map
    if 0.0 in occu_vals: # if any atoms in the selection fail occupancy check, the selection must not have electron density
        return False 
    else:
        return True
###########################################################################################################################################################
## Method for getting each amino acid position and using it to make a separate residue selection, and calculate and write SASA to output. (ALL RESIDUES) ##
###########################################################################################################################################################
def find_Allchain(seq, chain, start):

    print("Start of chain is position " + str(start) + " and this sequence of length " + str(len(seq)) + " is:\n" + seq + "\n\n")

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
            cmd.select("sele, resn " + AA_attributes[ltr][0] + " and resi " + count + " and sidechain and chain " + keyVal)
            sasa = cmd.get_area("sele", 1, 0)            
            rel_sasa = sasa / AA_attributes[ltr][1]  # Retrieves max SASA attribute from the tuple value at [ltr]
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
                    
        # at the end of each calculation, clear the current selection
        cmd.delete("sele")
    return

##################################################################################################################################################################
## Method for getting each lysine position and using it to make a separate residue selection, and calculate and write SASA to output. (ONLY SPECIFIED RESIDUES) ##
##################################################################################################################################################################
def find_res(seq, chain, ch, start):
    for count, ltr in enumerate(seq, start):    
        if ltr == ch:
            ## Typecasting count (index or aa position) to a string allows the script to use count in a PyMOL selection-expression.
            count = str(count)
            residue = "" + ltr + count
            
            # Method call: Qualify this selection as valid; then check if electron density is actually present in the structure model
            presence = occupancy(count)
            occu_vals.clear()

            # Calculations are not performed if the electron density is missing from the model; instead, their values are
            # defaulted to "Not present in structure model", "N/A", and "N/A"
            if presence == True:
                # using the reference residue side chain, calculate relative SASA by accessing from the dictionary of max SASA values.
                # Chain must be specified every time, otherwise the value will be multiplied by the number of identical chains.
                cmd.select("sele, resn " + AA_attributes[ltr][0] + " and resi " + count + " and sidechain and chain " + keyVal)
                sasa = cmd.get_area("sele", 1, 0)                
                rel_sasa = sasa / AA_attributes[ltr][1]
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
            
            # at the end of each calculation, clear the current selection
            cmd.delete("sele")
    return
    
################################################################################################################################################
##~~~~~~~~~~~~~~~~ ___DRIVER CODE___~~~~~~~~~~~~~~~~##  (Writes to the csv file with each method call).

## "Stopwatch" starts now
start_time = time.time()

header = ["SOLVENT ACCESSIBLE SURFACE AREAS OF TARGET PROTEOME"]

with open('SASA_' + query + '_' + depth + '.csv', 'w', newline = '') as file: # write output values into csv row by row, vals in separate columns
    # create a .csv file writer object                                                                         
    writer = csv.writer(file, delimiter = ',')
    
    # Begin writing into the csv output file with a master header describing the job
    writer.writerow(header)

    # start with a fresh pymol session
    cmd.reinitialize

    # SASA settings
    cmd.set('dot_solvent', 1)  ## 1 is for solvent surface area. 0 is for total molecular surface area [default]
    cmd.set('dot_density', 4)  ## 1-4; defines quality (accuracy) of the calculation, better=more CPU

    # Import structure, then remove unwanted (non-amino acid, or "het") objects
    cmd.fetch(query)
    cmd.remove("het")
    
    # detect all chains present in the file and get the full fasta string sequence for each unique chain; remove the first line (unwanted header), then the whitespace,
    # leaving only the AA single-letter characters in the string.
    chainID_list = []    # for recording chain ID's
    chainPlusFasta = {}   # dictionary for matching the current retrieved chain ID with its associated fasta

    for chain in cmd.get_chains():
        fasta = cmd.get_fastastr("Chain " + str(chain)) # gets the associated fasta with each chain ID     
        fasta = fasta.split("\n",1)[1]   # removes the first line from the /n delimited fasta string
        fasta = re.sub("\n", "", fasta)  # removes remaining /n
        chainID_list.append(chain)       # collects chain IDs

        # Build a dictionary containing the chain ID and its associated (unique) fasta sequence. When the current value of 'fasta' is unique,
        # it is added as the value and its associated chain ID (value of 'chain') is its key.
        if fasta not in chainPlusFasta.values():
            chainPlusFasta.update({chain:fasta})

    # print(chainPlusFasta)    #NOTE this print() is for testing

    # Begin writing into the csv output file. For a multichain protein, this writes each chain to the same file separated by subheaders.
    chain_count = 0

    for keyVal in chainPlusFasta:
        chain_count += 1
        subheader = ["Residue", "Absolute SASA", "Relative SASA", "Exposed or Buried?", "PDB ID: " + query, "Chain " + str(chain_count) + " FASTA:", chainPlusFasta[keyVal]]
        writer.writerow("")
        writer.writerow(subheader)   # subheaders are labeled by chain number with each iteration

        # For each unique chain, iterate the residue positions into a Set named stored_residues, move the set into a List (so elements can be accessed by index) named residues, 
        # and then use List[0] to set enumerator's startpoint for recording residue positions
        cmd.iterate("chain " + str(keyVal), 'stored_residues.add(resv)')
        for i in stored_residues:
            residues.append(i)

        # start should take the smallest value from the list, not the one in [0] because that could vary if set() elements are added randomly... unless they "store" the order they are iterated... (?)
        # in which case, [0] is always the true start value because iterate() works in order?  
        start = residues[0]
        
        # Run the job with the query (requested PDB ID) based on its type. If a single amino acid is requested, the input is compared to an amino acid dictionary
        # to ensure it's a real amino acid; otherwise the job is skipped and the user is given an error message.
        if len(depth) == 1:
            if depth in AA_attributes.keys():
                find_res(chainPlusFasta[keyVal], keyVal, depth, start)
            else:
                print("\n\n####################################################################################\n# ATTN user! Check your batch file: '" + depth + "' is not a valid single-letter residue code. #\n####################################################################################\n\n")
                                                                                                                        
        elif depth == "ALL":
            find_Allchain(chainPlusFasta[keyVal], keyVal, start)
        
        # clear container objects for the next chain (if applicable) in this protein
        stored_residues.clear()
        residues.clear()

    ## "Stopwatch" stops now; print runtime
    stop_time = time.time()
    print("\nTime (seconds) taken for SASA calculations: " + str(stop_time - start_time) + "\n\nOutput file: " + "SASA_" + query + "_" + depth + ".csv was saved in the working directory: " + os.getcwd())

    ## End