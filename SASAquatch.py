## Calculate the Solvent-Accessible Surface Area (SASA) of individual residues at a time or as a group in PyMOL,
## then print to an output file.

## NOTE This script works by starting a fresh PyMOL session and then loading each requested pymol session (or fetching each requested PDB ID). 
## To access individual residues one at a time, a string is made from PyMOL's fasta sequence command, and this string is copied and "cleaned" 
## so that a version containing only the amino acids is left. Then, the string is used to get the position (index) of each residue from python's enumerator.
## Whenever a residue is found in the string, its index is recorded where the first index of the string is set to '1'. Then, PyMOL selection and SASA calculation
## functions are called and written to a new csv output file. 

## NOTE This script can be run on a computing cluster for larger jobs

"""
The recommended way to run PyMOL-Python scripts is by using PyMOL as the interpreter. This is supported by all versions of PyMOL, 
including the pre-compiled bundles provided by Schrödinger.

Example from a shell:

shell> pymol -c script.py
With arguments (sys.argv becomes ["script.py", "foo", "bar"]):

shell> pymol -c script.py foo bar
(linux)> ./pymol -c script.py foo bar

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
#depth = str(sys.argv[4].upper()) # depth: assign query type ('ALL' or aa single letter code) to the second argument. 
                                 # NOTE This script must take arguments at [3] and [4] because the intrepreter argument "pymol", the "-c" flag, and the script name are being interpreted
                                 # by this script as initial arguments [0], [1], and [2]. Unclear why this is the case. TODO eventually fix that, but for now everything works!

# for testing, until a list (whether in text file or directly in submit file) of PDB IDs is ready
depth = "ALL"

print("query: " + query)
print("depth: " + depth)

## NOTE This Dictionary contains the maximum calculated SASA (dot_solvent = 1, dot_density = 4) value and 3-letter residue code as a tuple value (changes will only be coded directly)
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

nucleic_attributes = {' GUA ' : ('guanine', 0),      # TODO
                      ' CYT ' : ('cytosine', 0),
                      ' ADE ' : ('adenine', 0),
                      ' THY ' : ('thymine', 0)}     # TODO: complete this dictionary containing nucleic acid names and attributes; 0 is currently a placeholder for the actual SASA value
                                                     # of the nitrogenous base 'sidechain'; need to confirm exactly how this string appears
                                                    
small_molecules = { ' AlF4 ' : ('aluminum tetrafluoride', 0),    # TODO
                    ' PMSF ' : ('phenylmethylsulfonyl fluoride', 0)}   # TODO: make a dictionary of non-nucleic acid small molecule co-crystallants that are in-chain with fasta; 
                                                                # or, make a way to count them as "not aa and not nucleic". Current entries are just examples that may exist.


# The cutoff value for relative SASA. Residues with a value greater than 0.25 are considered solvent-exposed, otherwise are considered buried.
threshold = 0.25

# Initialize empty data stuctures for handling python's iterator stream and generally crappy API
stored_residues = set()
atmNum = 0

#############################################################################################################################################################
## Helper Method: evaluate if electron density exists at the current residue position is actually present in the structure model, and annotate accordingly ##
#############################################################################################################################################################
def occupancy(count):
    ##TODO implement an if branch to check for a valid selection; if invalid, print to error file and skip SASA calculation.
    ##TODO Each bad query will have its own file or entry in a file. Currently, invalid selections are treated as real objects, and default to 0
    ##TODO when calculations are performed. This leads to falsely labeling residues as buried and their SASA as 0, when that may not reflect reality.

    atmNum = cmd.count_atoms("resi " + count)
    if atmNum == 1:
        return False
    else:
        return True
###########################################################################################################################################################
## Method for getting each amino acid position and using it to make a separate residue selection, and calculate and write SASA to output. (ALL RESIDUES) ##
###########################################################################################################################################################
def find_Allchain(seq, chain, start):
    pattern = ""
    print("Start of chain is position " + str(start) + " and this sequence of length " + str(len(seq)) + " is:\n" + seq + "\n\n")
    if "nucleic_" in seq:
        pattern = "nucleic_"
        seq = [re.sub(pattern,'',i) for i in seq]  # trim off the annotation; then check for all instances of nucleic acids and treat them as units
        # TODO implement code to "enumerate" with substrings (nucleic acid codes)

    for count, ltr in enumerate(seq, start): 
        ## Typecasting count (index or aa position) to a string allows the script to use "count" in a PyMOL selection-expression.
        count = str(count)
        residue = "" + ltr + count
            
        # Method call: Qualify this selection is actually present in the structure model
        presence = occupancy(count)

        # Calculations are not performed if the electron density is missing from the model; instead, their values are
        # defaulted to "Not present in structure model", "N/A", and "N/A"
        if presence == True:
            # Calculate SASA of the current selected sidechain, and then use this to calculate relative SASA based on the dictionary values.
            # Chain identifier must be specified each time to avoid multiple selections occurring in identical chains.            
            cmd.select("sele, resn " + AA_attributes[ltr][0] + " and resi " + count + " and sidechain and chain " + keyVal)
            sasa = cmd.get_area("sele")            
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
    # TODO make it possible here to execute the script on pse or cif files that are saved in the current directory; currently the script attempts
    # to call the PDB with whatever first argument is supplied to it. This should check that the query is
    #query = "Full-length_AHA2 mixed xlinks_11-06-2019.pse" # NOTE temporary, remember to revert to fetching
    #cmd.fetch(query)
    cmd.load(query)  # NOTE temporary, remember to revert to fetching
    cmd.remove("het")
    
    # detect all chains present in the file and get the full fasta string sequence for each unique chain; remove the first line (unwanted header), then the whitespace,
    # leaving only the AA single-letter characters in the string.
    chainID_list = []     # List for recording chain ID's
    chainPlusFasta = {}   # Dictionary for matching the current retrieved chain ID with its associated fasta

    pattern = "a generic substring containing whitespace flanking ASCII characters - nucleic acids and small molecules are indicated this way in .cif files"

    chaincontent = ""       # a string variable that preppends to the current fasta. "" indicates amino acids only. If other elements are present in the chain,
                            # they are detected in the if branch below and the chain is modified to indicate

    for chain in cmd.get_chains():
        fasta = cmd.get_fastastr("Chain " + str(chain)) # gets the associated fasta with each chain ID     
        fasta = fasta.split("\n",1)[1]   # removes the first line from the /n delimited fasta string
        fasta = re.sub("\n", "", fasta)  # removes remaining /n
        if pattern in fasta:             # TODO make sure this checks fasta type (amino acid vs. nucleic acid) appropriately and annotates the string
            string = "this is pseudocode: check nucleic acids dictionary, then small molecules dictionary (impl. later) for a match"                             # checks that the fasta chain A has logical values (avoid selector errors)
            # fasta is annotated with a notation that is detectable in method calls
            if pattern in nucleic_attributes:
                chaincontent = "nucleic_"
                fasta = chaincontent + fasta                             # NOTE additional fasta processing is required here; currently this is blind to terms like "?" and " UNK " and nucleic acids.
            if pattern in small_molecules:   # or, not a nucleic or amino acid
                chaincontent = "smallMolecule_"                             # fasta is collected in aa chain oriented fashion, so it should be possible to treat each chain uniquely
                fasta = chaincontent + fasta                         # (ie., make another dictionary of nested lists for nucleic acids that can be called when a nucleic acid term is recognized)

        chainID_list.append(chain)       # collects chain IDs

        # Build a dictionary containing the chain ID and its associated (unique) fasta sequence. When the current value of 'fasta' is unique,
        # it is added as the value and its associated chain ID (value of 'chain') is its key. NOTE that this is not always the correct way to 
        # look at surface area - a biological unit may actually be multimeric vs. just being a crystallization artifact. TODO implements more code
        # to deal with this?
        if fasta not in chainPlusFasta.values():
            chainPlusFasta.update({chain:fasta})

    # Begin writing into the csv output file. For a multichain protein, this writes each chain to the same file separated by subheaders.
    chain_count = 0

    for keyVal in chainPlusFasta:
        chain_count += 1
        subheader = ["Residue", "Absolute SASA", "Relative SASA", "Exposed or Buried?", "PDB ID: " + query, "Chain " + str(chain_count) + " FASTA:", chainPlusFasta[keyVal]]
        writer.writerow("")
        writer.writerow(subheader)   # subheaders are labeled by chain number with each iteration

        # For each unique chain, iterate the residue positions into a Set named stored_residues
        cmd.iterate("chain " + str(keyVal), 'stored_residues.add(resv)')

        # The starting index of the fasta must be the minimum value in the set; this is defined as `start` point for calculating SASA. Accounts for negative resi values.
        start = min(stored_residues)
        
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

    ## "Stopwatch" stops now; print runtime
    stop_time = time.time()
    print("\nTime (seconds) taken for SASA calculations: " + str(stop_time - start_time) + "\n\nOutput file: " + "SASA_" + query + "_" + depth + ".csv was saved in the working directory: " + os.getcwd())

    ## End
