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

Example from a shell (to run locally):

shell> pymol -c script.py
With arguments (sys.argv becomes ["script.py", "foo", "bar"]):

shell> pymol -c script.py foo bar
(linux)> ./pymol -c script.py foo bar

Example from a running PyMOL instance:
PyMOL> run script.py

## NOTE - To run a script from the PyMOL command line, it needs to be saved in Pymol's working directory.
## NOTE - PyMOL's command line commands are Unix-like (same as Bash and Windows command line)
## NOTE - 'reinitialize' removes all objects from this pymol session's memory (it's a reset as if the program was restarted) 
## NOTE - the -c flag (or -cp) means CLASSPATH - it tells the shell that it needs to use the pymol 'class' to run
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

# FIXME FIXME FIXME Include a new method for calculating each chain in the context of the entire crystal unit (covers true oligomeric structures and also artifactual ones)
# FIXME FIXME FIXME This is a memory-intensive process because to get "relative" SASA values it would need to calculate each chain's per-res SASA in the absence of the other chains,
# FIXME FIXME FIXME ...AND to then re-do the calculation on each chain in the presence of the other chains. To minimize memory, the second step has to happen first; then, each chain would need to 
# FIXME FIXME FIXME be calculated on its own. Those values are then used as the "reference" for calculating the relative SASA of each chain.

# FIXME FIXME FIXME Also, using a "threshold" for relative SASA is sketchy. It may be better to use the actual value of water's SASA as the hard cutoff for solvent-accessibility.

# for testing, until a list (whether in text file or directly in submit file) of PDB IDs is ready
depth = "ALL"

print("query: " + query)
print("depth: " + depth)

## Dictionary of residue attributes needed for SASA calculations. 1 entry per amino acid; the single-letter AA code is the key,
## and indices [0] - [3] are, respectively, the 3-letter code for the amino acid, the total sidechain (R-group) SASA for the amino acid when free
## in solution, the total SASA for the whole amino acid in solution, and then the full amino acid's name.
## Dict. Name = { key: ( [0], [1], [2], [3] ) } 
AA_attributes = {'R' : ('arg', 249.4673309, 380.2754517, 'arginine'),
                 'H' : ('his', 194.6491394, 323.2107239, 'histidine'),
                 'K' : ('lys', 214.5114136, 344.1722412, 'lysine'),
                 'D' : ('asp', 142.6577301, 272.8369141, 'aspartic acid'),  
                 'E' : ('glu', 173.6901855, 303.7915955, 'glutamic acid'),
                 'S' : ('ser', 110.3298187, 240.5674744, 'serine'),
                 'T' : ('thr', 141.3123779, 268.1856995, 'threonine'),
                 'N' : ('asn', 155.2926178, 283.9153137, 'aspargine'),
                 'Q' : ('gln', 184.1813354, 313.8310242, 'glutamine'),   
                 'C' : ('cys', 131.9754333, 261.8157043, 'cysteine'),   
                 'G' : ('gly', 34.49541855, 200.8881836, 'glycine'),   
                 'P' : ('pro', 153.1377106, 258.9144592, 'proline'),   
                 'A' : ('ala', 91.39417267, 229.5201416, 'alanine'),   
                 'V' : ('val', 157.4189758, 279.6635132, 'valine'),   
                 'I' : ('ile', 183.812027, 306.1242676, 'isoleucine'),    
                 'L' : ('leu', 181.1264496, 308.8426819, 'leucine'),   
                 'M' : ('met', 192.694458, 323.5989075, 'methionine'),    
                 'F' : ('phe', 218.0111237, 342.800415, 'phenylalanine'),   
                 'Y' : ('tyr', 231.1389923, 357.9006042, 'tyrosine'),   
                 'W' : ('trp', 262.1514587, 387.6826172, 'tryptophan')}   
'''
nucleic_attributes = {' GUA ' : ('guanine', 0),      # TODO
                      ' CYT ' : ('cytosine', 0),
                      ' ADE ' : ('adenine', 0),
                      ' THY ' : ('thymine', 0)}     # TODO: complete this dictionary containing nucleic acid names and attributes; 0 is currently a placeholder for the actual SASA value
                                                     # of the nitrogenous base 'sidechain'; need to confirm exactly how this string appears
                                                    
small_molecules = { ' AlF4 ' : ('aluminum tetrafluoride', 0),    # TODO
                    ' PMSF ' : ('phenylmethylsulfonyl fluoride', 0)}   # TODO: make a dictionary of non-nucleic acid small molecule co-crystallants that are in-chain with fasta; 
                                                                # or, make a way to count them as "not aa and not nucleic". Current entries are just examples that may exist.
'''

# The cutoff value for relative SASA. Residues with a value greater than 0.25 are considered solvent-exposed, otherwise are considered buried.
threshold = 0.25

# Initialize empty data stuctures for handling pymol's iterator stream and generally crappy API
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
    if atmNum == 1:   # By default, pymol reports occupancy as '1' when it can't find electron density
        return False
    else:
        return True
###########################################################################################################################################################
## Method for getting each amino acid position and using it to make a separate residue selection, and calculate and write SASA to output. (ALL RESIDUES) ##
###########################################################################################################################################################
def find_Allchain(seq, chain, start):
    pattern = ""
    print("Start of chain is position " + str(start) + " and this sequence of length " + str(len(seq)) + " is:\n" + seq + "\n\n")
    '''
    if "nucleic_" in seq:
        pattern = "nucleic_"
        seq = [re.sub(pattern,'',i) for i in seq]  # trim off the annotation; then check for all instances of nucleic acids and treat them as units
        # TODO implement code to "enumerate" with substrings (nucleic acid codes)
    '''
    for count, ltr in enumerate(seq, start): 
        ## Typecasting the variable `count` (index, or aa position) to a string type allows the script to use `count` in a PyMOL selection-expression.
        count = str(count)
        residue = "" + ltr + count
            
        # Method call: Confirm that this selection is actually present in the structure model
        presence = occupancy(count)

        ## Calculations are performed if the density of the current residue can be explicity selected within the Pymol session.
        if presence == True:
            # Calculate the total SASA of the current selected residue (sidechain and backbone), and then use this to calculate relative SASA based on the dictionary values.
            # Chain identifier must be specified each time to avoid multiple selections occurring in identical chains.            
            
            # first, select the whole amino acid (current residue) and record its SASA
            cmd.select("sele, resn " + AA_attributes[ltr][0] + " and resi " + count + " and chain " + keyVal)
            tot_sasa = cmd.get_area("sele")

            # next, clear the current selector, re-select the current residue, and record the SASA of its sidechain
            cmd.delete("sele")
            cmd.select("sele, resn " + AA_attributes[ltr][0] + " and resi " + count + " and sidechain and chain " + keyVal)
            side_sasa = cmd.get_area("sele")

            # now, calculate the relative SASA of the sidechain and the whole residue separately, and assess "burial" status based on the current threshold.
            # Burial status defaults to "exposed" UNLESS the value falls below the threshold. Relative SASA values are calculated using the maximum SASA
            # values stored in the Dictionary named `AA_attributes`.
            totrel_sasa = tot_sasa / AA_attributes[ltr][2]
            totburial = "Exposed"
            siderel_sasa = side_sasa / AA_attributes[ltr][1]
            sideburial = "Exposed"

            if totrel_sasa <= threshold:
                totburial = "Buried"

            if siderel_sasa <= threshold:
                sideburial = "Buried"
            # before printing, typecast all float values to string types
            tot_sasa = str(tot_sasa)
            totrel_sasa = str(totrel_sasa)
            side_sasa = str(side_sasa)
            siderel_sasa = str(siderel_sasa)
        ## When a residue's density cannot be selected in the Pymol session, all values default to "Not present..." and "N/A"    
        else:
            side_sasa = tot_sasa = "Not present in structure model"
            siderel_sasa = totrel_sasa = "N/A"
            sideburial = totburial = "N/A"
        
        ## Each SASA printed to console and then to output
        print(residue + " | " + tot_sasa + " | " + totrel_sasa + " | " + totburial + " | " + side_sasa + " | " + siderel_sasa + " | " + sideburial)
        current_row = [residue, tot_sasa, totrel_sasa, totburial, side_sasa, siderel_sasa, sideburial]
        writer.writerow(current_row)            
                    
        # at the end of each calculation, clear the current selection from the pymol session (saves memory, allows calculations to be per-residue only)
        cmd.delete("sele")
    return

##################################################################################################################################################################
## Method for getting each lysine position and using it to make a separate residue selection, and calculate and write SASA to output. (ONLY SPECIFIED RESIDUES) ##
##################################################################################################################################################################
def find_res(seq, chain, ch, start):    ### FIXME FIXME FIXME This will do better as a subsidiary or helper method to find_Allchain(). Ie,. in this method create a mini Dictionary containing only the desired residue positions, and then cycle through those using find_Allchain()
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
                cmd.select("sele, resn " + AA_attributes[ltr][0] + " and resi " + count + " and chain " + keyVal)
                sasa = cmd.get_area("sele")                
                rel_sasa = sasa / AA_attributes[ltr][2]
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
    chainID_list = []     # List for recording chain ID's
    chainPlusFasta = {}   # Dictionary for matching the current retrieved chain ID with its associated fasta

    # pattern = "a generic substring containing whitespace flanking ASCII characters - nucleic acids and small molecules are indicated this way in .cif files"

    chaincontent = ""       # a string variable that preppends to the current fasta. "" indicates amino acids only. If other elements are present in the chain,
                            # they are detected in the if branch below and the chain is modified to indicate
    # Using a loop, iterate through each chain, storing both its ID
    for chain in cmd.get_chains():
        fasta = cmd.get_fastastr("Chain " + str(chain)) # gets the associated fasta with each chain ID     
        fasta = fasta.split("\n",1)[1]   # removes the first line from the /n delimited fasta string
        fasta = re.sub("\n", "", fasta)  # removes remaining /n
        '''    
        if pattern in fasta:             # TODO make sure this checks fasta type (amino acid vs. nucleic acid) appropriately and annotates the string
            string = "this is pseudocode: check nucleic acids dictionary, then small molecules dictionary (impl. later) for a match"                             # checks that the fasta chain A has logical values (avoid selector errors)
            # fasta is annotated with a notation that is detectable in method calls
            if pattern in nucleic_attributes:
                chaincontent = "nucleic_"
                fasta = chaincontent + fasta                             # NOTE additional fasta processing is required here; currently this is blind to terms like "?" and " UNK " and nucleic acids.
            if pattern in small_molecules:   # or, not a nucleic or amino acid
                chaincontent = "smallMolecule_"                             # fasta is collected in aa chain oriented fashion, so it should be possible to treat each chain uniquely
                fasta = chaincontent + fasta                         # (ie., make another dictionary of nested lists for nucleic acids that can be called when a nucleic acid term is recognized)
        '''
        chainID_list.append(chain)       # collects chain IDs  TODO double-check, is this list necessary? (maybe not for calling `find_Allchain()`, but maybe necessary for `find_res()`)

        # Build a dictionary containing the chain ID and its associated (unique) fasta sequence. When the current value of 'fasta' is unique,
        # it is added as the value and its associated chain ID (value of 'chain') is its key. NOTE that this is not always the correct way to 
        # look at surface area - a biological unit may actually be multimeric vs. just being a crystallization artifact. TODO implements more code
        # to deal with this?
        if fasta not in chainPlusFasta.values():
            chainPlusFasta.update({chain:fasta})

    # Begin writing into the csv output file. For a multichain protein containing unique chains, this writes each chain to the same file separated by subheaders.
    chain_count = 0

    for keyVal in chainPlusFasta:
        chain_count += 1                ## TODO create a dictionary to handle alphanumeric conversion here (ie, Chain 1 = Chain A, Chain 2 = Chain B, etc.)
        subheader = ["Residue", "Total SASA", "Total Relative SASA", "Total: Exposed or Buried?", "Sidechain SASA", "Sidechain Relative SASA", "Sidechain: Exposed or Buried?", "PDB ID: " + query, "Chain " + str(chain_count) + " FASTA:", chainPlusFasta[keyVal]]
        writer.writerow(subheader)   # subheaders are labeled by chain number with each iteration

        # FIXME Willing to bet that sets are not necessary here. It may be possible to pre-filter the iterator's resv parameter to only keep the minimum value, because
        # FIXME that's really the only reason we would need to get the residue values.
        # For each unique chain, iterate the residue positions into a Set named `stored_residues`
        cmd.iterate("chain " + str(keyVal), 'stored_residues.add(resv)')

        # The starting index of the fasta must be the minimum value in the set; this is defined as `start` point for calculating SASA.
        # (This line of code explicitly accounts for negative resi values).
        start = min(stored_residues)
        
        # Run the job with the query (requested PDB ID) based on its type. If a single amino acid is requested the value of depth will be the character representing that amino acid.
        # In that case the input is compared to an amino acid dictionary to ensure it's a real amino acid (user control); otherwise the job is skipped and the user is given an error message.
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
