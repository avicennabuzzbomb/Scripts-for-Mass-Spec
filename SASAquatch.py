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
including the pre-compiled bundles provided by Schrodinger.

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

print("Number of arguments entered when running this instance of `SASAquatch.py`:", len(sys.argv)-1) # check number of arguments entered
print("WARNING 10/20/2022: cmd.load is now supported but on Condor those files need to be available, explicitly listed in the batch list file and then sent to the cluster.\n This means `BulkSubmit.sh needs to be updated to check for the load condition and to then copy the listed files over to the cluster with everything else.")

# NOTE This script must take arguments at [3] and [4] because the intrepreter argument "pymol", the "-c" flag, and the script name are being interpred by this script as initial arguments [0], [1], and [2].
# QUERY: assign PDB ID to the first argument in the list of possible arguments, uppercase it (arg is case insensitive), and typecast to string
#query = str(sys.argv[3].upper())           # TODO is the upper() method call necessary? Forces all queries to be case insensitive BUT if a file is being identified for cmd.load,
                                            # TODO that is likely to be case-sensitive and may be causing the pdb file failures on the CHTC even though it works locally

query = str(sys.argv[3])                    # TODO testing query without uppercasing                                        

# DEPTH: assign query type ('ALL' as default, unless a 4th argument is explicitly entered: an amino acid single letter code or 3-letter code).
# checking first for another argument must be done before attempting to assign it; else a nonexistent argument will make the script close instead.
if len(sys.argv)-1 == 4:
    depth = str(sys.argv[4].upper())        # TODO is the upper() method call necessary? Forces all queries to be case insensitive BUT if a file is being identified for cmd.load,
                                            # TODO that is likely to be case-sensitive and may be causing the pdb file failures on the CHTC even though it works locally

else:
    depth = "ALL"

# check whether user is fetching a protein by its PDB ID or if the user wants to load an existing molecule file instead.
if "." in query:   #NOTE greedy for all file extensions - refine later
    mode = "load"
else:
    mode = "fetch"
    query = query.upper()

# finally, check whether user wants to include -het- chains in the structure (matters if SASA modeling is expected to be affected by presence of a biological ligand)
# TODO build another conditional - if het == someLigandCSVstring, import string into a list by default and call a remove-het helper method.
# TODO also requires converting the cmd.remove(het) codeline into a call to a new helper method for parsing a user-requested list of hets worth keeping.
# TODO example hets worth keeping: AMPPNP or AMPPCP in AHA2 or SERCA

#_____________________________________________________________________________________________________________________________________________________________________________________________________________
# FIXME FIXME FIXME Include a new method for calculating each chain in the context of the entire crystal unit (covers true oligomeric structures and also artifactual ones)
# FIXME FIXME FIXME This is a memory-intensive process because to get "relative" SASA values it would need to calculate each chain's per-res SASA in the absence of the other chains,
# FIXME FIXME FIXME ...AND to then re-do the calculation on each chain in the presence of the other chains. To minimize memory, the second step has to happen first; then, each chain would need to 
# FIXME FIXME FIXME be calculated on its own. Those values are then used as the "reference" for calculating the relative SASA of each chain.
# FIXME FIXME FIXME Also, using a "threshold" for relative SASA is sketchy. It may be better to use the actual value of water's SASA as the hard cutoff for solvent-accessibility.
#______________________________________________________________________________________________________________________________________________________________________________________________________________

## Dictionary of residue attributes needed for SASA calculations. 1 entry per amino acid; the single-letter AA code is the key,
## and indices [0] - [3] are, respectively, the 3-letter code for the amino acid, the total sidechain (R-group) SASA for the amino acid when free
## in solution, the total SASA for the whole amino acid in solution, and then the full amino acid's name.

#NOTE this Dictionary is used by all query types
## Dict. Name = { 3-letter code (key): ( [0, single-letter code], [1, maxSASA], [2, sideSASA], [3, full name] ) }
AA_attributes = {'ARG' : ('R', 249.4673309, 380.2754517, 'arginine'),
                 'HIS' : ('H', 194.6491394, 323.2107239, 'histidine'),
                 'LYS' : ('K', 214.5114136, 344.1722412, 'lysine'),
                 'ASP' : ('D', 142.6577301, 272.8369141, 'aspartic acid'),  
                 'GLU' : ('E', 173.6901855, 303.7915955, 'glutamic acid'),
                 'SER' : ('S', 110.3298187, 240.5674744, 'serine'),
                 'THR' : ('T', 141.3123779, 268.1856995, 'threonine'),
                 'ASN' : ('N', 155.2926178, 283.9153137, 'aspargine'),
                 'GLN' : ('Q', 184.1813354, 313.8310242, 'glutamine'),   
                 'CYS' : ('C', 131.9754333, 261.8157043, 'cysteine'),   
                 'GLY' : ('G', 34.49541855, 200.8881836, 'glycine'),   
                 'PRO' : ('P', 153.1377106, 258.9144592, 'proline'),   
                 'ALA' : ('A', 91.39417267, 229.5201416, 'alanine'),   
                 'VAL' : ('V', 157.4189758, 279.6635132, 'valine'),   
                 'ILE' : ('I', 183.812027, 306.1242676, 'isoleucine'),    
                 'LEU' : ('L', 181.1264496, 308.8426819, 'leucine'),   
                 'MET' : ('M', 192.694458, 323.5989075, 'methionine'),    
                 'PHE' : ('F', 218.0111237, 342.800415, 'phenylalanine'),   
                 'TYR' : ('Y', 231.1389923, 357.9006042, 'tyrosine'),   
                 'TRP' : ('W', 262.1514587, 387.6826172, 'tryptophan')}

# NOTE this Dictionary is used only by non-default queries
AA_letterCode = {'R' : ('ARG'), 'H' : ('HIS'), 'K' : ('LYS'), 'D' : ('ASP'), 'E' : ('GLU'), 'S' : ('SER'), 'T' : ('THR'), 'N' : ('ASN'), 'Q' : ('GLN'), 'C' : ('CYS'),   
                 'G' : ('GLY'), 'P' : ('PRO'), 'A' : ('ALA'), 'V' : ('VAL'), 'I' : ('ILE'), 'L' : ('LEU'), 'M' : ('MET'), 'F' : ('PHE'), 'Y' : ('TYR'), 'W' : ('TRP'),
                 'ARG' : ('ARG'), 'HIS' : ('HIS'), 'LYS' : ('LYS'), 'ASP' : ('ASP'), 'GLU' : ('GLU'), 'SER' : ('SER'), 'THR' : ('THR'), 'ASN' : ('ASN'), 'GLN' : ('GLN'), 'CYS' : ('CYS'),   
                 'GLY' : ('GLY'), 'PRO' : ('PRO'), 'ALA' : ('ALA'), 'VAL' : ('VAL'), 'ILE' : ('ILE'), 'LEU' : ('LEU'), 'MET' : ('MET'), 'PHE' : ('PHE'), 'TYR' : ('TYR'), 'TRP' : ('TRP'),} 
                 
# The cutoff value for relative SASA. Residues with a value greater than 0.25 are considered solvent-exposed, otherwise are considered buried.
# FIXME/NOTE Note that this is an *arbitrary* cutoff, generally accepted. There are arguments to be made against this value.
threshold = 0.25

## Initialize empty variables and data stuctures for handling pymol's iterator stream
stored_residues = set() # for capturing iterate's `resi` list
atmNum = 0
resname = set()         # for capturing iterate's `resn` list, called by the current `resi`, and clear duplicates; and a string variable to capture the `resn`

#############################################################################################################################################################
## Helper Method: evaluate if electron density exists at the current residue position is actually present in the structure model, and annotate accordingly ##
#############################################################################################################################################################
def occupancy(position, chain):

    #NOTE Currently, pymol treats invalid selections are treated as real objects, and default to 0 when calculations are performed. This leads to falsely labeling residues
    # as buried and their SASA as 0, when that may not reflect reality. May need to extend this method to check for bad selections.

    # By default, pymol reports occupancy as '1' when it can't find electron density. Check for absent density so it can be flagged in the output file.
    atmNum = cmd.count_atoms("resi " + position + " and chain " + chain)

    if atmNum == 1:   
        return False
    else:
        return True

#################################################################################################
## Helper Method: get and sort all values of `resi` so that SASA calculations can be performed ##
#################################################################################################
def extractResCode(selexpression, stored_residues):

    resPositions = []     # For storing the typecasted instance of the set `stored_residues` during each chain's calculations

    print("Gathering `resi` values with the selection expression ", selexpression)

    # For each unique chain, iterate the residue positions into a Set named `stored_residues`
    cmd.iterate(selexpression, 'stored_residues.add(resi)')
    resPositions = list(stored_residues)                 # re-casting the set to a list, then back to a set erases duplicates, because sets do not keep duplicate elements
    resPositions.sort()                                  # sort() first to apply numerical character sorting
    resPositions.sort(key=len)                           # sort(key=len) uses the string element's length as the key to sort against, so now both sort rules apply
    
    # Resets `stored_residues` list object to an empty set
    stored_residues.clear()
    print("\nContents of set `stored_residues` should be empty and are now:", stored_residues,"\nContents of `resPositions` are now:", resPositions)

    return resPositions

#######################################################################################################################
## Helper Method: call the current value of `resn` so that SASA calculations can reference the amino acid Dictionary ##
#######################################################################################################################
def getRESN(resi, chain):    # `resi` is the current position (from `currposition`)

    # use the current `resi` to get the corresponding three-letter code, `resn`      
    cmd.iterate("resi " + resi + " and chain " + chain, 'resname.add(resn)')
    print("`resname` contents:",resname)

    # pop the current `resn` from the set into `currRes`; this also empties the set, making it ready for the next pass in this loop.
    currRes = resname.pop()

    return currRes

##################################################################################################################################################################
#|  SASA METHOD: Uses a List `stored_residues` populated with `resi` values for all selection-expressions; `resi` retrieves the PSE residue position as shown.  |#
##################################################################################################################################################################
def find_Allchain_resi(seq, chain, resi, writer):   # seq is a string type, chain is a character type, resi is a List-type

    print("Start of chain is position " + resi[0] + " and this sequence of length " + str(len(seq)) + " is:\n" + seq + "\n\n")

    currRes = ""
    currposition = ""  
    residue = ""

    for i in range(len(resi)):
        # store the current `resi` and get its corresponding residue code, then call `occupancy()`: Confirm that this selection is actually present in the structure model
        currposition = resi[i]        

        # Method calls to check for the presence of electron density in the current selection, and then to collect the `resn` attached to each position.
        presence = occupancy(currposition, chain)
        currRes = getRESN(currposition, chain)

        residue = "" + AA_attributes[currRes][0] + currposition
        
        ## Calculations are performed if the density of the current residue can be explicity selected within the Pymol session.
        if presence == True:           
            
            # first, select the whole amino acid (current residue) and record its SASA
            tot_sasa = cmd.get_area("resi " + currposition + " and chain " + chain)

            # next, clear the current selector, re-select the current residue, and record the SASA of its sidechain
            side_sasa = cmd.get_area("resi " + currposition + " and chain " + chain + " and sidechain")

            # now, calculate the relative SASA of the sidechain and the whole residue separately, and assess "burial" status based on the current threshold.
            # Burial status defaults to "exposed" UNLESS the value falls below the threshold. Relative SASA values are calculated using the maximum SASA
            # values stored in the Dictionary named `AA_attributes`.
            totrel_sasa = tot_sasa / AA_attributes[currRes][2]    # relative SASA of the whole residue; divide calculated SASA by max SASA of naked residue
            totburial = "Exposed"
            siderel_sasa = side_sasa / AA_attributes[currRes][1]  # relative SASA of the sidechain only; divide calculared SASA by max SASA of the sidechain of the nake residue
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

    return    # DONE

####################################################################################################################################################################
#|  WRITER METHOD: Uses a List `stored_residues` populated with `resi` values for all selection-expressions; `resi` retrieves the PSE residue position as shown.  |#
####################################################################################################################################################################
def GO(query, header, requested, selexpression, stored_residues, mode, depth="ALL"):

    ## Job description
    print("Query: ", query,"\nResidue(s) requested:", depth)

    ## "Stopwatch" starts now
    start_time = time.time()

    with open('SASA_' + query + '_' + depth + '.csv', 'w', newline = '') as file: # write output values into csv row by row, vals in separate columns
        # create a .csv file writer object                                                                         
        writer = csv.writer(file, delimiter = ',')
    
        # Begin writing into the csv output file with a master header describing the job
        writer.writerow(header)

        # start with a fresh pymol session
        cmd.reinitialize()

        # SASA settings
        cmd.set('dot_solvent', 1)  ## 1 is for solvent surface area. 0 is for total molecular surface area [default]
        cmd.set('dot_density', 4)  ## 1-4; defines quality (accuracy) of the calculation, better=more CPU

        # Import structure, then remove unwanted (non-amino acid, or "het") objects. #FIXME this may remove biochemically important groups from proteins (ex, the `CHO` fluorophore in GFP; see 2B3P, resi #65-67 for details)
        if mode == "fetch":
            cmd.fetch(query)
            print("Downloaded query ID: " + query)
        else:
            cmd.load(query)
            print("Loaded query file: " + query)

        # TODO convert this into a remove-het helper method for parsing user-requested hets that should remain in the structure while SASA analysis is happening.
        # TODO example hets worth keeping: AMPPNP or AMPPCP in AHA2 or SERCA
        cmd.remove("het")
    
        # detect all chains present in the file and get the full fasta string sequence for each unique chain; remove the first line (unwanted header), then the whitespace,
        # leaving only the AA single-letter characters in the string.
        chainID_list = []     # List for recording chain ID's
        chainPlusFasta = {}   # Dictionary for matching the current retrieved chain ID with its associated fasta

        # Build a list of all chains in the structure file
        for chain in cmd.get_chains():
            fasta = cmd.get_fastastr("Chain " + str(chain)) # gets the associated fasta with each chain ID     
            fasta = fasta.split("\n",1)[1]   # removes the first line from the /n delimited fasta string
            fasta = re.sub("\n", "", fasta)  # removes remaining /n
            chainID_list.append(chain)       # Records chain IDs. TODO Important to have this in addition to `chainPlusFasta{}` for when all chains present are called. 

            # Build a dictionary containing the chain ID and its associated (unique) fasta sequence. When the current value of 'fasta' is unique,
            # it is added as the value and its associated chain ID (value of 'chain') is its key.
            # NOTE that this is not always the correct way to look at surface area - a biological unit may actually be multimeric vs. just being a crystallization artifact.
            # TODO need to implement more code to deal with this? Should keep all the chains by default, leave user to specify that unique chains are desired.

            # NOTE FIXME NOTE FIXME a simple implementation would be to 1) add another parameter that this method accepts, with a default value if not specified
            # and 2) to have this value passed here, and checked in an if-then branch that builds in the `if fasta not in ...` statement. This branch is taken if unique chains
            # are desired. 3) Finally, the alternative branch is all chains that are present. 4) Yet another branch may specify that chains are done either in isolation,
            # or in complex. The value in this type of calculation would be in comparing protein components in a complex before and after they adopt a multimeric or complex configuration.
            # NOTE FIXME NOTE FIXME if additional conditions are in consideration, switch statements may be preferable to/faster than multiple if-then-else branches.

            if fasta not in chainPlusFasta.values():
                chainPlusFasta.update({chain:fasta})

        # Begin writing into the csv output file. For a multichain protein containing unique chains, this writes each chain to the same file separated by subheaders.
        for keyVal in chainPlusFasta:

            # Write the subheader. Subheaders are labeled by chain number and updated with each iteration
            writer.writerow("")
            subheader = ["Residue", "Total SASA", "Total Relative SASA", "Total: Exposed or Buried?", "Sidechain SASA", "Sidechain Relative SASA", "Sidechain: Exposed or Buried?", "PDB ID: " + query, "Chain " + str(keyVal) + " FASTA:", chainPlusFasta[keyVal]]
            writer.writerow(subheader)

            # update the selection-expression string for iterating `resi`
            selexpression = "" + "chain " + str(keyVal) + requested
            print("Requesting next set of `resi` with selexpression =`" + selexpression + "`")  ## DEBUG
            print("Set `stored_residues` should be empty for this iteration before method call; `stored_residues` contains:",stored_residues)  ## DEBUG

            # Method call to `extractResCode()`: generate the List of `resi` based on `selexpression`'s value
            resPositions = extractResCode(selexpression, stored_residues)

            # Method call to `find_Allchain_resi()`: get and print the SASA values for each requested residue.
            find_Allchain_resi(chainPlusFasta[keyVal], keyVal, resPositions, writer)  # the sorted list of indices `resPositions[]` is passed to the counting method

    ## "Stopwatch" stops now; print runtime
    stop_time = time.time()
    print("\nTime (seconds) taken for SASA calculations: " + str(stop_time - start_time) + "\n\nOutput file: " + "SASA_" + query + "_" + depth + ".csv was saved in the working directory: " + os.getcwd())

    return    # DONE

#######################################################################################################################################################
##~~~~~~~~~~~~~~~~ ___DRIVER CODE (VROOM VROOM!)___~~~~~~~~~~~~~~~~##  (MAIN method; calls all others. Writes to the csv file with each method call).

# BEGIN ~ ~ ~
header = ["SOLVENT ACCESSIBLE SURFACE AREAS OF TARGET PROTEOME"]
requested = "" ## NOTE the `requested` variable stores the `selexpression` string's extension which affects which `resi` are iterated.
selexpression = ""

# check user's requested job parameters
if depth == "ALL":
    # call GO()
    GO(query, header, requested, selexpression, stored_residues, mode)    # NOTE: an empty string extends `selexpression` when set to default: 'ALL'

elif depth in AA_letterCode.keys():
    # call GO()
    # TODO convert this case into a string-splitting process, where "," is the delimiter. This would avoid requiring a separate data structure by letting me stick with strings.
    
    # TODO this will need to be able to accept a string-as-list, where "," separates elements and each are bundled into a list type (even if only one code is requested)
    # TODO then, the `requested=` needs to be modified to reflect the number of items requested when the number of aa codes > 1.
    depth = AA_letterCode[depth]        ## NOTE currently this only works with 1 requested lettercode at a time; converting depth into a list that accepts all valid, unique instances of an amino acid and then getting their resi values (and sorting in order) will allow multiple aa requests in the same run.
    requested = " and resn " + depth                                # NOTE: if `sys.argv[4]` is found to be valid, the request is run with the appropriate 3-letter code.
    GO(query, header, requested, selexpression, stored_residues, mode, depth)

else:
    # error - user must enter appropriate <args>
    err1 = "#   ATTN user! Check your batch file: `" + depth + "` is not a valid single or 3-letter residue code.   #"
    err2 = "#" * len(err1)
    print("\n\n",err2,"\n\n",err1,"\n\n",err2,"\n\n")

## END ~ ~ ~