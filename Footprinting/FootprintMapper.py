## 

## NOTE This script works by starting a fresh PyMOL session and then loading each requested pymol session (or fetching each requested PDB ID). 
## To access individual residues one at a time, a string is made from PyMOL's fasta sequence command, and this string is copied and "cleaned" 
## so that a version containing only the amino acids is left. Then, the string is used to get the position (index) of each residue from python's enumerator.
## Whenever a residue is found in the string, its index is recorded where the first index of the string is set to '1'. Then, PyMOL selection and aesthetic editing occurs
## to create a molecular image summary of structural mass spectrometry date [NOTE this is originally developed for mapping data from EDC/GEE labeling experiments but is usable
## for any covalent labeling experiment that is intended to measure SASA changes.] 

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
## NOTE: "Bad" residue selections are detected if cmd.count_atom("sele") == 0. The moment a residue is detected like this,
          that entire protein is logged in an error file as failed due to bad selection algebra, and the loop skips to the next query.
"""

## IMPORT LIBRARIES
from datetime import date   # Import date class from datetime module
from pymol import cmd       # PyMOL's methods and commands; this is required for the script to use PyMOL's functions.
import pymol                # this one might be unnecessary since we specifically need pymol.cmd()
import re                   # methods for string editing
import csv                  # methods for handling csv file i/o
import time                 # methods for tracking efficiency of the code (CPU time)
import os                   # methods for directory handling
import sys                  # methods for taking command line arguments. Script's name is sys.arg[0] by default when -c flag is used
import matplotlib           # for graphing; TODO check if this is how to import matplotlib

## IMPORTANT VARIABLES
fasta = ""                  # Treat as a constant. String literal that stores the entire amino acid fasta sequence of the current PDB file
modSite = ""                # string literal pulled from a tab delimited text result file from PD; stores current `Modification Site` field value
modSeq = ""                 # also a string literal pulled from same; stores current `Annotated Sequence` field value
PDindex = ""                # also same as above; stores current `Position in Protein` value
adjustedIndex = ""          # stores a value from `PDindex` that has been corrected to the value of a real sequence position (ex., StrepII-HA fusion tag adds 12 amino acids, which PD would define as native aa positions)
tagLength = 17              # length of the N-terminal StrepII::HA epitope tag on AHA2. Position values <= 0 should be considered physiologically irrelevant.

## 1) NOTE TODO Add `METHOD` code block here: Open text file, iterate through rows in the columns for `Annotated Sequence` and `Modification Site`, correct the mod site and return those values for use in Pymol
def getMods(resultFile):
    
    """
    Extend code here for gathering important values above. Populate a list of modification sites and positions, as described in TODO messages:
    ##TODO 1) Extract all data necessary from PSM result files and merge into one master file; save this master file "output file 1"
    ##TODO 2) Extract only the following columns to a new file: `Annotated Sequence`, `Modifications`, `Charge`, `Precursor Abundance`, `Spectrum File`, and `File ID`; save this file "output file 2"
    ##TODO 3) From "output file 2", store only the unique records (rows) into a 3rd output file. Into the same output file, add a new column marking the presence of GEE modifications.
    ##TODO 4) Using the 3rd output file, store each unique `Annotated Sequence` in a nested dictionary structure, as the key. Then,
              populate each key's value with a 4-item list: `File ID`, `Modded Abundance` per sequence, `Unmodded Abundance` per sequence, and `Percent Modded`.
              Print these data to a separate `output file #4`.
    ##TODO 5) 
    """
    
    return


## 2) NOTE TODO Add `METHOD` code block here: Use values returned by `getMods()` to map mod sites on a PDB file object
def mapModsToProtein():

    """
    Extend code here for selecting a residue position in the PDB file object and manipulating it to mark it as a modified site (aka., re-color it and show it as spheres or sticks), , as described in TODO messages:
    """

    return

##################################################################################################################################
## 3) Driver code. Opens file reader object, imports values from a text file and manipulates them (calls `getMods()`),
##    and loads a PDB session file to perform aesthetic operations on the object (calls `mapModsToProtein()`) 
"""
Opens pymol file of interest (or fetches structure); using values gathered by getMods() it sends to mapModsToProtein(), which makes selections and changes to
the pymol file and saves the changes to a new filename (such as original filename with date, or date only)
"""

## NOTE the code below works as desired - a pre-existing .pse file is loaded, specific residues of a named type are visually modified, and the changes
## NOTE are saved to a new .pse file in the current working directory. FIXME Proof of principle here is complete, now all instances of "resn Glu" need to be replaced by
## a variable that iterates its values through a list of actual label events (defined by current resi value and whether that resi meets a minimum abundance threshold to be marked.)
## NOTE probably the best way to do that is by building a 2D list or dictionary containing actual D and E resi values and to then associate those with resi from the output text file,
## rather than by crunching numbers to correct PD's resi values. Then they can be referenced directly to retrieve their real resi and their labeling abundance.
today = date.today()
today = str(today)
sessionFile = str(sys.argv[3].upper())
outFileName = "" + today + "_" + sessionFile
cmd.load(sessionFile)
cmd.show("spheres", "resn Glu")
cmd.color("cyan", "resn Glu")
cmd.save(outFileName)