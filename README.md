# **Scripts for analyzing protein structure files.** #

1. `SASAquatchy.py` calculates solvent accessible surface areas in a PyMOL subshell. Can be invoked with arguments in the terminal or PyMOL's command line.\
    *From the terminal:* run with `pymol -c SASAquatch.py <PDB_ID> <depth>`. For example, with structure file 4h1w and calculating for all residues in the protein, run with:\
     `pymol -c SASAquatch.py 4H1W all`.\\
    *From PyMOL's command line:* First, make sure PyMOL is running in the directory where `SASAquatch.py`. Directory navigation can be done in PyMOL's command line exactly the same way in Windows command line or Bash. Then, run with:\
    `run SASAquatch.py <PDB_ID> <depth>`.

2. `NunuDrives.py` Reads arguments from a batch file and passes them to a second script. Currently, used to run SASAquatch.py on large PDB ID lists. Run only from a terminal with `python NunuDrives.py`.

3. `batch.txt` a text file with arguments (delimiter = ','); changes periodically as I run tests.

**Directory folders** contain other scripts, input files and output files for handling chemical crosslinking and/or footprinting data from mass spectrometry experiments. These are currently unfinished; development will begin again in the near future.