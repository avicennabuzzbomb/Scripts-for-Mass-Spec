# **Scripts for analyzing protein structure files.** #
`SASAquatchy.py` calculates solvent accessible surface areas in a PyMOL subshell.   Can be invoked with arguments in the command line.
`NunuDrives.py` Reads arguments from a batch file and passes them to a second script. Currently, used to run SASAquatch.py on large PDB ID lists.
`batch.txt` a text file with arguments (delimiter = ','); changes periodically as I run tests.
**Directories** contain scripts, input files and output files for handling chemical crosslinking and/or footprinting data from mass spectrometry experiments. These are currently unfinished; development will begin again in the near future.