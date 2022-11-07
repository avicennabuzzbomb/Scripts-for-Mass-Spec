###################### BulkSubmit.sh: a script that writes a custom Submit file and executable shell script pair for any given set of jobs #########################################################################

## USAGE (Note that arg1 is a single PDB ID or a file containing a list of PDB IDs): $ ./BulkSubmit.sh <arg1>
## where <arg1> is the argument variable $1; write the name of a text file (list of PDB IDs to run) as the value of this variable
## writes a custom executable file for each PDB ID in a list of PDB IDs, then submits a job to the executive node with each iteration of the loop.

## Executable filename never changes, so that script is just automatically rewritten and then submitted with each job
## Submission filename does change depending on the job submitted and date

d=$(date +%Y-%m-%d)  

## Subsume everything below this comment into a for loop, which iterates through a list file of PDB ID

### CHECK FOR VALID ARGUMENTS and BEGIN WRITING EXECUTABLE SCRIPT ###
if [ $# -gt 0 ]; then

	echo "Modifying executable code..."  # informs user

	## Begin writing/editing the current executable file:
	echo '#!/bin/bash' > querySASA.sh

	# Codeline: untar your Python installation. Make sure you are using the right version!
	echo 'tar -xzf python37.tar.gz' >> querySASA.sh

	# Codelines: make sure the script will use your Python installation,
	# and the working directory as its home location
	echo "export PATH=\$PWD/python/bin:\$PATH" >> querySASA.sh
	echo "export PYTHONPATH=\$PWD" >> querySASA.sh
	echo "export HOME=\$PWD" >> querySASA.sh

	# Codeline: move to the software directory to use pymol as interpreter
	echo "cd /software/mrblackburn/pymol/" >> querySASA.sh

	# Codeline: run the SASA script
	echo "./pymol -c \$_CONDOR_SCRATCH_DIR/SASAquatch.py \$1 \$2 \$3" >> querySASA.sh

	# grant necessary permission (lets condor know that this is an executable file)
	chmod +x querySASA.sh
	### EXECUTABLE SCRIPT IS NOW FINISHED ###	

	### NOW WRITE THE SUBMIT FILE THAT WILL REFERENCE THE EXECUTABLE AND LIST THE JOBS TO RUN ###
	# "static" variables
	echo "Generating submit file..."  # informs user
	echo "# Name: SASA"$d".sub | written with list file"$1 > SASA$d.sub
	echo "# submit file for getting solvent accessible surface areas of the residues of target PDB IDs." >> SASA$d.sub
	echo "# This submit file was written with listname "$1 >> SASA$d.sub
	echo "#" >> SASA$d.sub
	echo "universe = vanilla" >> SASA$d.sub
	echo "requirements = (HasCHTCSoftware == true)" >> SASA$d.sub
	echo "executable = querySASA.sh" >> SASA$d.sub
	echo "should_transfer_files = YES" >> SASA$d.sub
	echo "when_to_transfer_output = ON_EXIT" >> SASA$d.sub
	echo "#" >> SASA$d.sub
	echo "#" >> SASA$d.sub
	echo "# transfer_input_files = file1,/absolute/pathto/file2,etc" >> SASA$d.sub

	# TODO This is where detection of pdb load files and listing of load files to send to SQUID needs to be extended!!
	files2transfer=$(echo "transfer_input_files = SASAquatch.py, http://proxy.chtc.wisc.edu/SQUID/chtc/python37.tar.gz")
	# Something like... if current line in list file contains the char ".", append its name to a string variable via echo, which stores the string that
	# lists the input files!

	## first check whether coordinate files are to be loaded with the job. Make sure the listed files are also present in the working dir where this is running!
	# NOTE if this does not work, comment out this block and switch back to the original echo statement directly below this comment:

	
	#echo "transfer_input_files = SASAquatch.py, http://proxy.chtc.wisc.edu/SQUID/chtc/python37.tar.gz" >> SASA$d.sub
	

	while IFS= read line; do
		
		name=$(echo $line)
		
		if [[ $line == *"."* ]]; then
			files2transfer=$(echo "$files2transfer, $line")
		fi

	done < "$1"

	# now echo the filename string to the submit file
	echo $files2transfer >> SASA$d.sub

	# request number of cpu's, memory and disk size to use (may change during the job run)
	echo "request_cpus = 2" >> SASA$d.sub
	echo "request_memory = 2GB" >> SASA$d.sub
	echo "request_disk = 300MB" >> SASA$d.sub

	### EXTEND THE SUBMIT FILE BY ADDING QUERIES AND QEUE COMMANDS (FOR-LOOP HERE FOR REPEATEDLY WRITING QUERIES FROM A LIST TO THE SUBMIT FILE) ###

	# specify the list file as an input to be read through line by line; correctly reads delimiters as part of the same line.
	while IFS= read line; do
		
		## to simplify downstream file handling, first convert entire argument string to uppercase
		#line=$(echo ${line^^})
		
		## write to the submit file
		echo "" >> SASA$d.sub
		echo 'arguments =' $line >> SASA$d.sub   # add query
		echo 'log = log/SASA_'$line'_$(Cluster).log' >> SASA$d.sub
		echo 'error = errors/SASA_'$line'_$(Cluster)_$(Process).err' >> SASA$d.sub
		echo 'output = output/SASA_'$line'_$(Cluster)_$(Process).out' >> SASA$d.sub
		echo "queue" >> SASA$d.sub	       # add queue request

	done < "$1"

	### SUBMIT FILE AND EXECUTABLE CODE (WHICH CALLS THE SASA SCRIPT) ARE NOW FINISHED AND READY FOR SUBMISSION ###

	echo "Done writing executable and submit files."

	### NOW SUBMIT THIS BATCH OF JOBS TO THE CLUSTER ###
	condor_submit SASA$d.sub
	echo "Done - submitted queries shown below:"
	condor_q

else
	echo "No query (ies) have been specified; this script must be invoked with at least 1 argument specifying one or more queries to submit as jobs."
fi
