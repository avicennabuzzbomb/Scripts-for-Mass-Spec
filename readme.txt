These scripts are backed up on GitHub; scripts relating to protein SASA including helper shell scripts live on my HTCondor Account.

TO DOWNLOAD AN EXISTING COPY OF AN HTCONDOR FILE, MOVE TO THE DESIRED DIRECTORY ON THIS MACHINE, THEN TYPE:
$ scp mrblackburn@submit-1.chtc.wisc.edu:/home/mrblackburn/*sh ./    # in this example, this command requests a download of every file ending in the extension `.sh`.
								     # full filenames will also work

# This will prompt a request for my password; enter it here, and the downloads will save in the current working directory.


CONVERSELY, TO UPLOAD FILES TO MY HTCONDOR DIRECTORY, TYPE THE FOLLOWING:
$ scp <myFile> mrblackburn@submit-1.chtc.wisc.edu:/home/mrblackburn/myPreferredDirectory   # in this example, a file named myFile will be uploaded to a folder named 
											   # myPreferredDirectory in my home directory.

# This will again prompt a request for my password as above, prior to upload from my local machine to my HTCondor home directory
