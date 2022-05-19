## Write a script here for mapping all-atom distances to a single selection, identifying the shortest point to point distance and storing that.
## Just 1 example of how this would be useful: Find the inter-residue distances between the conserved proton binding site of AHA2, D684, and all other sidechains within
## some distance. Generating a list of inter-residue distances, then pulling out the shortest ones associated with each other non-D684 residue, would
## rapidly highlight the contours of the aqueous proton-binding pocket of AHA2.

# TODO: algorithm: 
# 0) Open a read-wrote object to stream data into an output file (.csv)
# 1) Import arguments: PDB ID, fasta, and selected residue(s) to map nearest-neighbors.
# 2) Make requisite selection (loop, if multiple selections to be used)
# 3) Store all of the selected object's atoms as a list.
# 4) Looping through that list, use cmd.distance with wildcards on the selected atom
#    to find all the closest atoms within a specified range of Angstroms (ie, dont cover the entire protein); store in array.
# 5) Search the array for smallest distances; print to output with resn "name"