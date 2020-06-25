##### This script will produce a list of interesting proteins that match a desired set of functional/structural attributes
## Either 1) scrape the pdb directly, or 2) scrape Uniprot (more organized, better annotated) for a starting -ome, then search
## PDB for entries that have a match in the list. (H+-ATPases for all species -> ome list -> cross-reference PDB -> add fetchables to list
## if they meet a minimum resolution and are non-redundant. This list will be fed to the Cluster with SASAquatch.py)