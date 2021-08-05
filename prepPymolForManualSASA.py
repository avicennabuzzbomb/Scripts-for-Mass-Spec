##################### handles repetitive commands when fetching a new structure for SASA command testing ############################

##################### run this from within a pymol session, NOT in the background! Only useful when manually handling a session. ####

from pymol import cmd     # PyMOL's methods and commands; this is required for the script to use PyMOL's functions.
#import sys                # methods for taking command line arguments. Script's name is sys.arg[0] by default when -c flag is used

#query=str(sys.argv[1].upper())


# a set to capture iterate's `resn` list, called by the current `resi`, and clear duplicates; and a string variable to capture the `resn`
resPosition = set()
resname = set()
currRes = ""

# the usual chain isolation commands
cmd.reinitialize()
cmd.fetch("5ksd", async=0)
cmd.remove("Chain B")
cmd.remove("het")
cmd.set("seq_view", 1)

# prints the fasta chain(s) present in the session
for chain in cmd.get_chains():
    fasta = cmd.get_fastastr("Chain " + str(chain)) # gets the associated fasta with each chain ID     
    print(fasta)

# get `resi` string literals
cmd.iterate("chain A", 'resPosition.add(resi)')

# remove duplicates and sort to match fasta index 
resPosition = list(resPosition) # re-casting the set to a list, then back to a set erases duplicates, because sets do not keep duplicate elements
resPosition = list(set(resPosition)) 
resPosition.sort()          # sort() first to apply numerical character sorting
resPosition.sort(key=len)   # sort(key=len) uses the string element's length as the key to sort against, so now both sort rules apply

# use a resPosition element (a `resi` string literal) to call its corresponding `resn`
cmd.select("resi " + resPosition[0])    # NOTE in the full script `resi` is extracted the same way from a List that is passed
cmd.iterate("sele", 'resname.add(resn)')
print("\n`resname` as a raw set:")
print(resname)

# try popping the value out of the set
currRes=resname.pop()
print("\nthe `resn` popped from `resname`: " + currRes + "\n`resname` contents are:")
print(resname)

# prepare the session for SASA calculations
cmd.set("dot_solvent", 1)
cmd.set("dot_density", 4)
print("\nReady for SASA calculations: dot_solvent = TRUE, dot_density = 4")

# calculate SASA with the current setup
print("\nCalculating total SASA of first residue VAL using a `resi` selection method: ")
totSASA = cmd.get_area("resi " + resPosition[0] + " and Chain A")
print(totSASA)

print("\nCalculating sidechain SASA of first residue VAL using a `resi` selection method: ")
sideSASA = cmd.get_area("resi " + resPosition[0] + " and Chain A and sidechain")
print(sideSASA)