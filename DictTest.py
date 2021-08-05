import sys

AA_letterCode = {'R' : ('ARG'), 'H' : ('HIS'), 'K' : ('LYS'), 'D' : ('ASP'), 'E' : ('GLU'), 'S' : ('SER'), 'T' : ('THR'), 'N' : ('ASN'), 'Q' : ('GLN'), 'C' : ('CYS'),   
                 'G' : ('GLY'), 'P' : ('PRO'), 'A' : ('ALA'), 'V' : ('VAL'), 'I' : ('ILE'), 'L' : ('LEU'), 'M' : ('MET'), 'F' : ('PHE'), 'Y' : ('TYR'), 'W' : ('TRP'),
                 'ARG' : ('ARG'), 'HIS' : ('HIS'), 'LYS' : ('LYS'), 'ASP' : ('ASP'), 'GLU' : ('GLU'), 'SER' : ('SER'), 'THR' : ('THR'), 'ASN' : ('ASN'), 'GLN' : ('GLN'), 'CYS' : ('CYS'),   
                 'GLY' : ('GLY'), 'PRO' : ('PRO'), 'ALA' : ('ALA'), 'VAL' : ('VAL'), 'ILE' : ('ILE'), 'LEU' : ('LEU'), 'MET' : ('MET'), 'PHE' : ('PHE'), 'TYR' : ('TYR'), 'TRP' : ('TRP'),} 
                # NOTE this Dictionary is used only by non-"ALL" queries

arg = str(sys.argv[1].upper())

print("You asked for " + arg + " and you got " + AA_letterCode[arg] + " now get lost.")