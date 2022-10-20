##NOTE This is for testing split()

# create a List[] of amino acid ID's from a "," delimited string
aastring="D, E, glu, AsP"
print("The argument string is " + aastring)

# try splitting the string into a list of delimited elements (with specified delimiter removed)
aastring=aastring.split(", ")
for i in aastring:
    j = 0
    print("Element number " + str(j) + " in aastring is " + i.upper())
    j = j + 1

