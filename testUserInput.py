

Output = open('Predicted_mixedPeaks.csv', 'w+') # create a new file object Output, open an empty file with it to 'w'rite into.
print("HEADER:: Predicted m/z for light-light(LL), light-heavy(LH), heavy-light(HL), heavy-heavy(HH) dipeptides.", file = Output)  
print("", file = Output)
COUNT = 1

response = 'y'   #default in order to enter the loop

##This large codeblock is designed for manual entry of sequences by the user (slow, mistake prone)
while response == 'Y' or response == 'y':
    seqA = input("Enter peptide sequence A: ")
    print("seqA is ",seqA)
    seqB = input("Enter peptide sequence B: ")
    print("seqB is ",seqB)
    mZ = input("Enter mZ of seqA-B: ")
    mZ = float(mZ)          #have to typecast because input defaults to str
    print("mZ is ",mZ)      
    charge = input("Enter charge, z: ")
    charge = int(charge)    #have to typecast because input defaults to str
    print("Charge: z = ",charge)

    dipeptide = seqA + "-" + seqB
    print("The dipeptide is ",dipeptide,"and its charge state is +",charge)
    print("")

    labeledPeps = calc_Labeled(seqA,seqB,mZ,charge) #get list of predicted sequences from calc_Labeled
    pepNames = ["LL","LH","HL","HH"]

    #This block prints the contents of the array of lists: label names, and predicted m/z to the new file.
    print(COUNT,":","Dipeptide sequence is ",dipeptide,"and its charge state is +",charge, file = Output)
    for x in range(len(pepNames)):
        print(pepNames[x], labeledPeps[x], file = Output)
    print("", file=Output)
    COUNT += 1

    # To user: continue?
    response = input("Analysis done. Enter another dipeptide?\n" "(Any key to continue; to stop, N or n.)")
 