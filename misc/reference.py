# extracts cells of interest from csv files of MS data, and outputs to a new csv
# this is NOT portable in current form (directory handling); use instead as reference code.
# Original author: Kento Haruta

import glob
import csv
from operator import itemgetter


column = int(input("enter column"))
ooga = int(input("top how many???"))
top = []
bigcount = 0
count = 0
biggest = ""
biglist = []
foldername = input("enter a folder name")
files = glob.glob('/Users/kissmyface/Desktop/TabularData/'+foldername+'/*.csv')#put the directory here
finalist = []
def getUniqueItems(L):
   result = []
   for item in L:
       if [item[0],item[1].lower().strip()] not in result:
           result.append([item[0],item[1].lower().strip()])
   return result

for i in files:
   with open(i) as f:
       data = csv.reader(f)
       biglist=(biglist+ list(csv.reader(f)))
       
       
for b in biglist:
   count = 0
   for v in biglist:
       #print(b,"??!!??",v," |||||")
       if b[column-1].lower().strip()==v[column-1].lower().strip():
           count = count+1
   top.append([count,b[column-1]])
       
   if count>bigcount:
       bigcount = count
       biggest = b[column-1]
finalist = getUniqueItems(top)
finalist.sort(key=itemgetter(0),reverse=True)
for buy in range(ooga):
   print(finalist[buy])
