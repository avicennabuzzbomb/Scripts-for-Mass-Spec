# remove modification symbols from the string so they don't wrongly affect the amino acid calculation

import re

seq = "[R].FFDNHIFAScSDDNILR.[F]"
cleanedStr = ""

print(seq)

seq = re.sub('[a-z]', '', seq)

print(seq)


## TODO comme back to this later, needed for completing the peak predictor script for identifying valid m/z shifts in labeled peptides