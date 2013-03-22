#! /usr/bin/python

"""Write Plink .ped file from calls in GTC files; use to generate test data

Produces a single line, representing one sample, and appends to .ped file

Adapted from legacy zCall.py
Usage: writePEDfile.py [bpm path] [gtc path] [ped path]"""

import sys
sys.path.append(sys.path[0]+'/../zcall')
from BPM import *
from GTC import *

bpm = BPM(sys.argv[1])
gtc = GTC(sys.argv[2], bpm.normID)
outPath = sys.argv[3]
useManifest = False

### Parse sample name from input gtc file name
sampleName = sys.argv[2].split("/")
sampleName = sampleName[len(sampleName) - 1]
sampleName = sampleName.split(".")[0]

out = [sampleName, sampleName, "0", "0", "-9", "-9"] #output holder in python list; have no sample information so use "0" and "-9" for mid, pid, gender, case/control status
for i in range(gtc.getTotalSNPs()):
    if useManifest: 
        alleleA = bpm.A[i]
        alleleB = bpm.B[i]
    else:
        alleleA = "A"
        alleleB = "B"
    origCall = gtc.genotypes[i]     
    if origCall == 1: ## AA is original call
        out.append(alleleA)
        out.append(alleleA)
    elif origCall == 2: ## AB is original call
        out.append(alleleA)
        out.append(alleleB)
    elif origCall == 3: ## BB is original call
        out.append(alleleB)
        out.append(alleleB)
    else: ## NC is original call
        out.append("0")
        out.append("0") 
## Output to std out the new calls in PED format
outFile = open(outPath, 'a') # append to existing file
outFile.write(" ".join(out)+"\n")
outFile.close()
