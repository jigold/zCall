#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# May 8th, 2012

import sys
from GTC import *
from BPM import *
from optparse import OptionParser

### Parse Inputs from Command Line
parser = OptionParser()
parser.add_option("-B","--bpm",type="string",dest="bpm",action="store",help="bpm file with normID and chr, pos")
parser.add_option("-G","--gtc",type="string",dest="gtc",action="store",help="gtc file to recall")
parser.add_option("-T","--thresholds",type="string",dest="thresholds",action="store",help="file with threshold definitions")
(options, args) = parser.parse_args()

if options.bpm == None:
    print "specify BPM file path with -B"
    sys.exit()

if options.gtc == None:
    print "specify GTC file path with -G"
    sys.exit()

if options.thresholds == None:
    print "specify thresholds.txt file with -T"
    sys.exit()

### Parse sample name from input gtc file name
sampleName = options.gtc.split("/")
sampleName = sampleName[len(sampleName) - 1]
sampleName = sampleName.split(".")[0]

### Initialize GTC and BPM file classes
bpm = BPM(options.bpm)    
gtc = GTC(options.gtc, bpm.normID)

### Get Number of SNPs
numSNPs = gtc.numSNPs

### Parse thresholds
thresholdsX = []
thresholdsY = []
for line in open(options.thresholds, 'r'):
    line = line.replace("\n", "")
    if line.find("Tx") != -1:
        continue
    else:
        fields = line.split("\t")

        snp = fields[0]
        if fields[1] == "NA" or fields[2] == "NA":
            tx = fields[1]
            ty = fields[2]
        else:
            tx = float(fields[1])
            ty = float(fields[2])

        thresholdsX.append(tx) # python list with X thresholds
        thresholdsY.append(ty) # python list with Y thresholds
    
### Make PED file
out = [sampleName, sampleName, "0", "0", "-9", "-9"] #output holder in python list; have no sample information so use "0" and "-9" for mid, pid, gender, case/control status

for i in range(numSNPs):
    normX = gtc.normXintensities[i]
    normY = gtc.normYintensities[i]
    Tx = thresholdsX[i]
    Ty = thresholdsY[i]
    A = bpm.A[i] # Get what A allele is from bpm.csv file [A,G,T,C]
    B = bpm.B[i] # Get what B allele is from bpm.csv file [A,G,T,C]
    
    origCall = gtc.genotypes[i] # 0 - "No Call", 1 - AA, 2 - AB, 3 - BB

    if origCall == 0 and Tx != "NA" and Ty != "NA": ## ONLY RECALL NO CALLS!!!!
        if normX < Tx and normY < Ty: ## Lower left quadrant
            out.append("0")
            out.append("0")
            
        elif normX >= Tx and normY <= Ty: ## Lower right quadrant
            out.append("A")
            out.append("A")
#             out.append(A) # uncomment to output A and B alleles from manifest file; make sure to comment lines above
#             out.append(A)
            
        elif normX < Tx and normY >= Ty: ## Upper left quadrant
            out.append("B")
            out.append("B")
#             out.append(B) # uncomment to output A and B alleles from manifest file; make sure to comment lines above
#             out.append(B)
            
        else: ## Upper right quadrant
            out.append("A")
            out.append("B")
#             out.append(A) # uncomment to output A and B alleles from manifest file; make sure to comment lines above
#             out.append(B)

    else: ## Output original caller's call
        if origCall == 0: ## NC is original call
            out.append("0")
            out.append("0")
            
        if origCall == 1: ## AA is original call
            out.append("A")
            out.append("A")
#             out.append(A) # uncomment to output A and B alleles from manifest file; make sure to comment lines above
#             out.append(A)
            
        if origCall == 2: ## AB is original call
            out.append("A")
            out.append("B")
#             out.append(A) # uncomment to output A and B alleles from manifest file; make sure to comment lines above
#             out.append(B)
            
        if origCall == 3: ## BB is original call
            out.append("B")
            out.append("B")
#             out.append(B) # uncomment to output A and B alleles from manifest file; make sure to comment lines above
#             out.append(B)

## Output to std out the new calls in PED format
print " ".join(out)




