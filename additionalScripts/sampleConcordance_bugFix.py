#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# April 5th, 2012

import sys
from GTC import *
from BPM import *
from EGT import *
from optparse import OptionParser

### Parse Inputs from Command Line
parser = OptionParser()
parser.add_option("-B","--bpm",type="string",dest="bpm",action="store",help="bpm file with normID and chr, pos")
parser.add_option("-E","--egt",type="string",dest="egt",action="store",help="egt file")
parser.add_option("-G","--gtc",type="string",dest="gtc",action="store",help="gtc file to recall")
parser.add_option("-T","--thresholds",type="string",dest="thresholds",action="store",help="file with threshold definitions")
(options, args) = parser.parse_args()

if options.bpm == None:
    print "specify BPM file path with -B"
    sys.exit()

if options.egt == None:
    print "specify EGT file path with -E"
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
egt = EGT(options.egt)
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
        if fields[1] != "NA":
            tx = float(fields[1])
        else:
            tx = fields[1]
        if fields[2] != "NA":
            ty = float(fields[2])
        else:
            ty = fields[2]

        thresholdsX.append(tx)
        thresholdsY.append(ty)
    
### Count concordance between GenCall and zCall for common sites (MAF > 5%)
### order is (genCall, zCall)
counts = {(0,0):0,(0,1):0,(0,2):0,(0,3):0,
          (1,0):0,(1,1):0,(1,2):0,(1,3):0,
          (2,0):0,(2,1):0,(2,2):0,(2,3):0,
          (3,0):0,(3,1):0,(3,2):0,(3,3):0}

for i in range(numSNPs):

    chr = bpm.chr[i]

    nAA = egt.nAA[i]
    nAB = egt.nAB[i]
    nBB = egt.nBB[i]

    if nAA > nBB:
        maf = (nAB + 2 * nBB) / float(2*(nAA + nAB + nBB))
    elif nBB >= nAA:
        maf = (nAB + 2 * nAA) / float(2*(nAA + nAB + nBB))

    if maf < 0.05 or chr == "X" or chr == "Y": # Only want to look at common, autosomal sites!!!
        continue

    if nAA < 10 or nBB < 10: # Make sure at least 10 points in both homozygote clusters
        continue
    
    normX = gtc.normXintensities[i]
    normY = gtc.normYintensities[i]
    Tx = thresholdsX[i]
    Ty = thresholdsY[i]
    A = bpm.A[i]
    B = bpm.B[i]

    if Tx == "NA" or Ty == "NA":
        continue
    
    origCall = gtc.genotypes[i] # 0 - "No Call", 1 - AA, 2 - AB, 3 - BB

    if nBB > nAA: ## Flip genotype call so 1 is always the common allele homozygote and 3 is the minor allele homozygote
        if origCall == 1:
            origCall = 3
        elif origCall == 3:
            origCall = 1
            

    if normX < Tx and normY < Ty: ## Lower left quadrant
        newCall = 0
            
    elif normX >= Tx and normY <= Ty: ## Lower right quadrant
        if nAA >= nBB:
            newCall = 1
        else: ## Flip genotype call so 1 is always the common allele homozygote and 3 is the minor allele homozygote
            newCall = 3
            
    elif normX < Tx and normY >= Ty: ## Upper left quadrant
        if nAA >= nBB:
            newCall = 3
        else: ## Flip genotype call so 1 is always the common allele homozygote and 3 is the minor allele homozygote
            newCall = 1
            
    else: ## Upper right quadrant
        newCall = 2

    
    counts[(origCall, newCall)] += 1 # add one to counts for origCall, newCall combination


### Print output to std out

head = ["#SampleName", "0_0", "0_1", "0_2", "0_3", "1_0", "1_1", "1_2", "1_3", "2_0", "2_1", "2_2", "2_3", "3_0", "3_1", "3_2", "3_3"]
print "\t".join(head)

out = [sampleName, counts[(0,0)], counts[(0,1)], counts[(0,2)], counts[(0,3)],
       counts[(1,0)], counts[(1,1)], counts[(1,2)], counts[(1,3)],
       counts[(2,0)], counts[(2,1)], counts[(2,2)], counts[(2,3)],
       counts[(3,0)], counts[(3,1)], counts[(3,2)], counts[(3,3)]]

out = [str(o) for o in out]

print "\t".join(out)
