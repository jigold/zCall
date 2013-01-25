#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# May 8th, 2012

import sys
from optparse import OptionParser
from EGT import *

### Parse Inputs from Command Line
parser = OptionParser()
parser.add_option("-B","--betas",type="string",dest="betas",action="store",help="betas.txt file path")
parser.add_option("-E","--egt",type="string",dest="egt",action="store",help=".EGT file path")
parser.add_option("-I","--minint",type="string",dest="minIntensity",action="store",help="minimum mean int signal for comm hom cluster")
parser.add_option("-Z","--z",type="string",dest="z",action="store",help="z-score threshold")
(options, args) = parser.parse_args()

if options.egt == None:
    print "specify EGT file path with -E"
    sys.exit()

if options.z == None:
    options.z = 7

if options.betas == None:
    print "specify betas.txt file with -B"
    sys.exit()

if options.minIntensity != None:
    options.minIntensity = float(options.minIntensity)
else:
    options.minIntensity = 0.2

options.z = int(options.z)

### Print header line to std out
head = ["SNP", "Tx", "Ty"] 
print "\t".join(head)


### Initialize EGT and BPM class
egt = EGT(options.egt)


### Parse betas.txt file
### Order is as follows:
### meanY ~ meanX
### meanX ~ meanY
### sdY ~ sdX
### sdX ~ sdY
beta0 = [] # list container for beta intercept
beta1 = [] # list container for beta of slope
for line in open(options.betas, 'r'):
    line = line.replace("\n", "")
    if line.find("Beta0") == -1:
        fields = line.split("\t")
        beta0.append(float(fields[1]))
        beta1.append(float(fields[2]))


### Find thresholds for each site
### AA is always quadrant 4 by definition, BB is always quadrant 1 by definition
### because X always tags A and Y always tags B
        
numSNPs = egt.numCodes # extract number of snps from EGT class
z = options.z # z-score threshold to use

for i in range(numSNPs):

    # Get SNP Name
    snp = egt.names[i]
    
    # Get number of points in each homozygote genotype cluster
    nAA = egt.nAA[i]
    nBB = egt.nBB[i]

    # Extract mean and standard deviations from EGT class
    meanXAA = egt.meanXAA[i]
    meanXBB = egt.meanXBB[i]
    devXAA = egt.devXAA[i]
    devXBB = egt.devXBB[i]

    meanYAA = egt.meanYAA[i]
    meanYBB = egt.meanYBB[i]
    devYAA = egt.devYAA[i]
    devYBB = egt.devYBB[i]

    # Calculate Thresholds depending on which homozygote cluster is tagging the common allele for that SNP
    if nAA <= 2 and nBB <= 2: # Not enough points in common allele homozygote cluster
        Tx = "NA"
        Ty = "NA"

    else:
        if nAA >= nBB:
            if meanXAA < options.minIntensity: # site has less than min. intensity to recall
                Tx = "NA" # mark with "NA" so skip new genotype calls in zCall.py
                Ty = "NA"
            else:
                Ty = meanYAA + z * devYAA
                meanXBB = beta1[1]*meanYAA + beta0[1] # Solve for the mean of the minor allele hom. cluster based on betas and mean of common allele hom. cluster
                devXBB = beta1[3]*devYAA + beta0[3] # Solve for the sd of the minor allele hom. cluster based on betas and sd of common allele hom. cluster
                Tx = meanXBB + z * devXBB # Use inferred mean and sd to find Tx

        if nAA < nBB:
            if meanYBB < options.minIntensity: # site has less than min. intensity to recall
                Tx = "NA" # mark with "NA" so skip new genotype calls in zCall.py
                Ty = "NA"
            else:
                Tx = meanXBB + z * devXBB
                meanYAA = beta1[0] * meanXBB + beta0[0] # Solve for the mean of the minor allele hom. cluster based on betas and mean of common allele hom. cluster
                devYAA = beta1[2] * devXBB + beta0[2] # Solve for the sd of the minor allele hom. cluster based on betas and sd of common allele hom. cluster
                Ty = meanYAA + z * devYAA # Use inferred mean and sd to find Ty
            

    # Write thresholds to std out
    out = [snp, Tx, Ty]
    out = [str(o) for o in out]
        
    print "\t".join(out)

