#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# April 5th, 2012

import sys
from optparse import OptionParser
from EGT import *

### Parse Inputs from Command Line
parser = OptionParser()
parser.add_option("-E","--egt",type="string",dest="egt",action="store",help=".EGT file path")
(options, args) = parser.parse_args()

if options.egt == None:
    print "specify EGT file path with -E"
    sys.exit()

### Write header line
head = ["SNP", "meanX", "meanY", "sdX", "sdY", "nMinorHom", "nCommonHom"] 
print "\t".join(head) # Write header line

### Initialize EGT class
egt = EGT(options.egt)

## Iterate over each SNP in EGT object
numSNPs = egt.numCodes ## number of SNPs in EGT file
for i in range(numSNPs):

    # Get SNP Name
    snp = egt.names[i]
    
    # Get number of points in each genotype cluster
    nAA = egt.nAA[i]
    nAB = egt.nAB[i]
    nBB = egt.nBB[i]
    nTotal = nAA + nAB + nBB

    # Calculate Missing Rate (ignore SNPs that have less than 99% call rate)
    if float(nTotal) / float(egt.numPoints) < 0.99:
        continue

    # Make sure there are at least 10 points in each homozygote cluster
    if nAA < 10 or nBB < 10:
        continue
    
    # Calculate MAF
    if nAA > nBB:
        maf = (nAB + 2 * nBB) / float(2 * nTotal)
    elif nAA <= nBB:
        maf = (nAB + 2 * nAA) / float(2 * nTotal)

    # MAF check ( >5% MAF)
    if maf < 0.05:
        continue
    
    # Hardy-Weinberg Equilibrium Check (don't use site if p_hwe < 0.00001)
    chiCritical = 19.5 # p = 0.00001 for 1 DOF

    if nAA > nBB:
        p = 1.0 - maf
        q = maf        
        expAA = p**2 * nTotal
        expAB = 2 * p * q * nTotal
        expBB = q**2 * nTotal

    if nBB >= nAA:
        p = 1.0 - maf
        q = maf        
        expAA = q**2 * nTotal
        expAB = 2 * p * q * nTotal
        expBB = p**2 * nTotal
        
    chiSquare = ((nAA - expAA)**2 / float(expAA)) + ((nAB - expAB)**2 / float(expAB)) + ((nBB - expBB)**2 / float(expBB))
    if chiSquare > chiCritical:
        continue
    
    # Extract the mean and sd for each common allele homozygote clusters in the noise dimension
    meanXAA = egt.meanXAA[i]
    meanXBB = egt.meanXBB[i]
    devXAA = egt.devXAA[i]
    devXBB = egt.devXBB[i]

    meanYAA = egt.meanYAA[i]
    meanYBB = egt.meanYBB[i]
    devYAA = egt.devYAA[i]
    devYBB = egt.devYBB[i]

    if meanXAA >= meanYAA: ## AA is in the lower right quadrant
        meanY = meanYAA
        devY = devYAA
        meanX = meanXBB
        devX = devXBB

    elif meanXAA < meanYAA: ## AA is in the upper left quadrant; however, this should never be the case by definition X -> A, Y -> B
        meanY = meanYBB
        devY = devYBB
        meanX = meanXAA
        devX = devXAA

    if nAA >= nBB:
        out = [snp, meanX, meanY, devX, devY, nBB, nAA] # output array
    elif nBB > nAA:
        out = [snp, meanX, meanY, devX, devY, nAA, nBB] # output array
    out = [str(o) for o in out]

    print "\t".join(out)

