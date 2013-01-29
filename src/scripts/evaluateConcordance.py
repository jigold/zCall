#! /usr/bin/env python

# Iain Bancarz, ib5@sanger.ac.uk, January 2013

# Evaluate concordance for a list of GTC files

from GTC import *
from BPM import *
from EGT import *
from zCallBaseClass import zCallBase

class evaluator(zCallBase):

    def __init__(self, bpmPath, egtPath, threshPath):
        self.bpm = BPM(bpmPath)
        self.egt = EGT(egtPath)
        self.thresholds = readThresholds(threshPath)

    def evalSingle(self, gtcPath):
        gtc = GTC(gtcPath, self.bpm.normID)


    def call(self):
        pass
        

    def normalizeCall(self, call, nAA, nBB):
        ## Normalization:  Flip genotype call so 1 is always the common allele homozygote and 3 is the minor allele homozygote
        if nBB > nAA: 
            if call == 1:
                call = 3
            elif call == 3:
                call = 1
        return call

    def includeSNP(self, i, nAA=None, nBB=None):
        # should ith SNP be included in call rate calculation?
        # require autosomal snp with MAF>=5%
        # want at least 10 points in each homozygote cluster
        include = True
        chrom = self.bpm.chr[i]
        if nAA==None: nAA = self.egt.nAA[i]
        if nBB==None: nBB = self.egt.nBB[i]
        maf = self.findMAF(nAA, nBB)
        if maf < 0.05 or chrom == "X" or chrom == "Y" or nAA < 10 or nBB < 10:
            include = False
        return include

    def countCallTypes(self, gtc, thresholdsX, thresholdsY):
        # based on method in sampleConcordance.py
        # call codes: 0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        numSNPs = gtc.numSNPs
        counts = {}
        for i in range(4):
            for j in range(4): counts[(i,j)] = 0
        for i in range(numSNPs):
            nAA = self.egt.nAA[i]
            nAB = self.egt.nAB[i]
            nBB = self.egt.nBB[i]
            if not self.includeSNP(i, nAA, nBB): continue
            normX = gtc.normXintensities[i]
            normY = gtc.normYintensities[i]
            Tx = thresholdsX[i]
            Ty = thresholdsY[i]
            if Tx == "NA" or Ty == "NA": continue
            A = self.bpm.A[i]
            B = self.bpm.B[i]
            origCall = self.normalizeCall(gtc.genotypes[i], nAA, nBB)
            # re-call with zcall thresholds
            if normX < Tx and normY < Ty: ## Lower left quadrant
                newCall = 0
            elif normX >= Tx and normY <= Ty: ## Lower right quadrant
                newCall = self.normalizeCall(1, nAA, nBB)
            elif normX < Tx and normY >= Ty: ## Upper left quadrant
                newCall = self.normalizeCall(3, nAA, nBB)
            else: ## Upper right quadrant
                newCall = 2
            counts[(origCall, newCall)] += 1
        return counts
