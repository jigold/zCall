#! /usr/bin/env python

# Iain Bancarz, ib5@sanger.ac.uk, January 2013

from GTC import *
from BPM import *
from EGT import *
from thresholdContainer import ThresholdContainer

class zCallBase:

    """ 'Base' class containing useful methods for zcall subclasses """

    def __init__(self, threshPath, bpmPath, egtPath):
        self.thresholds = ThresholdContainer(threshPath)
        self.bpm = BPM(bpmPath)
        self.egt = EGT(egtPath)
        self.snpTotal = self.egt.getTotalSNPs()
        if self.bpm.getTotalSNPs() != self.snpTotal:
            raise ValueError("ERROR: SNP totals in .egt and .bpm inputs differ")

    def call(self, gtc, i):
        """ re-call ith SNP in GTC file, using zcall thresholds
        
        call codes: 0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        """
        normX = gtc.normXintensities[i]
        normY = gtc.normYintensities[i]
        Tx = self.thresholds.getX[i]
        Ty = self.thresholds.getY[i]
        call = None
        if normX < Tx and normY < Ty: ## Lower left quadrant
            call = 0
        elif normX >= Tx and normY <= Ty: ## Lower right quadrant
            call = 1
        elif normX < Tx and normY >= Ty: ## Upper left quadrant
            call = 3
        else: ## Upper right quadrant
            call = 2
        return call

    def findMAF(self, nAA, nBB, nAB):
        """ Find minor allele frequency """
        maf = None
        if nAA > nBB:
            maf = (nAB + 2 * nBB) / float(2*(nAA + nAB + nBB))
        else:
            maf = (nAB + 2 * nAA) / float(2*(nAA + nAB + nBB))
        return maf

    def normalizeCall(self, call, nAA, nBB):
        """Normalization:  Flip genotype call so 1 is always the common allele homozygote and 3 is the minor allele homozygote 

        Enforces convention that major allele is on X intensity axis
        Allele counts are taken from EGT object
        """
        if nBB > nAA: 
            if call == 1:
                call = 3
            elif call == 3:
                call = 1
        return call
