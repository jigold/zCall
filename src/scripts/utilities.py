#! /usr/bin/env python

# Iain Bancarz, ib5@sanger.ac.uk, January 2013

import json
from GTC import *
from BPM import *
from EGT import *

class SharedBase:
    """Constants and simple functions shared across multiple classes"""

    Z_KEY = 'BEST_Z'
    T_KEY = 'BEST_THRESHOLDS'
    M_KEY = 'SAMPLE_METRICS'

    def readSampleJson(self, inPath):
        """Read sample GTC paths from .json file used by genotyping pipeline"""
        samples = json.loads(open(inPath).read())
        gtc = []
        for sample in samples: gtc.append(sample['result'])
        return gtc
        

class CallingBase(SharedBase):
    """ 'Base' class containing useful methods for zcall subclasses """

    def __init__(self, bpmPath, egtPath, threshPath=None):
        self.bpm = BPM(bpmPath)
        self.egt = EGT(egtPath)
        self.snpTotal = self.egt.getTotalSNPs()
        if self.bpm.getTotalSNPs() != self.snpTotal:
            raise ValueError("ERROR: SNP totals in .egt and .bpm inputs differ")
        if threshPath!=None:
            self.thresholds = ThresholdContainer(threshPath)
        else:
            self.thresholds = None

    def call(self, gtc, i):
        """ re-call ith SNP in GTC file, using zcall thresholds
        
        call codes: 0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        """
        if self.thresholds==None:
            raise ValueError("Must specify thresholds before calling!")
        normX = gtc.normXintensities[i]
        normY = gtc.normYintensities[i]
        Tx = self.thresholds.getX(i)
        Ty = self.thresholds.getY(i)
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

    def setThresholds(self, thresholds):
        """Set thresholds to given ThresholdContainer object"""
        self.thresholds = thresholds


class ThresholdContainer:
    """Class to contain thresholds, read from thresholds.txt file"""

    def __init__(self, inPath):
        (self.x, self.y) = self.readThresholds(inPath)
        
    def getX(self, i):
        """Get x threshold for the ith SNP"""
        return self.x[i]

    def getY(self, i):
        """Get y threshold for the ith SNP"""
        return self.y[i]

    def readThresholds(self, inPath):
        """ Read a thresholds.txt file; return lists of x and y thresholds """
        thresholdsX = []
        thresholdsY = []
        for line in open(inPath, 'r'):
            line = line.replace("\n", "")
            if line.find("Tx") != -1:
                continue
            else:
                fields = line.split("\t")
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
        return (thresholdsX, thresholdsY)
