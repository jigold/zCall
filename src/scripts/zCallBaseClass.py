#! /usr/bin/env python

# Iain Bancarz, ib5@sanger.ac.uk, January 2013

# Define 'base' class containing useful methods for zcall scripts

class zCallBase:

    def __init__(self):
        pass

    def findMAF(self, nAA, nBB):
        # find minor allele frequency
        maf = None
        if nAA > nBB:
            maf = (nAB + 2 * nBB) / float(2*(nAA + nAB + nBB))
        else:
            maf = (nAB + 2 * nAA) / float(2*(nAA + nAB + nBB))
        return maf

    def readThresholds(self, inPath):
        # read a thresholds.txt file
        thresholdsX = []
        thresholdsY = []
        for line in open(inPath, 'r'):
            line = line.replace("\n", "")
            if line.find("Tx") != -1:
                continue
            else:
                fields = line.split("\t")
                #snp = fields[0]
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
