#! /usr/bin/python

"""Object to read and contain zcall thresholds

Author:  Iain Bancarz, ib5@sanger.ac.uk, January 2013
"""


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
