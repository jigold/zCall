#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# April 5th, 2012

#The Illumina provided Code was provided as-is and with no warranty as to performance and no warranty against it infringing any other party's intellectual property rights.

import struct
import math

class EGT:
    ''' Class to parse an Illumina .egt file'''
    def __init__(self, file):
        self.f = open(file, 'rb')
        self.readHeaderData()
        self.readFileData()
        self.numPoints = max(self.nTotal) # estimate the number of samples used for clustering by looking at the maximum number of points used for clustering was
        
    def readFileData(self):
        self.version = self.getInt()
        self.opa = self.getString()
        self.numCodes = self.getInt() # number of SNPs in file
        self.parseClusterPositions() # get mean and standard deviations of each cluster
        self.readSNPnames() # get names of each snp
        
    def readHeaderData(self):
        self.FileVersion = self.getInt()
        self.GcVersion = self.getString()
        self.ClusterVersion = self.getString()
        self.CallVersion = self.getString()
        self.NormalizationVersion = self.getString()
        self.DateCreated = self.getString()
        self.Mode = self.getByte()
        self.Manifest = self.getString()

    def readSNPnames(self):
        self.names = []

        buffer = self.f.read(13 * self.numCodes) # skip snp quality scores
        
        for i in range(self.numCodes): # skip genotype scores
            gt = self.getString()

        for i in range(self.numCodes): # parse snp names
            snp = self.getString()
            self.names.append(snp)

    def getTotalSNPs(self):
        return self.numCodes

    def parseClusterPositions(self):
        '''
        Parse cluster positions for each site and transform from r,theta to x,y
        '''
        self.nAA = []
        self.nAB = []
        self.nBB = []
        self.nTotal = []
        
        self.meanXAA = []
        self.meanXAB = []
        self.meanXBB = []

        self.meanYAA = []
        self.meanYAB = []
        self.meanYBB = []

        self.devXAA = []
        self.devXAB = []
        self.devXBB = []

        self.devYAA = []
        self.devYAB = []
        self.devYBB = []        

        buffer = self.f.read((30*4)*self.numCodes)
        
        for i in range(self.numCodes):
            nAA = self.parseInt(buffer,(i*30*4) + 0) # number of points in AA cluster
            nAB = self.parseInt(buffer,(i*30*4) + 4) # number of points in AB cluster
            nBB = self.parseInt(buffer,(i*30*4) + 8) # number of points in BB cluster
            
            devRAA = self.parseFloat(buffer,(i*30*4) + 12)
            devRAB = self.parseFloat(buffer,(i*30*4) + 16)
            devRBB = self.parseFloat(buffer,(i*30*4) + 20)

            meanRAA = self.parseFloat(buffer,(i*30*4) + 24)
            meanRAB = self.parseFloat(buffer,(i*30*4) + 28)
            meanRBB = self.parseFloat(buffer,(i*30*4) + 32)

            devThetaAA = self.parseFloat(buffer,(i*30*4) + 36)
            devThetaAB = self.parseFloat(buffer,(i*30*4) + 40)
            devThetaBB = self.parseFloat(buffer,(i*30*4) + 44)

            meanThetaAA = self.parseFloat(buffer,(i*30*4) + 48)
            meanThetaAB = self.parseFloat(buffer,(i*30*4) + 52)
            meanThetaBB = self.parseFloat(buffer,(i*30*4) + 56)


            meanXAA, meanYAA, devXAA, devYAA = self.polarToEuclidean(meanRAA, devRAA, meanThetaAA, devThetaAA)
            meanXAB, meanYAB, devXAB, devYAB = self.polarToEuclidean(meanRAB, devRAB, meanThetaAB, devThetaAB)
            meanXBB, meanYBB, devXBB, devYBB = self.polarToEuclidean(meanRBB, devRBB, meanThetaBB, devThetaBB)

            self.nAA.append(nAA)
            self.nAB.append(nAB)
            self.nBB.append(nBB)
            self.nTotal.append(nAA + nAB + nBB)
            
            self.meanXAA.append(meanXAA)
            self.meanXAB.append(meanXAB)
            self.meanXBB.append(meanXBB)

            self.meanYAA.append(meanYAA)
            self.meanYAB.append(meanYAB)
            self.meanYBB.append(meanYBB)

            self.devXAA.append(devXAA)
            self.devXAB.append(devXAB)
            self.devXBB.append(devXBB)

            self.devYAA.append(devYAA)
            self.devYAB.append(devYAB)
            self.devYBB.append(devYBB)


    def polarToEuclidean(self, meanR, devR, meanTheta, devTheta):
        '''
        Manhattan distance conversion from r,theta to x,y
        '''
        varTheta = devTheta**2
        varR = devR ** 2
        
        piFactor = math.pi / float(2)

        A = -1 * (piFactor * meanR) * (1 + math.tan(piFactor * meanTheta))**-2 * (1 / float(math.cos(piFactor * meanTheta)**2))
        B = 1 / float( 1 + math.tan(piFactor * meanTheta))
        varX = (A * A * varTheta) + (B * B * varR)
        C = -1 * A
        D = float(1 - B)
        varY = (C*C * varTheta) + (D*D * varR)

        meanX = (meanR / float( 1 + (math.tan(meanTheta * piFactor))))
        meanY = (meanR - (meanR / float(1 + math.tan(meanTheta * piFactor))))
        devX = varX**(0.5)
        devY = varY**(0.5)

        return (meanX, meanY, devX, devY)                                                                                                      
        
    def getByte(self):
        '''
        read byte from current place in file
        '''
        line = self.f.read(1)
        x = struct.unpack_from("b",line)[0]
        return x

    def getInt(self):
        '''
        read int from current place in file
        '''
        line = self.f.read(4)
        x = struct.unpack_from("i",line)[0]
        return x
    
    def getString(self):
        '''
        get string from current place in file
        '''
        line = self.f.read(1)        
        count = struct.unpack_from("b", line)[0]
        line = self.f.read(count)
        t = count * "s"
        x = "".join(list(struct.unpack(t, line)))
        return x

    def parseInt(self,buffer,position):
        '''
        parse int from buffer at given position
        '''
        y = buffer[position:position + 4]
        x = struct.unpack_from("i",y)[0]
        return x

    def parseFloat(self,buffer,position):
        '''
        parse float from input buffer at given position
        '''
        y = buffer[position:position + 4]
        x = struct.unpack_from("f",y)[0]
        return x

