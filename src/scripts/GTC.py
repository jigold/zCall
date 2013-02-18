#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# April 5th, 2012

# The Illumina provided Code was provided as-is and with no warranty as to performance and no warranty against it infringing any other party's intellectual property rights.

import struct
import math

class GTC:
    ''' Class for parsing GTC file'''

    def __init__(self, filename, bpmNormIDs):
        ''' Init function for class. Input is a .gtc file '''
        self.f = open(filename, 'rb') # open file handler for binary file
        self.BPMnormIDs = bpmNormIDs # list with norm ID for each snp
        
        self.TOC = self.parseTOC() # parse table of contents to get location of other information
        self.sampleName = self.parseString(10) # parse sample name
        self.samplePlate = self.parseString(11) # parse sample plate
        self.sampleWell = self.parseString(12) # parse sample well
        self.clusterFile = self.parseString(100) # parse what cluster file was used
        self.snpManifest = self.parseString(101) # parse what snp manifest was used
        self.imagingDate = self.parseString(200) # parse imaging date
        self.autoCallDate = self.parseString(201) # parse autocall date
        self.autoCallVersion = self.parseString(300) # parse autocall version
        self.rawXintensities = self.extractIntensities(1000) # parse raw x intensities into python list object
        self.rawYintensities = self.extractIntensities(1001) # parse raw y intensities into python list object
        self.normalizationTransformations = self.extractNormalizationTransformations(400) # parse normalization transformation arrays into python dictionary where key is the order they appeared in the gtc file and the value is a dictionary with keys offset_x, offset_y,scale_x, scale_y, shear, theta and values are floats        
        self.genotypes = self.extractGenotypes(1002) # parse genotypes (0,1,2,3) into python list object
        self.baseCalls = self.extractBaseCalls(1003) # parse basecalls (AT,TT,AT,--) into python list object

        self.normXintensities, self.normYintensities = self.normalizeIntensities()

            
    def parseTOC(self):
        '''Parse Table of Contents of GTC file
        No input
        Output is a dictionary where the ID for that entry is the key and the value is the offset for that variable in the GTC file
        '''
        self.f.seek(4,0)
        line = self.f.read(4)

        count = struct.unpack("i",line)[0] 
        TOC = {}

        for i in range(count):            
            line = self.f.read(2)
            id = struct.unpack("h",line)
            line = self.f.read(4)            
            offset = struct.unpack("I",line)
            TOC[id[0]] = offset[0]

        return TOC


    def parseString(self,id):
        '''
        Extract a string variable from GTC file such as SampleName.
        Input is ID for that variable in TOC
        Output is a string
        '''        
        offset = self.TOC[id]
        self.f.seek(offset,0)

        line = self.f.read(1)
        nbytes = struct.unpack("b",line)[0]
        line = self.f.read(nbytes)
        type = nbytes * "s"
        x = "".join(list(struct.unpack(type, line)))
        
        return x

    def extractIntensities(self, id):
        '''
        Extract intensity values (x or y depending on input ID).
        Input is ID for variable of interest in TOC
        Output is a list with integer intensity values in the order they were parsed
        '''        
        intensities = []
        offset = self.TOC[id]

        self.f.seek(offset,0)
        line = self.f.read(4)
        count = struct.unpack("i",line)[0]
        self.numSNPs = count
        
        for i in range(count):
            line = self.f.read(2)
            y = struct.unpack("H",line)
            intensities.append(y[0])

        return intensities

    def extractNormalizationTransformations(self, id):
        '''
        Extract normalization transformation arrays
        Input is ID for Normalization Transformations in TOC.
        Output is dictionary where keys are the order xForm array appears in gtc file (ex: 1,2,3...).
        The values of the dictionary are another dictionary
        where the keys are shear, offset_x, offset_y, theta, scale_x, scale_y and the values are floats
        '''
        normTransforms = {}
        offset = self.TOC[id]

        self.normIDlist = list(set(self.BPMnormIDs)) # ordered list of unique normIDs
        self.normIDlist.sort()
        
        self.f.seek(offset,0)
        line = self.f.read(4)
        count = struct.unpack("i",line)[0]

        for i in range(count):
            line = self.f.read(4)
            line = self.f.read(48)
            x = struct.unpack("<12f", line)
            normTransforms[self.normIDlist[i]] = {"offset_x":x[0],"offset_y":x[1],"scale_x":x[2],"scale_y":x[3],"shear":x[4],"theta":x[5]}

        return normTransforms
    
    def extractBaseCalls(self, id):
        '''
        Extract base calls.
        Input is id for BaseCalls in TOC
        Output is a list with one basecall for each SNP (ex: AT, GT,AA...)
        '''
        baseCalls = []
        offset = self.TOC[id]

        self.f.seek(offset,0)
        line = self.f.read(4)
        count = struct.unpack("i",line)[0]

        for i in range(count):
            line = self.f.read(2)
            calls = struct.unpack("ss",line)
            baseCalls.append(calls[0] + calls[1])

        return baseCalls

    def extractGenotypes(self, id):
        '''
        Extract genotypes.
        Input is ID for Genotypes in TOC
        Output is a list with one genotype per SNP (0,1,2,3)
        '''
        genotypes = []
        offset = self.TOC[id]

        self.f.seek(offset,0)
        line = self.f.read(4)
        count = struct.unpack("i",line)[0]

        for i in range(count):
            line = self.f.read(1)
            gt = struct.unpack("b",line)
            genotypes.append(gt[0])
            
        return genotypes


    def getTotalSNPs(self):
        return self.numSNPs

    def normalizeIntensities(self):
        '''
        Use Normalization transformations to convert raw intensities to normalized intensities
        No Input
        Outputs are normalized x and y intensities in python lists
        '''
        normXIntensities = []
        normYIntensities = []
        
        for i in range(self.numSNPs):
            xraw = self.rawXintensities[i]
            yraw = self.rawYintensities[i]
            normID = self.BPMnormIDs[i]

            offset_x = self.normalizationTransformations[normID]["offset_x"]
            offset_y = self.normalizationTransformations[normID]["offset_y"]
            scale_x = self.normalizationTransformations[normID]["scale_x"]
            scale_y = self.normalizationTransformations[normID]["scale_y"]
            theta = self.normalizationTransformations[normID]["theta"]
            shear = self.normalizationTransformations[normID]["shear"]

            tempx = xraw - offset_x
            tempy = yraw - offset_y

            tempx2 = math.cos(theta) * tempx + math.sin(theta) * tempy
            tempy2 = -1 * math.sin(theta) * tempx + math.cos(theta) * tempy

            tempx3 = tempx2 - (shear * tempy2)
            tempy3 = tempy2

            xn = tempx3 / float(scale_x)
            yn = tempy3 / float(scale_y)

            if xn < 0:
                xn = 0.0
            if yn < 0:
                yn = 0.0

            normXIntensities.append(xn)
            normYIntensities.append(yn)

        return (normXIntensities, normYIntensities)            
