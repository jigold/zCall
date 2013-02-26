#! /usr/bin/env python

import struct
from GTC import *
from utilities import CallingBase, ThresholdContainer

class SampleCaller(CallingBase):

    """Class to run zCall on given GTC files

    initialize with: bpmPath, egtPath, threshPath
    Write output in .bed format; use merge_bed in genotyping pipeline workflows
    """

    def callsToBinary(self, calls, outPath):
        """Translate genotype calls for one or more samples to Plink binary

        4 genotype calls are packed into one byte of output
        Returns a list of struct.pack strings corresponding to output bytes"""
        callTotal = len(calls)
        if callTotal % self.snpTotal != 0:
            raise ValueError("Number of calls is not a multiple of SNP total!")
        elif callTotal % 4 != 0:
            # if not an integer number of bytes, pad with no calls
            calls.extend([0]*(callTotal % 4)) 
        output = []
        i = 0
        while i < callTotal:
            byte = struct.pack(self.callsToByte(calls[i:i+4]))
            output.append(byte)
            i += 4
        return output

    def callsToByte(self, calls):
        """Convert list of 4 calls to an integer in Plink binary format 

        Create byte string of the form '01001101', convert to integer"""
        if len(calls) != 4:
            raise ValueError("Must have exactly 4 calls for byte conversion!")
        byte = []
        for call in calls:
            if call==1: bcall = '00' # major homozygote, 'AA'
            elif call==2: bcall = '01' # heterozygote, 'AB'
            elif call==3: bcall = '11' # minor homozygote, 'BB'
            else: bcall = '10' # missing genotype
            byte.append(bcall)
        byteString = ''.join(byte)
        byteString[::-1] # reverse order of string characters
        return int(byteString, 2)

    def makeCalls(self, gtc):
        """Apply zCall to 'no calls' from GTC input

        Return genotypes in numeric format
        0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        """
        calls = []
        for i in range(self.snpTotal):
            origCall = gtc.genotypes[i]
            if origCall == 0 \
                    and self.thresholds.getX(i) != "NA" \
                    and self.thresholds.getY(i) != "NA":
                calls.append(self.call(gtc, i))
            else:
                calls.append(origCall)
        return calls

    def run(self, sampleJsonPath, outPath):
        """Apply zCall to GTC files and write Plink .bed output"""
        gtcPaths = self.readSampleJson(sampleJsonPath)
        calls = []
        for gtcPath in gtcPaths:
            gtc = GTC(gtcPath, self.bpm.normID)
            calls.extend(self.makeCalls(gtc))
        self.writeBed(calls, outPath)
            
    def writeBed(self, calls, outPath):
        """Write output for one or more samples in Plink .bed format

        Input: List of lists of genotype call codes, and output path
        Output file:  First 2 bytes are Plink magic number
        3rd byte is flag for an individual-major file
        Subsequent bytes represent genotype calls
        """
        header = [0b01101100, 0b00011011, 0b00000000]
        output = []
        for byte in header: output.append(struct.pack('B', byte)
        output.extend(self.callsToBinary(calls))
        out = open(outPath, 'w')
        for byte in output: out.write(byte)
        out.close()
