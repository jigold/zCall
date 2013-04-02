#! /usr/bin/env python

"""Process Plink format genotyping data

See http://pngu.mgh.harvard.edu/~purcell/plink/
"""

import json, struct
from BPM import *

class PlinkHandler:
    """Class to handle Plink format data"""

    def __init__(self, bpm):
        """Initialize with a BPM object"""
        self.bpm = bpm
        self.snpTotal = self.bpm.getTotalSNPs()
        self.sortMap = self.snpSortMap()

    def callsToBinary(self, calls, reorder=True):
        """Translate genotype calls for one sample to Plink binary

        4 genotype calls are packed into one byte of output
        Returns a list of struct.pack strings corresponding to output bytes"""
        if len(calls) != self.snpTotal:
            msg = "Number of calls %s is not equal to SNP total %s" % \
                (len(calls), self.snpTotal)
            raise ValueError(msg)
        if reorder:
            sortedCalls = [None]*self.snpTotal
            for i in range(self.snpTotal):
                sortedCalls[self.sortMap[i]] = calls[i] 
            calls = sortedCalls
        if self.snpTotal % 4 != 0:
            # if not an integer number of bytes, pad with no calls
            calls.extend([0]*(self.snpTotal % 4)) 
        output = []
        i = 0
        while i < self.snpTotal:
            byte = struct.pack('B', self.callsToByte(calls[i:i+4]))
            output.append(byte)
            i += 4
        return output

    def callsToByte(self, calls):
        """Convert list of 4 calls to an integer in Plink binary format 

        Create byte string of the form '01001101', convert to integer
        See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
        """
        if len(calls) != 4:
            raise ValueError("Must have exactly 4 calls for byte conversion!")
        byte = []
        for call in calls:
            if call==1: bcall = '00' # major homozygote, 'AA'
            elif call==2: bcall = '01' # heterozygote, 'AB'
            elif call==3: bcall = '11' # minor homozygote, 'BB'
            else: bcall = '10' # missing genotype, call=0
            byte.append(bcall)
        byteString = ''.join(byte)
        byteString = byteString[::-1] # reverse order of string characters
        return int(byteString, 2)

    def getSampleFields(self, sample):
        """Get 6 sample metadata fields for .fam or .ped file

        Fields are:
        - Family ID
        - Individual ID
        - Paternal ID
        - Maternal ID
        - Sex (1=male; 2=female; other=unknown)
        - Phenotype
        Conventionally, set family/individual IDs to sample URI
        Set sex if known
        Other values set to -9 as placeholder"""
        fields = ['-9']*6
        fields[0] = sample['uri']
        fields[1] = sample['uri']
        fields[4] = str(sample['gender_code'])
        return fields

    def numericChromosomes(self, chroms):
        """Convert to numeric chromosome IDs used by Plink"""
        for i in range(len(chroms)):
            if chroms[i]=='X': chroms[i] = 23
            elif chroms[i]=='Y': chroms[i] = 24
            elif chroms[i]=='XY': chroms[i] = 25
            elif chroms[i]=='MT': chroms[i] = 26
            else: chroms[i] = int(chroms[i])
        return chroms

    def parseBed(self, bed):
        """Parse a single byte from a .bed file"""
        parsed = bin(ord(bed))[2:]
        gap = 8 - len(parsed)
        if gap > 0: parsed = ''.join(['0']*gap)+parsed
        elif gap < 0: raise ValueError
        # parsed is now a string of the form '01101100'
        return parsed

    def readGenotypes(self, parsedBed):
        """Read a block of 4 genotypes from a parsed .bed string 

        Return in numeric format: 0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        """
        parsedBed = parsedBed[::-1] # reverse character order
        i = 0
        gtypes = []
        while i < len(parsedBed):
            pair = parsedBed[i:i+2]
            if pair=='00': gtype = 1
            elif pair=='01': gtype = 2
            elif pair=='11': gtype = 3
            elif pair=='10': gtype = 0
            else: raise ValueError("Invalid genotype string")
            gtypes.append(gtype)
            i += 2
        return gtypes

    def readBedFile(self, bedPath, sampleTotal):
        """Read genotype calls from a Plink .bed file

        Input should be in (default) SNP-major order
        May have extra 'null' calls, to have integer number of bytes per sample
        Need total samples to identify and strip off null padding (if any)
        """
        bed = open(bedPath).read()
        total = 0
        if ord(bed[0])!=108 or ord(bed[1])!=27:
            raise ValueError("Header does not start with Plink 'magic number'!")
        elif ord(bed[2])!=1:
            raise ValueError("Plink .bed file must be in SNP-major order.")
        gtypes = []
        calls = 0 # calls for current SNP; should be one per sample
        for i in range(3, len(bed)): # skip first 3 bytes 
            gtBlock = self.readGenotypes(self.parseBed(bed[i])) # 4 genotypes
            gtypes.extend(gtBlock)
            calls += len(gtBlock)
            if calls >= sampleTotal: # remove padding (if any), reset count
                while len(gtypes) % sampleTotal != 0: gtypes.pop()
                calls = 0
        return gtypes

    def snpSortMap(self):
        """Sort snps into (chromosome, position) order

        Ensures compatibility with sorted .bim files generated by Plink
        Return a map from original position to sorted position"""
        chroms = self.numericChromosomes(self.bpm.getChromosomes())
        pos = self.bpm.getPositions()
        coords = [None]*self.snpTotal
        for i in range(self.snpTotal):
            coords[i] = (chroms[i], int(pos[i]), i)
        coords.sort()
        sortMap = {}
        for i in range(self.snpTotal):
            [chrom, pos, orig] = coords[i]
            sortMap[orig] = i
        return sortMap

    def writeBed(self, binaryCalls, outPath, verbose=False):
        """Write output for one or more samples in Plink .bed format

        Input: List of call bytes, and output path
        Output file:  First 2 bytes are Plink magic number
        3rd byte is flag for an individual-major file
        Subsequent bytes represent genotype calls
        """
        header = [0b01101100, 0b00011011, 0b00000000]
        output = []
        for byte in header: output.append(struct.pack('B', byte))
        output.extend(binaryCalls)
        out = open(outPath, 'w')
        for byte in output: out.write(byte)
        out.close()
        if verbose: print len(output), "bytes written."

    def writeBim(self, outPath):
        """Write a Plink .bim file to accompany .bed output

        Similar to Plink .map format, except:
        - 2 additional columns for allele names (use A and B as dummy values)
        - Entries are *sorted* into (chromosome, position) order
        - Chromosomes are given numeric codes (including for X, Y, etc.)
        """
        unsorted = [None]*self.snpTotal
        for i in range(self.snpTotal):
             snp = self.bpm.names[i]
             chr = self.bpm.chr[i]
             pos = self.bpm.pos[i]
             alleleA = self.bpm.A[i]
             alleleB = self.bpm.B[i]
             out = [chr, snp, "0", pos, alleleA, alleleB]
             unsorted[i] = out
        # sort manifest entries
        out = [None]*self.snpTotal
        for i in range(self.snpTotal):
            out[self.sortMap[i]] = unsorted[i]
        # write to file
        outFile = open(outPath, 'w')
        for i in range(self.snpTotal):
            words = []
            for item in out[i]: words.append(str(item))
            outFile.write("\t".join(words)+"\n")
        outFile.close()

    def writeFam(self, sampleJson, outPath):
        """Write a Plink .fam file to accompany .bed output"""
        samples = json.loads(open(sampleJson).read())
        outLines = []
        for sample in samples:
            outLines.append(' '.join(self.getSampleFields(sample))+"\n")
        outFile = open(outPath, 'w')
        outFile.write("".join(outLines))
        outFile.close()

    def writeMap(self, outPath):
        """Write Plink .map format file"""
        outLines = []
        for i in range(self.snpTotal):
            snp = self.bpm.names[i]
            chr = self.bpm.chr[i]
            pos = self.bpm.pos[i]
            out = [str(chr), str(snp), "0", str(pos)]
            outLines.append("\t".join(out)+"\n")
        outFile = open(outPath, 'w')
        outFile.write("".join(outLines))
        outFile.close()

    def writePed(self, calls, sampleJson, outPath):
        """Write Plink .ped format file.

        Each line represents one sample.
        First 6 fields are same as for .fam file, see writeFam() method.  
        Subsequent fields are allele pairs, eg. 'G C' or 'A B'. 
        If manifestNames==True then get allele symbols from manifest,
        otherwise use A and B respectively for major and minor alleles."""
        samples = json.loads(open(sampleJson).read())
        outLines = []
        if len(calls) % self.snpTotal !=0:
            msg = "Number of calls %s is not a multiple of SNP total %s" % \
                (len(calls), self.snpTotal)
            raise ValueError(msg)
        for i in range(len(samples)):
            fields = self.getSampleFields(samples[i])
            start = i * self.snpTotal
            for j in range(self.snpTotal):
                alleleA = self.bpm.A[j]
                alleleB = self.bpm.B[j]
                call = calls[start+j]
                if call == 1: symbol = alleleA+' '+alleleA
                elif call == 2: symbol = alleleA+' '+alleleB
                elif call == 3: symbol = alleleB+' '+alleleB
                else: symbol = "0 0"
                fields.append(symbol)
            outLines.append("\t".join(fields)+"\n")
        outFile = open(outPath, 'w')
        outFile.write("".join(outLines))
        outFile.close()
