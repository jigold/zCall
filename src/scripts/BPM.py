#! /usr/bin/python
import sys

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# April 5th, 2012

# The Illumina provided Code was provided as-is and with no warranty as to performance and no warranty against it infringing any other party's intellectual property rights.

class BPM:
    ''' Python class to parse a .bpm.csv file '''
    def __init__(self, bpmFile):
        self.names = []
        self.chr = []
        self.pos = []
        self.normID = []
        self.A = []
        self.B = []
        
        for line in open(bpmFile, 'r'):
            line = line.replace("\n", "")
            line = line.replace("\r", "")

            fields = line.split(",")

            if line.find("Chromosome") != -1:
                continue

            else:
                self.names.append(fields[1])
                self.chr.append(fields[2])
                self.pos.append(fields[3])
                alleles = fields[5].replace("[", "")
                alleles = alleles.replace("]", "")
                alleles = alleles.split("/")
                self.A.append(alleles[0]) # allele A
                self.B.append(alleles[1]) # allele B
                self.normID.append(int(fields[8])) # normalization ID for that snp
        self.numSNPs = len(self.names)

    def getTotalSNPs(self):
        return self.numSNPs
