#! /usr/bin/env python

"""Evaluate concordance and call gain for a list of GTC files.

Author:  Iain Bancarz, ib5@sanger.ac.uk, January 2013
"""

import os, re
try: 
    import argparse     # optparse is deprecated, using argparse instead
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)
from GTC import *
from BPM import *
from EGT import *
from zCallBase import zCallBase

class evaluator(zCallBase):
    """Class to evaluate multiple GTC files for given thresholds"""

    def concordanceRate(self, counts):
        """Find concordance rate between original and new call counts

        Ignores SNPs where original was a 'no call'"""
        [match, total] = [0,0]
        for i in range(1,4):
            for j in range(0,4):
                count = counts[(i,j)]
                total += count
                if i==j: match += count
        concord = float(match)/float(total)
        return concord

    def countCallTypes(self, gtc):
        """based on method in sampleConcordance.py

        call codes: 0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        also return rate of inclusion for SNPs"""
        numSNPs = gtc.numSNPs
        included = 0
        counts = {}
        for i in range(4):
            for j in range(4): counts[(i,j)] = 0
        for i in range(numSNPs):
            nAA = self.egt.nAA[i]
            nBB = self.egt.nBB[i]
            nAB = self.egt.nAB[i]
            if not self.includeSNP(i, nAA, nBB, nAB): continue
            origCall = self.normalizeCall(gtc.genotypes[i], nAA, nBB)
            newCall = self.normalizeCall(self.call(gtc, i), nAA, nBB)
            counts[(origCall, newCall)] += 1
            included += 1
        return (included, numSNPs, counts)

    def evaluate(self, inPath, outPath, verbose=True):
        """inPath to file listing multiple GTC input paths, one per line"""
        self.writeConcordances(inPath, outPath, verbose)

    def findConcordance(self, gtcPath):
        gtc = GTC(gtcPath, self.bpm.normID)
        (includedSNPs, totalSNPs, counts) = self.countCallTypes(gtc)
        concord = self.concordanceRate(counts)
        gain = self.gainRate(counts)
        return (includedSNPs, totalSNPs, concord, gain)

    def findMultipleConcordances(self, gtcListPath, verbose=True):
        gtcPaths = []
        for line in open(gtcListPath).readlines():
            # read input, ignoring comments and blank lines
            if re.match('#', line): continue
            line = line.strip()
            if line!='': gtcPaths.append(line)
        results = []
        gtcTotal = len(gtcPaths)
        # snp inclusion is the same for all GTC files (depends only on EGT)
        for i in range(gtcTotal):
            if verbose: print "Evaluating GTC path %s of %s" % (i+1, gtcTotal)
            (includedSNPs, totalSNPs, concord, gain) = \
                self.findConcordance(gtcPaths[i])
            results.append([gtcPaths[i], concord, gain])
        return (includedSNPs, totalSNPs, results)

    def gainRate(self, counts):
        """Find rate of call gain

        Defined as no calls in original which are called by zcall"""
        [gain, total] = [0,0]
        for i in range(4):
            count = counts[(0,i)] # no call in original GTC
            total += count
            if i!=0: gain += count
        gainRate = float(gain)/float(total)
        return gainRate

    def includeSNP(self, i, nAA, nBB, nAB):
        """Should ith SNP be included in concordance calculation?

        Require autosomal SNP with MAF>=5%
        Want at least 10 points in each homozygote cluster
        Also exclude SNPs without defined zcall thresholds"""
        include = True
        chrom = self.bpm.chr[i]
        maf = self.findMAF(nAA, nBB, nAB)
        if maf < 0.05 or chrom == "X" or chrom == "Y" or nAA < 10 or nBB < 10 \
                or self.thresholdsX[i]=="NA" or self.thresholdsY[i]=="NA":
            include = False
        return include        

    def writeConcordances(self, gtcListPath, outPath, verbose=True, digits=3):
        """Find concordances/gains and write results to file"""
        (includedSNPs, totalSNPs, results) = \
            self.findMultipleConcordances(gtcListPath, verbose)
        includeRate = float(includedSNPs)/totalSNPs
        headers = [
            '# evaluateConcordance.py results',
            '# BPM '+self.bpmPath,
            '# EGT '+self.egtPath,
            '# THRESHOLDS '+self.threshPath,
            '# INCLUDED_SNP '+str(includedSNPs),
            '# TOTAL_SNP '+str(totalSNPs),
            '# INCLUDE_RATE_SNP '+str(round(includeRate, digits)),
            '# [Input] [Concordance on original calls] [Gain]'
            ]
        out = open(outPath, 'w')
        for header in headers: 
            out.write(header+"\n")
        for result in results:
            [inPath, concord, gain] = result
            concord = round(concord, digits)
            gain = round(gain, digits)
            out.write("%s\t%s\n" % (inPath, concord, gain))
        out.close()
        if verbose: print "Finished.\n"

def main():
    """Method to run as script from command line.  Run with --help for usage."""
    description = "Evaluate concordance on multiple GTC files."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--thresholds', required=True, metavar="PATH", 
                        help="zCall thresholds .txt file")
    parser.add_argument('--bpm', required=True, metavar="PATH", 
                        help="BPM .csv manifest file")
    parser.add_argument('--egt', required=True, metavar="PATH", 
                        help="EGT input file")
    parser.add_argument('--gtc_list', required=True, metavar="PATH", 
                        help="File listing GTC input paths")
    parser.add_argument('--out', required=True, metavar="PATH", 
                        help="Output file")
    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Print status information to standard output")
    args = vars(parser.parse_args())
    inputKeys = ['thresholds', 'bpm', 'egt', 'gtc_list']
    for key in inputKeys:
        if not os.access(args[key], os.R_OK):
            raise OSError("Cannot read input path \""+args[key]+"\"")
        else:
            args[key] = os.path.abspath(args[key])
    (dirName, fileName) = os.path.split(os.path.abspath(args['out']))
    if fileName=='' or not os.access(dirName, os.R_OK):
        raise OSError("Invalid output path \""+args['out']+"\"")
    myEvaluator = evaluator(args['thresholds'], args['bpm'], args['egt'])
    myEvaluator.evaluate(args['gtc_list'], args['out'], args['verbose'])

if __name__ == "__main__":
    main()
