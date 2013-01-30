#! /usr/bin/env python

# Iain Bancarz, ib5@sanger.ac.uk, January 2013

# Evaluate concordance for a list of GTC files

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

    def concordanceRate(self, counts):
        # find concordance rate between original and new call counts
        # ignore SNPs where original was a "no call"
        [match, total] = [0,0]
        for i in range(1,4):
            for j in range(0,4):
                count = counts[(i,j)]
                total += 1
                if i==j: match += 1
        concord = float(match)/float(total)
        return concord

    def countCallTypes(self, gtc):
        # based on method in sampleConcordance.py
        # call codes: 0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        # also return rate of inclusion for SNPs
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
        includeRate = float(included)/numSNPs
        return (includeRate, included, counts)

    def evaluate(self, inPath, outPath, verbose=True):
        # input file lists multiple GTC input paths, one per line
        self.writeConcordances(gtcListPath, outPath, verbose)

    def findConcordance(self, gtcPath):
        gtc = GTC(gtcPath, self.bpm.normID)
        (includeRate, included, counts) = self.countCallTypes(gtc)
        concord = self.concordanceRate(counts)
        return (includeRate, included, concord)

    def findMultipleConcordances(self, gtcListPath, verbose=True, digits=3):
        gtcPaths = []
        for line in open(gtcListPath).readlines():
            # read input, ignoring comments and blank lines
            if re.match('#', line): continue
            line = line.strip()
            if line!='': gtcPaths.append(line)
        results = []
        total = len(gtcPaths)
        for i in range(total):
            if verbose: print "Evaluating GTC path %s of %s" % (i+1, total)
            (includeRate, included, concord) = self.findConcordance(gtcPaths[i])
            results.append([gtcPaths[i], round(includeRate, digits), included,
                            round(concord, digits)])
        return results

    def includeSNP(self, i, nAA, nBB, nAB):
        # should ith SNP be included in concordance calculation?
        # require autosomal SNP with MAF>=5%
        # want at least 10 points in each homozygote cluster
        # also exclude SNPs without defined zcall thresholds
        include = True
        chrom = self.bpm.chr[i]
        maf = self.findMAF(nAA, nBB, nAB)
        if maf < 0.05 or chrom == "X" or chrom == "Y" or nAA < 10 or nBB < 10 \
                or self.thresholdsX[i]=="NA" or self.thresholdsY[i]=="NA":
            include = False
        return include        

    def writeConcordances(self, gtcListPath, outPath, verbose=True):
        # evaluate thresholds and write results for various GTC paths
        results = self.findMultipleConcordances(gtcListPath, verbose)
        headers = [
            '# evaluateConcordance.py results',
            '# BPM '+self.bpmPath,
            '# EGT '+self.egtPath,
            '# THRESHOLDS '+self.threshPath,
            '# [Input] [Filter pass rate] [Filter pass total] '+\
                '[Concordance on original called SNPs]'
            ]
        out = open(outPath, 'w')
        for header in headers: 
            out.write(header+"\n")
        for result in results:
            out.write("%s\t%s\t%s\t%s\n" % result)
        out.close()
        if verbose: print "Finished.\n"

def main():

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
