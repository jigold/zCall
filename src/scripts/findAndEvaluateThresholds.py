#! /usr/bin/env python

# Iain Bancarz, ib5@sanger.ac.uk, January 2013

# Find thresholds for zcall and evaluate concordance on a list of GTC files
# Combines functionality of findThresholds.py and calibrationFromEGT.py

import os, sys, time
try: 
    import argparse     # optparse is deprecated, using argparse instead
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)
from evaluateConcordance import evaluator
from calibrationFromEGT import calibration

class multiEvaluator:

    def __init__(self, egt, bpm, configPath):
        self.egt = os.path.abspath(egt)
        self.bpm = os.path.abspath(bpm)
        self.cal = calibration(configPath)

    def findAndEvaluate(self, gtcList, zStart, zTotal, outDir, outName, 
                        verbose=True, force=False):
        z = zStart
        allResults = []
        for i in range(zTotal):
            threshPath = self.cal.run(self.egt, z, outDir, verbose, force)
            eva = evaluator(threshPath, self.bpm, self.egt)
            (includedSNPs, totalSNPs, results) = \
                eva.findMultipleConcordances(gtcList, verbose)
            for result in results: result.append(z)
            allResults.append(results)
            z += 1
        # includedSNPs, totalSNPs do not depend on GTC or z score
        (bestZ, bestZType) = self.findBestZ(allResults, verbose)
        outPath = os.path.join(outDir, outName)
        self.writeResults(outPath, includedSNPs, totalSNPs, 
                          bestZ, bestZType, allResults)

    def findBestZ(self, allResults, verbose=True):
        # find 'best' zscore from multiple GTC files and thresholds
        # defined as the smallest z s.t. mean concordance > mean gain
        # if none exists, return z with minimum of (gain - concordance)
        concords = {}
        gains = {}
        for results in allResults: # for each zscore
            for result in results: # for each gtc file
                [inPath, concordance, gain, z] = result
                try: concords[z].append(concordance)
                except KeyError: concords[z] = [concordance,]
                try: gains[z].append(gain)
                except KeyError: gains[z] = [gain,]
        concordanceGreaterThanGain = []
        gainMinusConcord = {}
        for z in concords.keys():
            cMean = sum(concords[z])/len(concords[z])
            gMean = sum(gains[z])/len(gains[z])
            if cMean > gMean: concordanceGreaterThanGain.append(z)
            gainMinusConcord[z] = gMean - cMean
        bestType = None
        if len(concordanceGreaterThanGain)>0: 
            bestType = 0
            best = min(concordanceGreaterThanGain)
        else: 
            bestType = 1
            leastDiff = min(gainMinusConcord.values())
            for z in gainMinusConcord.keys():
                if gainMinusConcord[z]==leastDiff:
                    best = z
                    break
        if verbose:
            print "BEST_Z", best
            print "BEST_Z_TYPE", bestType
        return (best, bestType)

    def writeResults(self, outPath, includedSNPs, totalSNPs, 
                     bestZ, bestZType, allResults, digits=3):
        includeRate = round(float(includedSNPs)/totalSNPs, digits)
        out = open(outPath, 'w')
        out.write("# EGT "+self.egt+"\n")
        out.write("# BPM "+self.bpm+"\n")
        out.write("# SNP_INCLUDE_RATE "+str(includeRate)+"\n")
        out.write("# BEST_Z "+str(bestZ)+"\n")
        out.write("# BEST_Z_TYPE "+str(bestZType)+"\n")
        headers = ['# Input', 'Concordance', 'Gain', 'Zscore']
        out.write("\t".join(headers)+"\n")
        for results in allResults:
            for result in results:
                [inPath, concord, gain, z] = result
                concord = round(concord, digits)
                gain = round(gain, digits)
                words = []
                for term in [inPath, concord, gain, z]: words.append(str(term))
                out.write("\t".join(words)+"\n")
        out.close()
        

def main():
    start = time.time()
    description = "Evaluate concordance and gain on multiple GTC files,"+\
        " for multiple Z scores.  Concordance and gain are respectively"+\
        " defined as calls on previous caller which match under zCall, and"+\
        " no-calls on previous caller which are called by zCall."+\
        " Also records a 'best' Z score, defined as the least z s.t."+\
        " mean concordance > mean gain (type 0); or if no such z exists, the "+\
        " z with the least value of mean gain - mean concordance (type 1)."
    configDefault = os.path.join(sys.path[0], '../etc/config.ini')
    configDefault = os.path.abspath(configDefault)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--bpm', required=True, metavar="PATH", 
                        help="BPM .csv manifest file")
    parser.add_argument('--egt', required=True, metavar="PATH", 
                        help="EGT input file")
    parser.add_argument('--gtc_list', required=True, metavar="PATH", 
                        help="File listing GTC input paths")
    parser.add_argument('--config', metavar="PATH", default=configDefault,
                        help="Path to .ini config file. Default = <installation directory>/etc/config.ini")
    parser.add_argument('--zstart', metavar="INT", default=7, type=int,
                        help='Starting z score. Default = %(default)s')
    parser.add_argument('--ztotal', metavar="INT", default=1, type=int,
                        help='Total number of integer z scores to generate. Default = %(default)s')
    parser.add_argument('--outdir', metavar="DIR", default=".",
                        help="Directory for output; defaults to current working directory. Threshold filename(s) will be of the form <prefix>_z<zscore>_thresholds.txt, for an input file of the form <prefix>.egt. Existing threshold files will not be overwritten unless --force is in effect.")
    parser.add_argument('--outname', metavar="STRING", required=True,
                        help='Output filename for concordance by zscore')
    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Print status information to standard output")
    parser.add_argument('--force', action='store_true', default=False,
                        help="Force overwrite of existing threshold files (if any)")
    args = vars(parser.parse_args())
    egt = args['egt']
    bpm = args['bpm']
    gtcList = args['gtc_list']
    outDir = args['outdir']
    outName = args['outname']
    config = args['config']
    zStart = args['zstart']
    zTotal = args['ztotal']
    verbose = args['verbose']
    force = args['force']
    multiEval = multiEvaluator(egt, bpm, config)
    multiEval.findAndEvaluate(gtcList, zStart, zTotal, 
                              outDir, outName, verbose, force)
    duration = time.time() - start
    if verbose: print "Finished.  Elapsed time: "+str(round(duration,0))+" s"

if __name__ == "__main__":
    main()
