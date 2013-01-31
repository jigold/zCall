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

    def __init__(self):
        pass

    def findAndEvaluate(self, egt, bpm, gtcList, configPath, 
                        zStart, zTotal, outDir, outName, verbose=True):
        z = zStart
        cal = calibration(configPath)
        allResults = []
        for i in range(zTotal):
            threshPath = cal.run(egt, z, outDir, verbose)
            eva = evaluator(threshPath, bpm, egt)
            results = eva.findMultipleConcordances(gtcList, verbose)
            for result in results: result.append(z)
            allResults.append(results)
            z += 1
        headers = ['# Input', 'Include_rate', 'Include_total', 'Concordance',
                   'Zscore']
        out = open(os.path.join(outDir, outName), 'w')
        out.write("\t".join(headers)+"\n")
        for results in allResults:
            for result in results:
                words = []
                for term in result: words.append(str(term))
                out.write("\t".join(words)+"\n")
        out.close()

def main():

    start = time.time()
    description = "Evaluate concordance on multiple GTC files,"+\
        " for multiple Z scores."
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
                        help="Path to .ini config file. Default = etc/config.ini")
    parser.add_argument('--zstart', metavar="INT", default=7, type=int,
                    help='Starting z score. Default = %(default)s')
    parser.add_argument('--ztotal', metavar="INT", default=1, type=int,
                        help='Total number of integer z scores to generate. Default = %(default)s')
    parser.add_argument('--outdir', metavar="DIR", default=".",
                        help="Directory for output; defaults to current working directory.  Threshold filename(s) will be of the form <prefix>_z<zscore>_thresholds.txt, for an input file of the form <prefix>.egt")
    parser.add_argument('--outname', metavar="STRING",
                        help='Output filename for concordance by zscore')
    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Print status information to standard output")
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

    multiEvaluator().findAndEvaluate(egt, bpm, gtcList, config, 
                                     zStart, zTotal, outDir, outName, verbose)

    duration = time.time() - start
    print "Finished.  Elapsed time: "+str(round(duration,0))+" s\n"

if __name__ == "__main__":
    main()
