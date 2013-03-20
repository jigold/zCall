#! /usr/bin/env python

import os, struct
from GTC import *
try: 
    import argparse, json
    from calibration import ThresholdFinder
    from utilities import CallingBase, ThresholdContainer
    from plink import PlinkHandler
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)


class SampleCaller(CallingBase):

    """Class to run zCall on given GTC files

    initialize with: bpmPath, egtPath, threshPath
    Write output in .bed format; use merge_bed in genotyping pipeline workflows
    Construct Plink .bim manifest file from BPM object
    Optionally, use sample JSON file to construct Plink .fam file
    """

    def makeCalls(self, gtc, null=False):
        """Apply zCall to 'no calls' from GTC input

        Return genotypes in numeric format
        0 - "No Call", 1 - AA, 2 - AB, 3 - BB

        If null, output the unchanged original calls (use for testing)
        """
        calls = []
        for i in range(self.snpTotal):
            origCall = gtc.genotypes[i]
            if null==False \
                    and origCall == 0 \
                    and self.thresholds.getX(i) != "NA" \
                    and self.thresholds.getY(i) != "NA":
                calls.append(self.call(gtc, i))
            else:
                calls.append(origCall)
        return calls

    def run(self, samplesPath, outStem, verbose=False, null=False, 
            reorder=True):
        """Apply zCall to GTC files and write Plink .bed output"""
        gtcPaths = self.readSampleJson(samplesPath)
        calls = []
        ph = PlinkHandler(self.bpm)
        for gtcPath in gtcPaths:
            if verbose: print "Calling GTC file", gtcPath
            gtc = GTC(gtcPath, self.bpm.normID)
            calls.extend(ph.callsToBinary(self.makeCalls(gtc, null), reorder))
        ph.writeBed(calls, outStem+'.bed', verbose)
        ph.writeBim(outStem+'.bim')
        ph.writeFam(samplesPath, outStem+'.fam')


def main():
    """Method to run as script from command line"""
    description = "Apply zCall to no-calls with given threshold and samples."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--thresholds', required=True, metavar="PATH", 
                        help="Path to zCall thresholds.txt file")
    parser.add_argument('--bpm', required=True, metavar="PATH", 
                        help="BPM .csv manifest file")
    parser.add_argument('--egt', required=True, metavar="PATH", 
                        help="EGT input file")    
    parser.add_argument('--samples', required=True, metavar="PATH", 
                        help="Path to .json file containing sample URI's, genders, and .gtc input paths")
    parser.add_argument('--out', required=True, metavar="PATH", 
                        help="Path stem for Plink binary output, without .bed, .bim, .fam suffix")
    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Print status information to standard output")
    parser.add_argument('--null', action='store_true', default=False,
                        help="Do not apply zcall. Instead output GTC calls unchanged to an individual-major Plink binary file. Used for testing.")
    args = vars(parser.parse_args())
    inputKeys = ['thresholds', 'bpm', 'egt', 'samples']
    for key in inputKeys:
        if not os.access(args[key], os.R_OK):
            raise OSError("Cannot read input path \""+args[key]+"\"")
        else:
            args[key] = os.path.abspath(args[key])
    (dirName, fileName) = os.path.split(os.path.abspath(args['out']))
    if not os.access(dirName, os.R_OK) or not os.path.isdir(dirName):
        raise OSError("Invalid output path \""+args['out']+"\"")
    else:
        args['out'] = os.path.join(os.path.abspath(dirName), fileName)
    if args['null']:
        print "WARNING: Null option in effect, input calls will not be changed"
    caller = SampleCaller(args['bpm'], args['egt'], args['thresholds'])
    caller.run(args['samples'], args['out'], args['verbose'], args['null'])

if __name__ == "__main__":
    main()
