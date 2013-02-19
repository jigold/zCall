#! /usr/bin/env python

"""Evaluate concordance and gain for a single GTC file and Z score.

Author:  Iain Bancarz, ib5@sanger.ac.uk, January 2013
"""

import os, re
from calibration import MetricFinder
from GTC import *
try: 
    import argparse     # optparse is deprecated, using argparse instead
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)

class evaluator:
    """Class to evaluate single GTC file on multiple z scores

    Inputs:
    Path to .json file with z scores and corresponding threshold paths
    Path to GTC input file

    Output:
    File with concordance and gain for each z score"""

    def __init__(self, thresholdsPath, gtcPath):

        self.gtc = GTC(gtcPath)

    def evaluate(self):
        pass



def main():
    """Method to run as script from command line.  Run with --help for usage."""
    description = "Evaluate concordance/gain for given zCall thresholds on a single GTC file."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--thresholds', required=True, metavar="PATH", 
                        help="Path to .json file containing thresholds paths")
    parser.add_argument('--bpm', required=True, metavar="PATH", 
                        help="BPM .csv manifest file")
    parser.add_argument('--egt', required=True, metavar="PATH", 
                        help="EGT input file")
    parser.add_argument('--gtc', required=True, metavar="PATH", 
                        help="Path to .gtc input file")
    parser.add_argument('--out', required=True, metavar="PATH", 
                        help="Output file")
    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Print status information to standard output")
    args = vars(parser.parse_args())
    inputKeys = ['thresholds', 'bpm', 'egt', 'gtc']
    for key in inputKeys:
        if not os.access(args[key], os.R_OK):
            raise OSError("Cannot read input path \""+args[key]+"\"")
        else:
            args[key] = os.path.abspath(args[key])
    (dirName, fileName) = os.path.split(os.path.abspath(args['out']))
    if fileName=='' or not os.access(dirName, os.R_OK):
        raise OSError("Invalid output path \""+args['out']+"\"")
    mf = MetricFinder(args['thresholds'], args['bpm'], args['egt'])
    mf.evaluate(args['gtc'], args['out'], args['verbose'])



if __name__ == "__main__":
    main()
