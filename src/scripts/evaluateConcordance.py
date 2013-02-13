#! /usr/bin/env python

"""Evaluate concordance and call gain for a list of GTC files.

Author:  Iain Bancarz, ib5@sanger.ac.uk, January 2013
"""

# TODO exclude failing samples from concordance evaluation
# either assess CR from GTC files, or use qc_results.json

import os
from calibration import ConcordanceGainFinder
try: 
    import argparse     # optparse is deprecated, using argparse instead
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)

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
    cgf = ConcordanceGainFinder(args['thresholds'], args['bpm'], args['egt'])
    cgf.evaluate(args['gtc_list'], args['out'], args['verbose'])

if __name__ == "__main__":
    main()
