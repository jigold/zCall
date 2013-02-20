#! /usr/bin/env python

"""Combine concordance and gain metric results for multiple samples and zscores.

Use to find 'best' z score for calling.

Author:  Iain Bancarz, ib5@sanger.ac.uk, January 2013
"""

import json, os, re
from calibration import MetricEvaluator
try: 
    import argparse    
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)

def main():
    """Method to run as script from command line.  Run with --help for usage."""
    description = "Evaluate concordance/gain results for multiple thresholds and samples; find the 'best' z score for subsequent use of zCall."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--metrics', required=True, metavar="PATH", 
                        help="Path to .json file containing paths to .json metrics files")
    parser.add_argument('--thresholds', required=True, metavar="PATH", 
                        help="Path to .json file containing threshold .txt paths indexed by z score")
    parser.add_argument('--out', required=False, metavar="PATH", 
                        help="Directory for .json output", default=".")
    args = vars(parser.parse_args())
    
    MetricEvaluator().writeBest(args['metrics'], args['thresholds'], 
                                args['out'])

if __name__ == "__main__":
    main()
