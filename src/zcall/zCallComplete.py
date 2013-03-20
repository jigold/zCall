#! /usr/bin/env python

"""Standalone script to run the complete zcall process.

zcall steps:
1. Generate thresholds
2. Evaluate thresholds
3. Merge evaluations to find best threshold
4. Apply zcall with given threshold to no-calls in input data

Standalone script is somewhat inefficient.
Does not allow parallelization of evaluation or calling, or reuse of thresholds.
However it may be useful as a convenience script and illustration of the zcall method.
"""
import os, sys
try: 
    import argparse, json
    from calibration import ThresholdFinder, SampleEvaluator, MetricEvaluator
    from runZCall import SampleCaller
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)


class ZCallComplete:

    def __init__(self):
        pass



def main():

    zcall = ZCallComplete()
    zcall.run()

def parseArgs():
    description = "Standalone script to run the complete zcall process: Generate and evaluate thresholds, and apply zcall to no-calls in the input data."
    parser = argparse.ArgumentParser(description=description)
    configDefault = os.path.join(sys.path[0], '../etc/config.ini')
    configDefault = os.path.abspath(configDefault)
    parser.add_argument('--egt', required=True, metavar="PATH", 
                        help="Path to .egt input file.")
    parser.add_argument('--config', metavar="PATH", default=configDefault,
                        help="Path to .ini config file. Default = etc/config.ini")
    parser.add_argument('--out', metavar="DIR", default=".",
                        help="Directory for output; defaults to current working directory.")
    parser.add_argument('--zstart', metavar="INT", default=7, type=int,
                    help='Starting z score. Default = %(default)s')
    parser.add_argument('--ztotal', metavar="INT", default=1, type=int,
                        help='Total number of integer z scores to generate. Default = %(default)s')
    
    parser.add_argument('--samples', required=True, metavar="PATH", 
                        help="Path to .json file containing sample URIs (unique identifiers), gender codes, and .gtc data paths")
    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Print status information to standard output")
    

    args = vars(parser.parse_args())
    return args

if __name__ == "__main__":
    main()
