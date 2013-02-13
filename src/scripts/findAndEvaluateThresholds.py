#! /usr/bin/env python

"""Find thresholds for zcall and evaluate for a list of GTC files.

Combines functionality of evaluateConcordance.py and calibrationFromEGT.py.
Author: Iain Bancarz, ib5@sanger.ac.uk, January 2013
"""

import os, sys, time
from calibration import ZScoreEvaluator
try: 
    import argparse     # optparse is deprecated, using argparse instead
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)

def main():
    """Method to run as script from command line.  Run with --help for usage."""
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
    eva = ZScoreEvaluator(egt, bpm, config)
    eva.findAndEvaluate(gtcList, zStart, zTotal, 
                        outDir, outName, verbose, force)
    duration = time.time() - start
    if verbose: print "Finished.  Elapsed time: "+str(round(duration,0))+" s"

if __name__ == "__main__":
    main()
