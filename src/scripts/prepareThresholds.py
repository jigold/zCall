#! /usr/bin/env python

""" Find thresholds for given .egt file and Z score(s)
Combines findMeanSD.py, findBetas.r, findThresholds.py from original zCall

Iain Bancarz, ib5@sanger.ac.uk
January 2013
"""

import os, sys
try: 
    import argparse, json
    from calibration import ThresholdFinder
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)

"""
Calibration procedure:
1. Run findMeanSD.py on given EGT file.  Outputs mean_sd.txt
2. Run findBetas.r on output from (1). Outputs betas.txt
3. Run findThresholds.py on EGT file and betas.txt, with given Z score(s).
Outputs from (1) and (2) are written to a temporary directory, deleted on exit.

Recommended default Z score = 7.  Suggested range of alternatives = 3 to 15.
"""

def main():
    # 'main' method to run script from command line
    args = validate_args()
    egt = os.path.abspath(args['egt'])
    out = os.path.abspath(args['out'])
    indexPath = os.path.join(out, args['index_name'])
    z = args['zstart']
    tf = ThresholdFinder(os.path.abspath(args['config']))
    thresholdPaths = {}
    for i in range(args['ztotal']):
        thresholdPath = tf.run(egt, z, out, args['verbose'], args['force'])
        thresholdPaths[str(z)] = thresholdPath # .json requires string as key
        z += 1
    index = open(indexPath, 'w')
    index.write(json.dumps(thresholdPaths))
    index.close()
    
def validate_args():
    # parse command-line arguments and return dictionary of params
    description = "Generates threshold files for use with the zCall genotype caller.  Inputs are an .egt file and one or more Z score values.  The .egt file is a proprietary Illumina binary file format, containing typical means and standard deviations for intensity clusters.  An .egt file is supplied by Illumina for its own genotyping chips, or it may be generated using the GenomeStudio software for custom probe sets."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--egt', required=True, metavar="PATH", 
                        help="Path to .egt input file.")
    parser.add_argument('--out', metavar="DIR", default=".",
                        help="Directory for output; defaults to current working directory.  Filename(s) will be of the form <prefix>_z<zscore>_thresholds.txt, for an input file of the form <prefix>.egt")
    configDefault = os.path.join(sys.path[0], '../etc/config.ini')
    configDefault = os.path.abspath(configDefault)
    parser.add_argument('--config', metavar="PATH", default=configDefault,
                        help="Path to .ini config file. Default = etc/config.ini")
    parser.add_argument('--zstart', metavar="INT", default=7, type=int,
                    help='Starting z score. Default = %(default)s')
    parser.add_argument('--ztotal', metavar="INT", default=1, type=int,
                        help='Total number of integer z scores to generate. Default = %(default)s')
    parser.add_argument('--index_name', metavar="STRING", 
                        default="threshold_index.json",
                        help='Name for .json index file with paths to thresholds.txt output, written to output directory')
    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Print status information to standard output")
    parser.add_argument('--force', action='store_true', default=False,
                        help="Force overwrite of existing threshold files (if any)")
    args = vars(parser.parse_args())
    # validate arguments
    egt = args['egt']
    out = args['out']
    config = args['config']
    if not os.access(egt, os.R_OK):
        raise OSError("Cannot read .egt input path \""+egt+"\"")
    if not os.path.exists(out):
        raise OSError("Output path \""+out+"\" does not exist.")
    elif not os.path.isdir(out):
        raise OSError("Output path \""+out+"\" is not a directory.")
    elif not os.access(out, os.W_OK):
        raise OSError("Cannot write to output directory \""+out+"\"")
    if not os.access(config, os.R_OK):
        raise OSError("Cannot read config path \""+config+"\"")
    if args['ztotal']<1 or args['zstart']<1:
        raise ValueError("Invalid zstart or ztotal option.")
    return args

if __name__ == "__main__":
    main()
