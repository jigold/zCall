#! /usr/bin/env python

"""Standalone, self-contained script to run the complete zcall process.

zcall steps:
1. Generate thresholds
2. Evaluate thresholds
3. Merge evaluations to find best threshold
4. Apply zcall with given threshold to no-calls in input data

This standalone script does not allow parallelization of evaluation or calling, 
or reuse of thresholds. It is intended as a convenience script for small
datasets where parallelization is not required, and an illustration of the 
zcall method.
"""
import os, sys, time
try: 
    import argparse, json
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)
from calibration import ThresholdFinder, SampleEvaluator, MetricEvaluator
from runZCall import SampleCaller

class ZCallComplete:

    EVALUATION = 'threshold_evaluation.json'
    MERGED = 'merged_evaluation.json'
    LOG = 'zcall_log.json'

    def __init__(self, args):
        self.args = args

    def call(self, bpm, egt, thresholdPath, sJson, outDir, prefix, verbose):
        if verbose: print "Running zcall with thresholds", thresholdPath
        logPath = os.path.join(outDir, self.LOG)
        caller = SampleCaller(bpm, egt, thresholdPath)
        caller.run(sJson, outDir, prefix, logPath, verbose)

    def evaluate(self, bpm, egt, thresholds, gtc, start, end, outDir, verbose):
        """Evaluate threshold.txt files by concordance and gain metrics"""
        outPath = os.path.join(outDir, self.EVALUATION)
        eva = SampleEvaluator(bpm, egt)
        eva.run(thresholds, gtc, start, end, outPath, verbose)
        return outPath

    def merge(self, metricPath, thresholdJson, outDir, verbose):
        outPath = os.path.join(outDir, self.MERGED)
        eva = MetricEvaluator()
        results = eva.writeBest([metricPath,], thresholdJson, outPath, verbose)
        thresholdPath = results[eva.getBestThresholdKey()]
        return thresholdPath

    def prepare(self, zstart, ztotal, egtPath, outDir, config, verbose=False):
        """ Prepare threshold.txt files for given range of z scores"""
        tf = ThresholdFinder(config)
        return tf.runMultiple(zstart, ztotal, egtPath, outDir, verbose)
        

    def run(self):
        tJson = self.prepare(self.args['zstart'],
                             self.args['ztotal'],
                             self.args['egt'],
                             self.args['out'],
                             self.args['config'],
                             self.args['verbose'])
        mJson = self.evaluate(self.args['bpm'],
                              self.args['egt'],
                              tJson,
                              self.args['samples'],
                              self.args['gtc_start'],
                              self.args['gtc_end'],
                              self.args['out'],
                              self.args['verbose'])
        thresholdPath = self.merge(mJson, 
                                   tJson, 
                                   self.args['out'], 
                                   self.args['verbose'])
        self.call(self.args['bpm'],
                  self.args['egt'],
                  thresholdPath,
                  self.args['samples'],
                  self.args['out'],
                  self.args['plink'],
                  self.args['verbose'])
        

def main():
    args = parseArgs()
    start = time.time()
    ZCallComplete(args).run()
    if args['verbose']==True:
        duration = time.time() - start
        print "zCall finished. Duration:", round(duration, 2), "seconds."

def parseArgs():
    description = "Standalone script to run the complete zcall process: Generate and evaluate thresholds, and apply zcall to no-calls in the input data."
    parser = argparse.ArgumentParser(description=description)
    configDefault = os.path.join(sys.path[0], '../etc/config.ini')
    configDefault = os.path.abspath(configDefault)
    parser.add_argument('--bpm', required=True, metavar="PATH", 
                        help="BPM .csv manifest file")
    parser.add_argument('--egt', required=True, metavar="PATH", 
                        help="Path to .egt input file.")
    parser.add_argument('--config', metavar="PATH", default=configDefault,
                        help="Path to .ini config file. Default = etc/config.ini")
    parser.add_argument('--out', metavar="DIR", default=".",
                        help="Directory for output; defaults to current working directory.")
    parser.add_argument('--plink', default='zcall', metavar="STRING", 
                        help="Prefix for Plink output files")
    parser.add_argument('--zstart', metavar="INT", default=7, type=int,
                    help='Starting z score. Default = %(default)s')
    parser.add_argument('--ztotal', metavar="INT", default=1, type=int,
                        help='Total number of integer z scores to generate. Default = %(default)s')
    parser.add_argument('--gtc_start', metavar="INT", default = 0,
                        help="Starting index in GTC .json file for threshold evaluation")
    parser.add_argument('--gtc_end', metavar="INT", default = -1,
                        help="Ending index in GTC .json file for threshold evaluation")
    parser.add_argument('--samples', required=True, metavar="PATH", 
                        help="Path to .json file containing sample URIs (unique identifiers), gender codes, and .gtc data paths")
    parser.add_argument('--verbose', action='store_true', default=False,
                        help="Print status information to standard output")
    args = vars(parser.parse_args())
    inputKeys = ('bpm', 'egt', 'config')
    for key in inputKeys:
        if not os.access(args[key], os.R_OK):
            msg = "Cannot read path: \""+args[key]+"\"\n"
            sys.stderr.write(msg)
            sys.exit(1)
        else:
            args[key] = os.path.abspath(args[key])
    if not os.path.isdir(args['out']) or  not os.access(args['out'], os.W_OK):
        msg = "Output path \""+args['out']+"\" is not a writable directory!\n"
        sys.stderr.write(msg)
        sys.exit(1)
    

    return args

if __name__ == "__main__":
    main()
