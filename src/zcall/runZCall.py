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
        zcalls = 0
        gains = 0
        for i in range(self.snpTotal):
            origCall = gtc.genotypes[i]
            if null==False \
                    and origCall == 0 \
                    and self.thresholds.getX(i) != "NA" \
                    and self.thresholds.getY(i) != "NA":
                zcall = self.call(gtc, i)
                calls.append(zcall)
                zcalls += 1
                if zcall!=0: gains += 1
            else:
                calls.append(origCall)
        return (calls, zcalls, gains)

    def run(self, samplesPath, outDir, prefix, logPath=None, verbose=False, 
            null=False, reorder=True):
        """Apply zCall to GTC files and write Plink .bed output"""
        gtcPaths = self.readSampleJson(samplesPath)
        calls = []
        ph = PlinkHandler(self.bpm)
        zcallTotal = 0
        gainsTotal = 0
        for gtcPath in gtcPaths:
            if verbose: print "Calling GTC file", gtcPath
            gtc = GTC(gtcPath, self.bpm.normID)
            (callsRaw, zcalls, gains) = self.makeCalls(gtc, null)
            zcallTotal += zcalls
            gainsTotal += gains
            calls.extend(ph.callsToBinary(callsRaw, reorder))
        outStem = os.path.join(outDir, prefix)
        ph.writeBed(calls, outStem+'.bed', verbose)
        ph.writeBim(outStem+'.bim')
        ph.writeFam(samplesPath, outStem+'.fam')
        logData = { 'total_samples': len(gtcPaths),
                    'sample_json': samplesPath,
                    'plink_output': outStem,
                    'total_snps': self.snpTotal,
                    'total_calls': self.snpTotal*len(gtcPaths),
                    'zcall_attempts': zcallTotal,
                    'zcall_gains': gainsTotal, }
        if logPath!=None:
            log = open(logPath, 'w')
            log.write(json.dumps(logData, sort_keys=True,
                      indent=4, separators=(',', ': '))+"\n")
            log.close()


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
    parser.add_argument('--out', required=True, metavar="DIR", 
                        help="Directory for output data")
    parser.add_argument('--plink', default='zcall', metavar="STRING", 
                        help="Prefix for Plink output files")
    parser.add_argument('--log', metavar="PATH", default=None,
                        help="Path for .json log output. Defaults to zcall_log.json in same directory as Plink output.")
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
    dirName = os.path.abspath(args['out'])
    if not os.access(dirName, os.R_OK) or not os.path.isdir(dirName):
        raise OSError("Invalid output directory \""+args['out']+"\"")
    else:
        args['out'] = dirName
    if args['log'] == None:
        args['log'] = os.path.join(args['out'], 'zcall_log.json')
    else:
        args['log'] = os.path.abspath(args['log'])
    if args['null']:
        print "WARNING: Null option in effect, input calls will not be changed"
    caller = SampleCaller(args['bpm'], args['egt'], args['thresholds'])
    caller.run(args['samples'], args['out'], args['plink'], args['log'], 
               args['verbose'], args['null'])

if __name__ == "__main__":
    main()
