#! /usr/bin/env python

import os, struct
from GTC import *
try: 
    import argparse, json
    from calibration import ThresholdFinder
    from utilities import CallingBase, ThresholdContainer
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)


class SampleCaller(CallingBase):

    """Class to run zCall on given GTC files

    initialize with: bpmPath, egtPath, threshPath
    Write output in .bed format; use merge_bed in genotyping pipeline workflows
    """

    def callsToBinary(self, calls):
        """Translate genotype calls for one or more samples to Plink binary

        4 genotype calls are packed into one byte of output
        Returns a list of struct.pack strings corresponding to output bytes"""
        callTotal = len(calls)
        if callTotal % self.snpTotal != 0:
            raise ValueError("Number of calls is not a multiple of SNP total!")
        elif callTotal % 4 != 0:
            # if not an integer number of bytes, pad with no calls
            calls.extend([0]*(callTotal % 4)) 
        output = []
        i = 0
        while i < callTotal:
            byte = struct.pack('B', self.callsToByte(calls[i:i+4]))
            #print calls[i:i+4], byte
            output.append(byte)
            i += 4
        return output

    def callsToByte(self, calls):
        """Convert list of 4 calls to an integer in Plink binary format 

        Create byte string of the form '01001101', convert to integer
        See http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml
        """
        if len(calls) != 4:
            raise ValueError("Must have exactly 4 calls for byte conversion!")
        byte = []
        for call in calls:
            if call==1: bcall = '00' # major homozygote, 'AA'
            elif call==2: bcall = '01' # heterozygote, 'AB'
            elif call==3: bcall = '11' # minor homozygote, 'BB'
            else: bcall = '10' # missing genotype
            byte.append(bcall)
        byteString = ''.join(byte)
        byteString[::-1] # reverse order of string characters
        return int(byteString, 2)

    def makeCalls(self, gtc):
        """Apply zCall to 'no calls' from GTC input

        Return genotypes in numeric format
        0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        """
        calls = []
        for i in range(self.snpTotal):
            origCall = gtc.genotypes[i]
            if origCall == 0 \
                    and self.thresholds.getX(i) != "NA" \
                    and self.thresholds.getY(i) != "NA":
                calls.append(self.call(gtc, i))
            else:
                calls.append(origCall)
        return calls

    def run(self, sampleJsonPath, outPath, verbose=False):
        """Apply zCall to GTC files and write Plink .bed output"""
        gtcPaths = self.readSampleJson(sampleJsonPath)
        calls = []
        for gtcPath in gtcPaths:
            if verbose: print "Calling GTC file", gtcPath
            gtc = GTC(gtcPath, self.bpm.normID)
            calls.extend(self.makeCalls(gtc))
        self.writeBed(calls, outPath)
            
    def writeBed(self, calls, outPath):
        """Write output for one or more samples in Plink .bed format

        Input: List of lists of genotype call codes, and output path
        Output file:  First 2 bytes are Plink magic number
        3rd byte is flag for an individual-major file
        Subsequent bytes represent genotype calls
        """
        header = [0b01101100, 0b00011011, 0b00000000]
        output = []
        for byte in header: output.append(struct.pack('B', byte))
        output.extend(self.callsToBinary(calls))
        out = open(outPath, 'w')
        for byte in output: out.write(byte)
        out.close()


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
    parser.add_argument('--gtc', required=True, metavar="PATH", 
                        help="Path to .json file containing .gtc input paths")
    parser.add_argument('--out', required=True, metavar="PATH", 
                        help="Path for Plink .bed binary output")
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
    else:
        args['out'] = os.path.abspath(args['out'])
    caller = SampleCaller(args['bpm'], args['egt'], args['thresholds'])
    caller.run(args['gtc'], args['out'])

if __name__ == "__main__":
    main()
