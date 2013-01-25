#! /usr/bin/env python

# Script to calibrate zCall
# Finds thresholds for given .egt file and Z score(s)

# Iain Bancarz, ib5@sanger.ac.uk
# January 2013

import os, re, sys, getopt, tempfile

# optparse is deprecated; using getopt for compatibility with Python < 2.7

"""
Calibration procedure:
1. Run findMeanSD.py on given EGT file, eg. myprobeset.egt
2. Run findBetas.r on output from (1). Outputs betas.txt
3. Run findThresholds.py on EGT file and betas.txt, with given Z score(s). Outputs [EGT_prefix]_z[zscore]_thresholds.txt, eg. myprobeset_z07_thresholds.txt

Recommended default Z score = 7.  Suggested range of alternatives = 3 to 15.

TODO
Modify findMeanSD.py and findThresholds.py so they can be imported, instead of being run in a subshell
Read location of R script executable from a config file (don't want to hard-code in R script shebang)
"""


class calibration:

    def __init__(self):
        self.rScript = '/software/R-2.14.1/bin/Rscript'

    def cleanup(self, tempDir):
        # remove contents of temporary directory
        # does *not* remove subdirectories, throws exception if any are present
        for name in os.listdir(tempDir):
            os.remove(os.path.join(tempDir, name))
        os.rmdir(tempDir)

    def thresholdFileName(self, egtPath, zScore):
        egtName = re.split('/', egtPath).pop()
        items = re.split('\.', egtName)
        items.pop()
        name = '.'.join(items)
        return 'thresholds_'+name+'_z'+str(zScore).zfill(2)+'.txt'

    def run(self, egtPath, zScore=7, outDir='/tmp'):

        scriptDir = os.path.abspath(sys.path[0])
        tempDir = tempfile.mkdtemp(prefix='zcall_')
        sys.stderr.write("Writing temporary files to "+tempDir+"\n")
        meanSd = tempDir+'/mean_sd.txt'
        betas = tempDir+'/betas.txt'
        thresholds = self.thresholdFileName(egtPath, zScore)
        cmd = scriptDir+'/findMeanSD.py -E '+egtPath+' > '+meanSd
        sys.stderr.write(cmd+"\n")
        os.system(cmd)
        cmd = self.rScript+' '+scriptDir+'/findBetas.r '+meanSd+' '+betas+' 1'
        sys.stderr.write(cmd+"\n")
        os.system(cmd)
        cmd = scriptDir+'/findThresholds.py -B '+betas+' -E '+egtPath+' -Z '+\
            str(zScore)+' > '+outDir+'/'+thresholds
        sys.stderr.write(cmd+"\n")
        os.system(cmd)
        sys.stderr.write("Cleaning up temporary directory.\n")
        self.cleanup(tempDir)
        sys.stderr.write("Finished.\n")
        


def main(egtPath):

    calibration().run(egtPath)

if __name__ == "__main__":
    main(sys.argv[1])
