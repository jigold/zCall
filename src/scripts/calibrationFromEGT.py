#! /usr/bin/env python

# Find thresholds for given .egt file and Z score(s)
# Convenience script to combine findMeanSD.py, findBetas.r, findThresholds.py

# Iain Bancarz, ib5@sanger.ac.uk
# January 2013

import os, re, sys, getopt, tempfile
from ConfigParser import ConfigParser

# optparse is deprecated; using getopt for compatibility with Python < 2.7

defaultZstart = 7
defaultZincr = 1
defaultZtotal = 1

help = """Usage: calibration.py [options]

Generates threshold files for use with the zcall genotype caller.  Inputs are 
a .egt file and one or more Z score parameters.  The .egt file is a 
proprietary Illumina binary file format, containing typical means and 
standard deviations for intensity clusters.  An .egt file is supplied by 
Illumina for its own genotyping chips, or it may be generated using the 
GenomeStudio software for custom probe sets.

Options:
--egt      Path to .egt input file.
--out      Directory for output file(s).  Filename(s) will be of the form 
           <prefix>_z<zscore>_thresholds.txt, for an input file of the form
           <prefix>.egt 
--zstart   Starting z score.  Default = """+str(defaultZstart)+"""
--ztotal   Total number of z scores to generate.  Default = """+\
    str(defaultZtotal)+"""
--zincr    Increment between z scores, if thresholds are to be generated for 
           more than one score.  Default = """+str(defaultZincr)+"""
--verbose  Print additional status information to stderr

Calibration procedure:
1. Run findMeanSD.py on given EGT file.  Outputs mean_sd.txt
2. Run findBetas.r on output from (1). Outputs betas.txt
3. Run findThresholds.py on EGT file and betas.txt, with given Z score(s).
Outputs from (1) and (2) are written to a temporary directory, deleted on exit.
"""

"""
Recommended default Z score = 7.  Suggested range of alternatives = 3 to 15.

TODO
Modify findMeanSD.py and findThresholds.py so they can be imported, instead of being run in a subshell
Get EGT file, zscore, and output directory from command line options; also print help text
Run with multiple z scores using start, increment, total?
"""


class calibration:

    def __init__(self, configPath=None):
        if configPath==None:
            configPath = os.path.join(sys.path[0], '../etc/config.ini')
            configPath = os.path.abspath(configPath)
        config = ConfigParser()
        config.readfp(open(configPath))
        self.rScript = config.get('zcall', 'rscript')

    def thresholdFileName(self, egtPath, zScore):
        egtName = re.split('/', egtPath).pop()
        items = re.split('\.', egtName)
        items.pop()
        name = '.'.join(items)
        return 'thresholds_'+name+'_z'+str(zScore).zfill(2)+'.txt'

    def run(self, egtPath, zScore=7, outDir='/tmp', verbose=True):
        scriptDir = os.path.abspath(sys.path[0])
        tempDir = tempfile.mkdtemp(prefix='zcall_')
        if verbose:
            sys.stderr.write("Writing temporary files to "+tempDir+"\n")
        meanSd = tempDir+'/mean_sd.txt'
        betas = tempDir+'/betas.txt'
        thresholds = self.thresholdFileName(egtPath, zScore)
        cmdList = [scriptDir+'/findMeanSD.py -E '+egtPath+' > '+meanSd,
                   self.rScript+' '+scriptDir+'/findBetas.r '+meanSd+' '+\
                       betas+' 1',
                   scriptDir+'/findThresholds.py -B '+betas+' -E '+egtPath+\
                       ' -Z '+str(zScore)+' > '+outDir+'/'+thresholds,
                   ]
        for cmd in cmdList:
            if verbose: sys.stderr.write(cmd+"\n")
            os.system(cmd)
        if verbose: sys.stderr.write("Cleaning up temporary directory.\n")
        os.system('rm -Rf '+tempDir)
        if verbose: sys.stderr.write("Finished.\n")


def main(egtPath, zStart, zIncr, zTotal, verbose=True):
    z = zStart
    cal = calibration()
    for i in range(zTotal):
        cal.run(egtPath, z, verbose)
        z += zIncr

if __name__ == "__main__":
    print help
    sys.exit(0)

    #main(sys.argv[1])
