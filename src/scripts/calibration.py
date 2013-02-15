#! /usr/bin/env python

"""Classes to find thresholds and evaluate z scores, as parameters for zCall.

Threshold finding combines: findMeanSD.py, findBetas.r, findThresholds.py
Evaluation looks for optimal concordance and gain metrics on given dataset.
Import classes to front-end scripts for calibration and calling
Author: Iain Bancarz, ib5@sanger.ac.uk, January 2013
"""

import json, os, re, sys, tempfile, time
from ConfigParser import ConfigParser
from GTC import *
from BPM import *
from EGT import *
from zCallBase import zCallBase

class ThresholdFinder:
    """Class to write threshold.txt files for given EGT input and z score(s).

    Threshold finding procedure:
    1. Run findMeanSD.py on given EGT file.  Outputs mean_sd.txt
    2. Run findBetas.r on output from (1). Outputs betas.txt
    3. Run findThresholds.py on EGT file and betas.txt, with given Z score(s).
    Outputs from (1) and (2) are written to a temporary directory, deleted on exit.

    Recommended default Z score = 7.  Suggested range of alternatives = 3 to 15.

    TODO
    Modify findMeanSD.py and findThresholds.py so they can be imported, instead of being run in a subshell
    """

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

    def run(self, egtPath, zScore=7, outDir='/tmp', verbose=True, force=False):
        outPath = os.path.join(outDir, self.thresholdFileName(egtPath, zScore))
        if os.path.exists(outPath) and force==False:
            if verbose: print outPath+" already exists; omitting calibration."
            return outPath
        scriptDir = os.path.abspath(sys.path[0])
        tempDir = tempfile.mkdtemp(prefix='zcall_')
        if verbose:
            msg = "Calibrating zCall: zscore = "+str(zScore)+"\n"+\
                "Writing temporary files to "+tempDir
            print msg
        meanSd = tempDir+'/mean_sd.txt'
        betas = tempDir+'/betas.txt'
        cmdList = [scriptDir+'/findMeanSD.py -E '+egtPath+' > '+meanSd,
                   'bash -c "'+self.rScript+' '+scriptDir+'/findBetas.r '+\
                       meanSd+' '+betas+' 1 " &> '+tempDir+'/rscript.log',
                   scriptDir+'/findThresholds.py -B '+betas+' -E '+egtPath+\
                       ' -Z '+str(zScore)+' > '+outPath,
                   ] # findBetas.r command uses bash to redirect stderr
        commandsOK = True
        for cmd in cmdList:
            if verbose: print cmd
            status = os.system(cmd)
            if status!=0: 
                if verbose: print "WARNING: Non-zero exit status."
                commandsOK = False
        if commandsOK:
            if verbose: print "Cleaning up temporary directory."
            os.system('rm -Rf '+tempDir)
        elif verbose: print "Possible error, retaining temporary directory."
        if verbose: print "Finished calibration."
        return outPath


class ConcordanceGainFinder(zCallBase):
    """Class to evaluate multiple GTC files by concordance and gain metrics"""

    # TODO instead of plain text list, use .json from sample_intensities.pl

    def concordanceRate(self, counts):
        """Find concordance rate between original and new call counts

        Ignores SNPs where original was a 'no call'"""
        [match, total] = [0,0]
        for i in range(1,4):
            for j in range(0,4):
                count = counts[(i,j)]
                total += count
                if i==j: match += count
        concord = float(match)/float(total)
        return concord

    def countCallTypes(self, gtc):
        """based on method in sampleConcordance.py

        call codes: 0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        also return rate of inclusion for SNPs"""
        included = 0
        counts = {}
        for i in range(4):
            for j in range(4): counts[(i,j)] = 0
        for i in range(gtc.numSNPs):
            nAA = self.egt.nAA[i]
            nBB = self.egt.nBB[i]
            nAB = self.egt.nAB[i]
            if not self.includeSNP(i, nAA, nBB, nAB): continue
            origCall = self.normalizeCall(gtc.genotypes[i], nAA, nBB)
            newCall = self.normalizeCall(self.call(gtc, i), nAA, nBB)
            counts[(origCall, newCall)] += 1
            included += 1
        return (included, gtc.numSNPs, counts)

    def evaluate(self, inPath, outPath, verbose=True):
        """alias for writeResults function"""
        self.writeResults(inPath, outPath, verbose)

    def findMultipleConcordances(self, inPath, verbose=True):
        """Method to find concordance for multiple GTC files

        inPath = .json file created by sample_intensities.pl
        Returns included SNP count, total SNPs, and "results" list
        List contains GTC paths, concordance, and gain
        """
        gtcPaths = []
        for sample in json.loads(open(inPath).read()):
            gtcPaths.append(sample["result"])
        results = []
        gtcTotal = len(gtcPaths)
        # snp inclusion is the same for all GTC files (depends only on EGT)
        for i in range(gtcTotal):
            if verbose: print "Evaluating GTC path %s of %s" % (i+1, gtcTotal)
            gtc = GTC(gtcPaths[i], self.bpm.normID)
            (includedSNPs, totalSNPs, counts) = self.countCallTypes(gtc)
            concord = self.concordanceRate(counts)
            gain = self.gainRate(counts)
            results.append([gtcPaths[i], concord, gain])
        return (includedSNPs, totalSNPs, results)

    def gainRate(self, counts):
        """Find rate of call gain

        Defined as no calls in original which are called by zcall"""
        [gain, total] = [0,0]
        for i in range(4):
            count = counts[(0,i)] # no call in original GTC
            total += count
            if i!=0: gain += count
        gainRate = float(gain)/float(total)
        return gainRate

    def includeSNP(self, i, nAA, nBB, nAB):
        """Should ith SNP be included in concordance calculation?

        Require autosomal SNP with MAF>=5%
        Want at least 10 points in each homozygote cluster
        Also exclude SNPs without defined zcall thresholds"""
        include = True
        chrom = self.bpm.chr[i]
        maf = self.findMAF(nAA, nBB, nAB)
        if maf < 0.05 or chrom == "X" or chrom == "Y" or nAA < 10 or nBB < 10 \
                or self.thresholdsX[i]=="NA" or self.thresholdsY[i]=="NA":
            include = False
        return include        

    def writeResults(self, gtcListPath, outPath, verbose=True, digits=3):
        """Find concordances/gains and write results to file"""
        (includedSNPs, totalSNPs, results) = \
            self.findMultipleConcordances(gtcListPath, verbose)
        includeRate = float(includedSNPs)/totalSNPs
        headers = [
            '# evaluateConcordance.py results',
            '# BPM '+self.bpmPath,
            '# EGT '+self.egtPath,
            '# THRESHOLDS '+self.threshPath,
            '# INCLUDED_SNP '+str(includedSNPs),
            '# TOTAL_SNP '+str(totalSNPs),
            '# INCLUDE_RATE_SNP '+str(round(includeRate, digits)),
            '# [Input] [Concordance on original calls] [Gain]'
            ]
        out = open(outPath, 'w')
        for header in headers: 
            out.write(header+"\n")
        for result in results:
            [inPath, concord, gain] = result
            concord = round(concord, digits)
            gain = round(gain, digits)
            out.write("%s\t%s\n" % (inPath, concord, gain))
        out.close()
        if verbose: print "Finished.\n"

class ZScoreEvaluator:
    """Find thresholds; evaluate for multiple z scores and GTC files."""

    def __init__(self, egt, bpm, configPath):
        """Constructor arguments:  EGT path, BPM path, .ini path"""
        self.egt = os.path.abspath(egt)
        self.bpm = os.path.abspath(bpm)
        self.tf = ThresholdFinder(configPath)

    def findAndEvaluate(self, gtcJson, zStart, zTotal, outDir, outName, 
                        verbose=True, force=False):
        """Main method to find and evaluate thresholds."""
        z = zStart
        allResults = []
        for i in range(zTotal):
            threshPath = self.tf.run(self.egt, z, outDir, verbose, force)
            cgf = ConcordanceGainFinder(threshPath, self.bpm, self.egt)
            (includedSNPs, totalSNPs, results) = \
                cgf.findMultipleConcordances(gtcJson, verbose)
            for result in results: result.append(z)
            allResults.append(results)
            z += 1
        # includedSNPs, totalSNPs do not depend on GTC or z score
        (bestZ, bestZType) = self.findBestZ(allResults, verbose)
        outPath = os.path.join(outDir, outName)
        self.writeResults(outPath, includedSNPs, totalSNPs, 
                          bestZ, bestZType, allResults)

    def findBestZ(self, allResults, verbose=True):
        """Find 'best' zscore from multiple GTC files and thresholds.

        The 'best' is defined as the smallest z s.t. mean concordance > mean gain (type 0); or if none exists, return z with minimum of mean gain - mean concordance (type 1)."""
        concords = {}
        gains = {}
        for results in allResults: # for each zscore
            for result in results: # for each gtc file
                [inPath, concordance, gain, z] = result
                try: concords[z].append(concordance)
                except KeyError: concords[z] = [concordance,]
                try: gains[z].append(gain)
                except KeyError: gains[z] = [gain,]
        concordanceGreaterThanGain = []
        gainMinusConcord = {}
        for z in concords.keys():
            cMean = sum(concords[z])/len(concords[z])
            gMean = sum(gains[z])/len(gains[z])
            if cMean > gMean: concordanceGreaterThanGain.append(z)
            gainMinusConcord[z] = gMean - cMean
        bestType = None
        if len(concordanceGreaterThanGain)>0: 
            bestType = 0
            best = min(concordanceGreaterThanGain)
        else: 
            bestType = 1
            leastDiff = min(gainMinusConcord.values())
            for z in gainMinusConcord.keys():
                if gainMinusConcord[z]==leastDiff:
                    best = z
                    break
        if verbose:
            print "BEST_Z", best
            print "BEST_Z_TYPE", bestType
        return (best, bestType)

    def writeResults(self, outPath, includedSNPs, totalSNPs, 
                     bestZ, bestZType, allResults, digits=3):
        """Write results to file.  Header includes summary stats."""
        includeRate = round(float(includedSNPs)/totalSNPs, digits)
        out = open(outPath, 'w')
        out.write("# EGT "+self.egt+"\n")
        out.write("# BPM "+self.bpm+"\n")
        out.write("# SNP_INCLUDE_RATE "+str(includeRate)+"\n")
        out.write("# BEST_Z "+str(bestZ)+"\n")
        out.write("# BEST_Z_TYPE "+str(bestZType)+"\n")
        headers = ['# Input', 'Concordance', 'Gain', 'Zscore']
        out.write("\t".join(headers)+"\n")
        for results in allResults:
            for result in results:
                [inPath, concord, gain, z] = result
                concord = round(concord, digits)
                gain = round(gain, digits)
                words = []
                for term in [inPath, concord, gain, z]: words.append(str(term))
                out.write("\t".join(words)+"\n")
        out.close()
