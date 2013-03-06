#! /usr/bin/env python

"""Classes to find thresholds and evaluate z scores, as parameters for zCall.

Threshold finding combines: findMeanSD.py, findBetas.r, findThresholds.py
Evaluation looks for optimal concordance and gain metrics on given dataset.
Import classes to front-end scripts for calibration and calling.

Contents:
- ThresholdFinder
- MetricEvaluator
- MetricFinder
- SampleEvaluator

Author: Iain Bancarz, ib5@sanger.ac.uk, January 2013
"""

import json, os, re, sys, tempfile
from ConfigParser import ConfigParser
from GTC import *
from BPM import *
from EGT import *
from utilities import CallingBase, ThresholdContainer, SharedBase

class ThresholdFinder:
    """Class to write threshold.txt files for given EGT input and z score(s).

    Threshold finding procedure:
    1. Run findMeanSD on given EGT file.  Outputs mean_sd.txt
    2. Run findBetas.r on output from (1). Outputs betas.txt
    3. Run findThresholds on EGT file and betas.txt, with given Z score(s).
    Outputs from (1) and (2) are written to a temporary directory, deleted on successful exit.

    Recommended default Z score = 7.  Suggested range of alternatives = 3 to 15.
"""

#    TODO  Modify findMeanSD and findThresholds so they can be imported, instead of being run in a subshell
 
    def __init__(self, configPath=None):
        if configPath==None:
            configPath = os.path.join(sys.path[0], '../etc/config.ini')
            configPath = os.path.abspath(configPath)
        config = ConfigParser()
        config.readfp(open(configPath))
        self.rScript = config.get('zcall', 'rscript')
        self.digits = int(config.get('zcall', 'digits'))

    def thresholdFileName(self, egtPath, zScore):
        egtName = re.split('/', egtPath).pop()
        items = re.split('\.', egtName)
        items.pop() # remove .egt suffix
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
        cmdList = [scriptDir+'/findMeanSD -E '+egtPath+' > '+meanSd,
                   'bash -c "'+self.rScript+' '+scriptDir+'/findBetas.r '+\
                       meanSd+' '+betas+' 1 " &> '+tempDir+'/rscript.log',
                   scriptDir+'/findThresholds -B '+betas+' -E '+egtPath+\
                       ' -Z '+str(zScore)+' -D '+str(self.digits)+\
                       ' > '+outPath,
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


class MetricEvaluator(SharedBase):
    """Class to assess concordance/gain metrics and choose best z score"""
    
    def __init__(self):
        pass

    def findBestZ(self, concords, gains):
        """Find best z score from mean concordance/gain values

        The 'best' is defined as the smallest z s.t. mean concordance > mean gain; or if none exists, return z with minimum of gain - concordance"""
        concordanceGreaterThanGain = []
        gainMinusConcord = {}
        for z in concords.keys():
            concord = concords[z]
            gain = gains[z]
            if concord > gain: concordanceGreaterThanGain.append(z)
            gainMinusConcord[z] = gain - concord
        best = None
        bestType = 0
        if len(concordanceGreaterThanGain)>0:
            best = min(concordanceGreaterThanGain)
        else:
            leastDiff = min(gainMinusConcord.values())
            for z in gainMinusConcord.keys():
                if gainMinusConcord[z] == leastDiff:
                    best = z
                    bestType = 1
                    break
        return best

    def findMeans(self, inPaths, outPath=None):
        """Read JSON result paths, find mean concordance/gain by z score

        Return concatenation of results and concordance/gain metrics
        """
        rows = []
        for inPath in inPaths:
            rows.extend(json.loads(open(inPath).read()))
        zCounts = {}
        concords = {}
        gains = {}
        for row in rows:
            [gtc, z, concord, gain, counts] = row
            try: 
                zCounts[z] += 1
                concords[z] += concord
                gains[z] += gain
            except KeyError: 
                zCounts[z] = 1
                concords[z] = concord
                gains[z] = gain
        for z in zCounts.keys():
            concords[z] = concords[z] / float(zCounts[z])
            gains[z] = gains[z] / float(zCounts[z])
        if outPath!=None:
            out = open(outPath, 'w')
            out.write(json.dumps(rows))
            out.close()
        return (rows, concords, gains)

    def writeBest(self, resultsPath, thresholdPath, outPath):
        """Find best z score & thresholds.txt, write to file for later use

        Arguments:
        - Text file listing SampleEvaluator output paths, one per line
        - JSON file with hash of thresholds.txt paths by z score
        - Output directory
        """
        inPaths = []
        lines = open(resultsPath).readlines()
        for line in lines:
            line = line.strip()
            if line!='': inPaths.append(line)
        (metrics, concords, gains) = self.findMeans(inPaths)
        best = self.findBestZ(concords, gains)
        thresholdPaths = json.loads(open(thresholdPath).read())
        results = { self.Z_KEY:best, self.T_KEY:thresholdPaths[best],
                    self.M_KEY:metrics}
        out = open(outPath, 'w')
        out.write(json.dumps(results))
        out.close()
        

class MetricFinder(CallingBase):
    """Class to evaluate GTC objects by concordance and gain metrics

    Initialize with egt path, bpm path
    Inherits common "calling" functions from CallingBase
"""

    def concordanceRate(self, counts):
        """Find concordance rate between original and new calls

        Ignores SNPs where original was a 'no call'"""
        [match, total] = [0,0]
        for i in range(1,4):
            for j in range(0,4):
                count = counts[(i,j)]
                total += count
                if i==j: match += count
        concord = float(match)/float(total)
        return concord

    def countCallTypes(self, thresholds, gtc):
        """Count call types (0, AA, AB, BB) for given GTC and thresholds.

        call codes: 0 - "No Call", 1 - AA, 2 - AB, 3 - BB
        based on method in sampleConcordance.py
        """
        self.setThresholds(thresholds)
        counts = {}
        for i in range(4):
            for j in range(4): counts[(i,j)] = 0
        includedSNPs = self.findIncludedSNPs(thresholds)
        for i in includedSNPs:
            nAA = self.egt.nAA[i]
            nBB = self.egt.nBB[i]
            nAB = self.egt.nAB[i]
            origCall = self.normalizeCall(gtc.genotypes[i], nAA, nBB)
            newCall = self.normalizeCall(self.call(gtc, i), nAA, nBB)
            counts[(origCall, newCall)] += 1
        return counts

    def getMetrics(self, thresholds, gtc):
        """Find call types, concordance and gain for given threshold and GTC

        Arguments are ThresholdContainer and GTC objects"""
        counts = self.countCallTypes(thresholds, gtc)
        concord = self.concordanceRate(counts)
        gain = self.gainRate(counts)
        return (counts, concord, gain)

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

    def includeSNP(self, i, nAA, nBB, nAB, thresholds):
        """Should ith SNP be included in concordance calculation?

        Require autosomal SNP with MAF>=5%
        Want at least 10 points in each homozygote cluster
        Also exclude SNPs without defined zcall thresholds"""
        include = True
        chrom = self.bpm.chr[i]
        maf = self.findMAF(nAA, nBB, nAB)
        if maf < 0.05 or chrom == "X" or chrom == "Y" or nAA < 10 or nBB < 10 \
                or thresholds.getX(i)=="NA" or thresholds.getY(i)=="NA":
            include = False
        return include        

    def findIncludedSNPs(self, thresholds):
        """ Find set of included SNP indices """
        included = []
        for i in range(self.snpTotal):
            nAA = self.egt.nAA[i]
            nBB = self.egt.nBB[i]
            nAB = self.egt.nAB[i]
            if self.includeSNP(i, nAA, nBB, nAB, thresholds): 
                included.append(i)
        return included

class SampleEvaluator(SharedBase):
    """Evaluate z scores and thresholds for one or more GTC files."""

    def __init__(self, bpmPath, egtPath):
        self.bpmPath = bpmPath
        self.egtPath = egtPath
        self.bpm = BPM(bpmPath)
        self.metricFinder = MetricFinder(bpmPath, egtPath)

    def convertCountKeys(self, counts):
        """Convert keys in counts dictionary to string; required for JSON output

        Counts are indexed by (original call, final call)
        """
        output = {}
        for key in counts.keys():
            (i,j) = key
            output[str(i)+':'+str(j)] = counts[key]
        return output

    def evaluate(self, thresholds, gtc, verbose=False):
        """Evaluate z thresholds for given thresholds & sample GTC

        Inputs:
        - dictionary of threshold paths indexed by z score
        - GTC object
        """
        gtcName = os.path.split(gtc.getInputPath())[1]
        if verbose: print "Evaluating z scores for sample", gtcName
        zList = thresholds.keys()
        zList.sort()
        results = {}
        for z in zList:
            if verbose: print "Finding metrics for z score", z
            results[z] = self.metricFinder.getMetrics(thresholds[z], gtc)
        output = []
        for z in zList:
            (counts, concord, gain) = results[z]           
            converted = self.convertCountKeys(counts)
            output.append([gtcName, z, concord, gain, converted])
        return output

    def run(self, thresholdPath, sampleJson, start, end, outPath,verbose=False):
        """Evaluate thresholds for given list of sample GTC paths

        Inputs:
        - Path to .json file with hash of paths to threshold.txt files
        - Path to .json file with paths of sample GTC files
        - Start index in GTC .json file (use to split GTC for batch processing)
        - End index in GTC .json file
        - Output path

        Output:
        - JSON file with GTC filename, z score, metrics, and call type counts

        GTC input is in sample JSON format used by genotyping pipeline """
        metrics = []
        if verbose: print "Reading thresholds."
        thresholdPaths = json.loads(open(thresholdPath).read())
        thresholds = {}
        for z in thresholdPaths.keys(): 
            thresholds[z] = ThresholdContainer(thresholdPaths[z])
        if verbose: print "Evaluating samples."
        output = []
        gtcPaths = self.readSampleJson(sampleJson)
        if start>=end:
            raise ValueError("Must have GTC start index < end index")
        elif start < 0:
            raise ValueError("Must have GTC start index > 0")
        elif end > len(gtcPaths):
            raise ValueError("Must have GTC end index <= total GTC paths")
        gtcPaths = gtcPaths[start:end]
        for gtcPath in gtcPaths:
            gtc = GTC(gtcPath, self.bpm.normID)
            output.extend(self.evaluate(thresholds, gtc, verbose))
        # write results
        out = open(outPath, 'w')
        out.write(json.dumps(output))
        out.close()
