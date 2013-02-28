#! /usr/bin/env python


"""Unit tests for zcall module

Test cases:
- prepareThresholds.py creates correct thresholds.txt files
- evaluateThresholds.py creates correct evaluation in .json format
- mergeEvaluation.py correctly merges results and identifies best Z score
- call.py produces Plink binary data

Required input data:
- Example BPM, EGT files
- Example GTC files -- may have privacy issues
- (TODO Test on public GTC data, with corresponding BPM/EGT)

Author:  Iain Bancarz, ib5@sanger.ac.uk
"""

import json, os, sys, unittest
from hashlib import md5
from tempfile import mkdtemp

class TestScripts(unittest.TestCase):

    """Test command-line python scripts used in WTSI genotyping pipeline"""

    def getMD5(self, inPath):
        """Get MD5 checksum for contents of given file"""
        m = md5()
        m.update(open(inPath).read())
        checksum = m.digest()
        return checksum

    def setUp(self):
        """Check for input/output directory called 'data' 

        Must have appropriate input files"""
        self.dataDir = 'data'
        if not os.access(self.dataDir, os.W_OK) \
                or not os.path.isdir(self.dataDir):
            msg = "Invalid test directory! Should run from: zCall/src\n"
            sys.stderr.write(msg)
            sys.exit(1)
        self.outDir = mkdtemp(dir=self.dataDir)
        print "Created output directory", self.outDir
        self.bpmPath = os.path.join(self.dataDir, 'HumanExome-12v1_A.bpm.csv')
        self.egtPath = os.path.join(self.dataDir, 'HumanExome-12v1.egt')
        self.gtcPaths = []
        for name in ('gtc00.json', 'gtc01.json'): 
            self.gtcPaths.append(os.path.join(self.dataDir, name))
        inPaths = [self.bpmPath, self.egtPath]
        inPaths.extend(self.gtcPaths)
        for inPath in inPaths:
            if not os.access(inPath, os.R_OK):
                sys.stderr.write("Cannot access test input '"+inPath+"'!\n")
                sys.exit(1)
        self.thresholdJsonName = 'thresholds.json'
        self.thresholdJson = os.path.join(self.dataDir, self.thresholdJsonName)

    def test_prepareThresholds(self):
        """Prepare thresholds.txt files"""
        zstart = 6
        ztotal = 3
        outPaths = []
        for i in range(zstart, zstart+ztotal):
            name = 'thresholds_HumanExome-12v1_z0'+str(i)+'.txt'
            outPaths.append(os.path.join(self.outDir, name))
        args = ['zcall/prepareThresholds.py',
                '--egt', self.egtPath,
                '--out', self.outDir,
                '--config etc/config.ini',
                '--zstart', str(zstart),
                '--ztotal', str(ztotal),
                '--index_name', self.thresholdJsonName]
        self.assertEqual(os.system(' '.join(args)), 0) # run script
        jsonOut = os.path.join(self.outDir, self.thresholdJsonName)
        self.assertTrue(os.access(jsonOut, os.R_OK))
        oldIndex = json.loads(open(self.thresholdJson).read())
        outPath = os.path.join(self.outDir, self.thresholdJsonName)
        newIndex = json.loads(open(outPath).read())
        for z in newIndex.keys():
            self.assertEqual(self.getMD5(oldIndex[z]), self.getMD5(newIndex[z]))

    def test_evaluateThresholds(self):
        """Evaluate thresholds for collections of sample GTC files"""
        argsBase = ('zcall/evaluateThresholds.py',
                    '--egt', self.egtPath,
                    '--bpm', self.bpmPath,
                    '--thresholds', self.thresholdJson)
        for i in range(len(self.gtcPaths)):
            args = list(argsBase)
            name = 'metrics0'+str(i)+'.json'
            outPath = os.path.join(self.outDir, name)
            args.extend(['--gtc', self.gtcPaths[i], '--out', outPath ])
            self.assertEqual(os.system(' '.join(args)), 0) # run script
            metricsNew = json.loads(open(outPath).read())
            oldPath = os.path.join(self.dataDir, name)
            metricsOld = json.loads(open(oldPath).read())
            self.assertEqual(metricsOld, metricsNew)

    def test_mergeEvaluation(self):
        """Merge evaluation results and find best Z score"""
        outPath = os.path.join(self.outDir, 'zEvaluation.json')
        args = ['zcall/mergeEvaluation.py',
                '--metrics', 'data/metrics.txt',
                '--thresholds', self.thresholdJson,
                '--out', outPath ]
        print ' '.join(args)
        self.assertEqual(os.system(' '.join(args)), 0) # run script
        resultsNew = json.loads(open(outPath).read())
        oldPath = os.path.join(self.dataDir, 'zEvaluation.json')
        resultsOld = json.loads(open(oldPath).read())
        self.assertEqual(resultsOld, resultsNew)

    def test_call(self):
        """Re-call GTC files using zCall"""
        #TODO read output data into PLINK to verify it is well-formatted
        outPath = os.path.join(self.outDir, 'test.bed')
        tPath = os.path.join(self.dataDir, 'thresholds_HumanExome-12v1_z07.txt')
        args = ['zcall/call.py',
                '--thresholds', tPath,
                '--bpm', self.bpmPath,
                '--egt', self.egtPath,
                '--gtc', os.path.join(self.dataDir, 'gtc.json'),
                '--out', outPath,
            ]
        self.assertEqual(os.system(' '.join(args)), 0) # run script
        checksum = self.getMD5(outPath)
        expected = '\x95A\xec\xbd\x87K\x80\xa3\x01\xb2"\x95\xea\x05\xe8\x88'
        self.assertEqual(checksum, expected)
        

if __name__ == "__main__":
    unittest.main(verbosity=2)
