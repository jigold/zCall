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
from ConfigParser import ConfigParser
from hashlib import md5
from tempfile import mkdtemp

class TestScripts(unittest.TestCase):

    """Test command-line python scripts used in WTSI genotyping pipeline"""

    def getMD5hex(self, inPath):
        """Get MD5 checksum for contents of given file, in hex format"""
        m = md5()
        m.update(open(inPath).read())
        checksum = m.hexdigest()
        return checksum

    def readConfig(self, configPath=None):
        """Read local params from config file

        - bigdata: Directory for test files too big to upload to github"""
        if configPath==None:
            configPath = os.path.abspath('etc/config.ini')
        if not os.access(configPath, os.R_OK):
            raise IOError("Cannot read config path '"+configPath+"'")
        config = ConfigParser()
        config.readfp(open(configPath))
        bigData = config.get('test', 'bigdata')
        return bigData

    def validateThresholds(self, jsonPath):
        """Check that thresholds.txt files are correct

        Use for test_prepareThresholds, and to validate input for other tests
        """
        index = json.loads(open(jsonPath).read())
        for z in index.keys():
            self.assertEqual(self.getMD5hex(index[z]), self.expectedT[z])

    def setUp(self):
        """Check for valid input/output files and directories"""
        self.dataDir = 'data'
        self.bigData = self.readConfig()
        for d in (self.dataDir, self.bigData):
            if not os.path.exists(d) or not os.path.isdir(d):
                msg = "Invalid test directory: \""+d+"\"\n"
                sys.stderr.write(msg)
                sys.exit(1)
        self.expectedT = {"6":"1a53e8cbba6750d43d5ff607cf616beb",
                          "7":"a8d8b62be728b62fc986230da13f4ef7",
                          "8":"1f14419d0053841cfa8ab3fb994de1c1"}
        self.outDir = mkdtemp(dir=self.dataDir)
        print "Created output directory", self.outDir
        self.bpmPath = os.path.join(self.bigData, 'HumanExome-12v1_A.bpm.csv')
        self.egtPath = os.path.join(self.bigData, 'HumanExome-12v1.egt')
        self.gtcPaths = []
        for name in ('gtc00.json', 'gtc01.json'): 
            self.gtcPaths.append(os.path.join(self.dataDir, name))
        self.thresholdJsonName = 'thresholds.json'
        self.thresholdJson = os.path.join(self.bigData, self.thresholdJsonName)
        if os.path.exists(self.thresholdJson):
            self.validateThresholds(self.thresholdJson)
        else:
            msg = "WARNING: Must generate zcall thresholds with .json index: "+\
                self.thresholdJson+"\nSee output from test_prepareThresholds."
            sys.stderr.write(msg)

    def test_prepareThresholds(self):
        """Prepare thresholds.txt files

        Run as part of normal test suite
        Can also be used to generate input thresholds for other tests
        Checksums should ensure that generated thresholds are correct
        """
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
        self.validateThresholds(jsonOut)

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
        self.assertEqual(os.system(' '.join(args)), 0) # run script
        resultsNew = json.loads(open(outPath).read())
        oldPath = os.path.join(self.dataDir, 'zEvaluation.json')
        resultsOld = json.loads(open(oldPath).read())
        for key in ('BEST_Z', 'SAMPLE_METRICS'):
            self.assertEqual(resultsOld[key], resultsNew[key])
        newT = resultsNew['BEST_THRESHOLDS'] # threshold.txt path
        self.assertEqual(self.expectedT["7"], self.getMD5hex(newT))

    def test_call(self):
        """Re-call GTC files using zCall"""
        #TODO read output data into PLINK to verify it is well-formatted
        outPath = os.path.join(self.outDir, 'test.bed')
        tPath = os.path.join(self.bigData, 'thresholds_HumanExome-12v1_z07.txt')
        args = ['zcall/call.py',
                '--thresholds', tPath,
                '--bpm', self.bpmPath,
                '--egt', self.egtPath,
                '--gtc', os.path.join(self.dataDir, 'gtc.json'),
                '--out', outPath,
            ]
        self.assertEqual(os.system(' '.join(args)), 0) # run script
        checksum = self.getMD5hex(outPath)
        expected = '9541ecbd874b80a301b22295ea05e888'
        self.assertEqual(checksum, expected)
        

if __name__ == "__main__":
    unittest.main(verbosity=2)
