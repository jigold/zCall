#! /usr/bin/env python


"""Unit tests for zcall module

Test cases:
- prepareThresholds.py creates correct thresholds.txt files (use checksum?)
- evaluateThresholds.py creates correct evaluation in .json format
- mergeEvaluation.py correctly merges results and identifies best Z score
- call.py produces valid Plink binary data
- sub-tests for BPM, EGT, GTC object creation

Required input data:
- Example BPM, EGT files
- Example GTC files -- may have privacy issues
- (TODO Test on public GTC data, with corresponding BPM/EGT)

Author:  Iain Bancarz, ib5@sanger.ac.uk
"""

import json, os, sys, unittest
from hashlib import md5

class TestScripts(unittest.TestCase):

    """Test command-line python scripts used in WTSI genotyping pipeline"""

    def setUp(self):
        """Check for input/output directory called 'data' 

        Must have appropriate input files"""
        self.testDir = 'data'
        self.refDir = self.testDir+'/ref'
        if not os.access(self.testDir, os.W_OK) \
                or not os.path.isdir(self.testDir):
            msg = "Invalid test directory! Should run from: zCall/src\n"
            sys.stderr.write(msg)
            sys.exit(1)
        self.bpmPath = os.path.join(self.testDir, 'HumanExome-12v1_A.bpm.csv')
        self.egtPath = os.path.join(self.testDir, 'HumanExome-12v1.egt')
        self.gtcPaths = []
        for name in ('gtc00.json', 'gtc01.json'): 
            self.gtcPaths.append(os.path.join(self.testDir, name))
        inPaths = [self.bpmPath, self.egtPath]
        inPaths.extend(self.gtcPaths)
        for inPath in inPaths:
            if not os.access(inPath, os.R_OK):
                sys.stderr.write("Cannot access test input '"+inPath+"'!\n")
                sys.exit(1)
        self.thresholdJsonName = 'thresholds.json'
        self.thresholdJson = os.path.join(self.testDir, self.thresholdJsonName)
        self.thresholdJsonRef = os.path.join(self.refDir,self.thresholdJsonName)

    def test_prepareThresholds(self):
        """Prepare thresholds.txt files"""
        zstart = 6
        ztotal = 3
        if os.path.exists(self.thresholdJson): os.remove(self.thresholdJson)
        outPaths = []
        for i in range(zstart, zstart+ztotal):
            name = 'thresholds_HumanExome-12v1_z0'+str(i)+'.txt'
            outPath = os.path.join(self.testDir, name)
            outPaths.append(outPath)
            if os.path.exists(outPath): os.remove(outPath)
        args = ['zcall/prepareThresholds.py',
                '--egt', self.egtPath,
                '--out', self.testDir,
                '--config etc/config.ini',
                '--zstart', str(zstart),
                '--ztotal', str(ztotal),
                '--index_name', self.thresholdJsonName]
        self.assertEqual(os.system(' '.join(args)), 0) # run script
        self.assertTrue(os.access(self.thresholdJson, os.R_OK))
        index = {}
        for i in range(ztotal):
            index[unicode(i+zstart)] = unicode(os.path.abspath(outPaths[i]))
        self.assertEqual(json.loads(open(self.thresholdJson).read()), index)
        # md5 checksums for thresholds.txt files
        # can't use md5 on json output, as order of key/value pairs may vary
        expected = ['\xaaD\x05\xf0r\xa12\x9a\xc9?\xd9\xfa\xaf\xd8\xb6z', 
                    "oL^;\x87r'\xa2\xd9%'8/5g\xa7", 
                    'U\x19\x98:\x92\xb6\x04@A\xe0\x84bE\xd5\r\xf4']
        for i in range(len(outPaths)):
            self.assertTrue(os.access(outPaths[i], os.R_OK))
            m = md5()
            m.update(open(outPaths[i]).read())
            checksum = m.digest()
            self.assertEqual(checksum, expected[i])

    def test_evaluateThresholds(self):
        """Evaluate thresholds for collections of sample GTC files"""
        argsBase = ('zcall/evaluateThresholds.py',
                    '--egt', self.egtPath,
                    '--bpm', self.bpmPath,
                    '--thresholds', self.thresholdJsonRef)
        for i in range(len(self.gtcPaths)):
            args = list(argsBase)
            outPath = os.path.join(self.testDir, 'metrics0'+str(i)+'.json')
            args.extend(['--gtc', self.gtcPaths[i], '--out', outPath ])
            self.assertEqual(os.system(' '.join(args)), 0) # run script
            metricsNew = json.loads(open(outPath).read())
            oldPath = os.path.join(self.testDir,'ref','metrics0'+str(i)+'.json')
            metricsOld = json.loads(open(oldPath).read())
            self.assertEqual(metricsOld, metricsNew)

    def test_mergeEvaluation(self):
        """Merge evaluation results and find best Z score"""
        outPath = os.path.join(self.testDir, 'zEvaluation.json')
        args = ['zcall/mergeEvaluation.py',
                '--metrics', 'data/metrics.txt',
                '--thresholds', self.thresholdJsonRef,
                '--out', outPath ]
        self.assertEqual(os.system(' '.join(args)), 0) # run script
        resultsNew = json.loads(open(outPath).read())
        oldPath = os.path.join(self.testDir,'ref','zEvaluation.json')
        resultsOld = json.loads(open(oldPath).read())
        self.assertEqual(resultsOld, resultsNew)

    def test_call(self):
        """Re-call GTC files using zCall"""

        #TODO read output data into PLINK to verify it is well-formatted
        outPath = os.path.join(self.testDir, 'test.bed')
        tPath = os.path.join(self.testDir, 'ref', 
                             'thresholds_HumanExome-12v1_z07.txt')
        args = ['zcall/call.py',
                '--thresholds', tPath,
                '--bpm', self.bpmPath,
                '--egt', self.egtPath,
                '--gtc', os.path.join(self.testDir, 'gtc.json'),
                '--out', outPath,
            ]
        self.assertEqual(os.system(' '.join(args)), 0) # run script
        m = md5()
        m.update(open(outPath).read())
        checksum = m.digest()
        expected = '\x95A\xec\xbd\x87K\x80\xa3\x01\xb2"\x95\xea\x05\xe8\x88'
        self.assertEqual(checksum, expected)
        

if __name__ == "__main__":
    unittest.main(verbosity=2)
