#! /usr/bin/env python

"""Convenience script to generate HTML documentation

Use instead of command-line pydoc, to ensure use of Python version 2.7
"""

import os, pydoc, re, sys

scriptDir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(scriptDir+"/..")) # get parent dir for import
os.chdir(scriptDir) # write html to createDocs.py directory
# Hard-coded list of modules is a workaround for limitations of module importing
import zcall, zcall.BPM, zcall.EGT, zcall.GTC, \
    zcall.calibration, zcall.call, zcall.evaluateThresholds, \
    zcall.mergeEvaluation, zcall.prepareThresholds, zcall.utilities
for module in (zcall, zcall.BPM, zcall.EGT, zcall.GTC, 
               zcall.calibration, zcall.call, zcall.evaluateThresholds, 
               zcall.mergeEvaluation, zcall.prepareThresholds, zcall.utilities):
    pydoc.writedoc(module)

# NB findThresholds and findMeanSD files are only ever run as scripts, not imported.  Omitting the .py extension from these files prevents pydoc from creating broken links in the main zcall page.
