#! /usr/bin/env python

"""Convenience script to generate HTML documentation

Use instead of command-line pydoc. Allows user to specify Python version 2.7, even if installation of pydoc is compiled with an earlier Python version.
"""

import os, pydoc, re, sys
from importlib import import_module
from modulefinder import ModuleFinder
try: 
    import argparse    
except ImportError: 
    sys.stderr.write("ERROR: Requires Python 2.7 to run; exiting.\n")
    sys.exit(1)

description = "Convenience script to generate HTML documentation for zCall"
parser = argparse.ArgumentParser(description=description)
parser.add_argument('--recursive', action='store_true', default=False,
                    help="Recursively import documentation for dependencies. If not recursive, zcall documents will contain broken links to standard modules.")
parser.add_argument('--verbose', action='store_true', default=False,
                    help="Write pydoc information to stdout.")
args = vars(parser.parse_args())
recursive = args['recursive']
verbose = args['verbose']
if not verbose: # suppress stdout chatter from pydoc.writedoc
    sys.stdout = open('/dev/null', 'w')

localDir = os.path.dirname(os.path.realpath(__file__))
zcallDir = os.path.abspath(localDir+"/../zcall")
sys.path.append(os.path.abspath(localDir+"/..")) # allows import from zcall dir
os.chdir(localDir) # write html to createDocs.py directory

import zcall
pydoc.writedoc(zcall)
modules = set()
mf = ModuleFinder()
for script in os.listdir(zcallDir):
    if re.search("\.py$", script) and script!="__init__.py":
        words = re.split("\.", script)
        words.pop()
        name = (".".join(words))
        modules.add("zcall."+name)
        if recursive:
            mf.run_script(os.path.join(zcallDir, script))
            for name, mod in mf.modules.iteritems(): modules.add(name)
for module in modules:
    pydoc.writedoc(import_module(module))

# NB findThresholds and findMeanSD files are only ever run as scripts, not imported.  Omitting the .py extension from these files prevents pydoc from creating broken links in the main zcall page.
