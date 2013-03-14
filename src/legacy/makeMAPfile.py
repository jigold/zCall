#! /usr/bin/python

import sys
sys.path.append(sys.path[0]+'/../zcall')
from BPM import *
from optparse import OptionParser

### Parse Inputs from Command Line
parser = OptionParser()
parser.add_option("-B","--bpm",type="string",dest="bpm",action="store",help="bpm.csv file with normID and chr, pos")
(options, args) = parser.parse_args()

if options.bpm == None:
    print "specify BPM file path with -B"
    sys.exit()

### Initialize BPM file class
bpm = BPM(options.bpm)    
numSNPs = len(bpm.names)

### Make MAP file
for i in range(numSNPs):
    snp = bpm.names[i]
    chr = bpm.chr[i]
    pos = bpm.pos[i]

    out = [chr, snp, "0", pos]
    print "\t".join(out)
