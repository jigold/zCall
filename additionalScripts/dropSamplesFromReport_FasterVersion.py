#! /usr/bin/python

# Based on version from Josh Randall

import sys

reportFileName = sys.argv[1]
dropSamplesFileName = sys.argv[2]

# read sample exclusion file and populate samples set
dropSamples = set([i.strip() for i in open(dropSamplesFileName)])

# open file => returns an iterator
reportFile = open(reportFileName)

tab = "\t"
# read first line and fill okColumns with column indices to be kept (the columns not present in the dropSamples)
headers = reportFile.next().strip().split(tab) 
nfields = len(headers)
okColumns = [0,1,2] + [i for i in range(3, nfields) if headers[i].split(".")[0] not in dropSamples]

# rewind file for inclusion of header in uniform filtering
reportFile.seek(0)
# read the whole report file, printing only the important columns (okColumns) 
for line in reportFile:
	fields = line.strip().split(tab)
	print tab.join([fields[i] for i in okColumns])

