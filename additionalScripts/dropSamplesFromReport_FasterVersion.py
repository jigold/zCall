#! /usr/bin/python

# Thank you to Josh Randall (Sanger) for providing this script!

import sys

reportFileName = sys.argv[1]
dropSamplesFileName = sys.argv[2]

# read sample exclusion file and populate samples set 
dropSamples = set() 
for line in open(dropSamplesFileName, 'r'):
	line = line.strip()
	dropSamples.add(line)
 
# open file => returns an iterator
reportFile = open(reportFileName, 'r')

okColumns = []
# read first line and fill okColumns with column indices to be kept (the columns not present in the dropSamples)
headerline = reportFile.next() 
headers = headerline.strip().split("\t") 
headerout = headers[:3]
for i in range(3, len(headers)):
	id = headers[i].split(".")[0]
	if id not in dropSamples:
		okColumns.append(i)
		headerout.append(headers[i])

# print new header
print "\t".join(headerout)
 
# read the rest of the report file lines, printing only the important columns (okColumns) 
for line in reportFile:
	fields = line.strip().split("\t")
	out = fields[:3]
	for i in okColumns:
		out.append(fields[i])
	print "\t".join(out)

