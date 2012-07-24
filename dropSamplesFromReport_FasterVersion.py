#! /usr/bin/python

# Thank you to Josh Randall (Sanger) for providing this script!

import sys

reportFileName = sys.argv[1]
dropSamplesFileName = sys.argv[2]

# read sample exclusion file and populate samples dictionary 
dropSamples = dict() 
for line in open(dropSamplesFileName, 'r'):
	line = line.replace("\n", "")
	dropSamples[line] = 1
 
# open file and connect to iterator
reportIter = iter(open(reportFileName, 'r').readline, '')

dropColumns = dict()
# read first line and fill dropColumns with column indices to drop 
headerline = reportIter.next() 
headerline = headerline.replace("\n", "") 
headers = headerline.split("\t") 
headerout = [headers[0],headers[1],headers[2]]
for i in range(3, len(headers)):
	id = headers[i].split(".")[0]
	if id in dropSamples:
		dropColumns[i] = 1;
	else:
		headerout.append(headers[i])

# print new header
print "\t".join(headerout)
 
# read the rest of the report file lines, dropping columns present in dropColumns dict
for line in reportIter:
	line = line.replace("\n", "")
	fields = line.split("\t")
	out = [fields[0],fields[1],fields[2]]
	for i in range(3, len(fields)):
		if i not in dropColumns:
			out.append(fields[i])
	print "\t".join(out)

