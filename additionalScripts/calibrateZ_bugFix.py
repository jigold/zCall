#! /usr/bin/python

# zCall: A Rare Variant Caller for Array-based Genotyping
# Jackie Goldstein
# jigold@broadinstitute.org
# May 8th, 2012

import sys
from optparse import OptionParser

### Parse Inputs from Command Line
parser = OptionParser()
parser.add_option("-R","--report",type="string",dest="report",action="store",help="genome studio report to call from")
parser.add_option("-T","--thresholds",type="string",dest="thresholds",action="store",help="thresholds file")
(options, args) = parser.parse_args()

if options.report == None:
    print "specify GenomeStudio report file path with -R"
    sys.exit()

if options.thresholds == None:
    print "specify thresholds file with -T"
    sys.exit()

### Make a Dictionary with SnpName: (nAA, nBB) with all SNPs that meet MAF > 5%, min 10 points in each cluster, HWE p > 1e-5.
clusterCounts = {}
for line in open(options.report, 'r'):
    line = line.replace("\n", "")
    line = line.replace("\r", "")
    
    if line.find("Name") != -1: # skip header row
        continue

    fields = line.split("\t")
    snp = fields[0]
    genotypes = [fields[i] for i in range(3,len(fields),3)]
    nAA = genotypes.count("AA")
    nAB = genotypes.count("AB")
    nBB = genotypes.count("BB")
    nTotal = nAA + nAB + nBB

    # Calculate Missing Rate (ignore SNPs that have less than 99% call rate)
    if float(nTotal) / float(len(genotypes)) < 0.99:
        continue

    # Make sure there are at least 10 points in each homozygote cluster
    if nAA < 10 or nBB < 10:
        continue
    
    # Calculate MAF
    if nAA > nBB:
        maf = (nAB + 2 * nBB) / float(2 * nTotal)
    elif nAA <= nBB:
        maf = (nAB + 2 * nAA) / float(2 * nTotal)

    # MAF check ( >5% MAF)
    if maf < 0.05:
        continue
    
    # Hardy-Weinberg Equilibrium Check (don't use site if p_hwe < 0.00001)
    chiCritical = 19.5 # p = 0.00001 for 1 DOF

    if nAA > nBB:
        p = 1.0 - maf
        q = maf        
        expAA = p**2 * nTotal
        expAB = 2 * p * q * nTotal
        expBB = q**2 * nTotal

    if nBB >= nAA:
        p = 1.0 - maf
        q = maf        
        expAA = q**2 * nTotal
        expAB = 2 * p * q * nTotal
        expBB = p**2 * nTotal
        
    chiSquare = ((nAA - expAA)**2 / float(expAA)) + ((nAB - expAB)**2 / float(expAB)) + ((nBB - expBB)**2 / float(expBB))
    if chiSquare > chiCritical:
        continue


    clusterCounts[snp] = (nAA,nBB)
    

# Get thresholds for each site
thresholdsX = {}
thresholdsY = {}
for line in open(options.thresholds, 'r'):
    line = line.replace("\n", "")
    if line.find("SNP") != -1:
        continue
    else:
        fields = line.split("\t")
        snp = fields[0]

        if fields[1] != "NA":
            Tx = float(fields[1])
        else:
            Tx = fields[1]

        if fields[2] != "NA":
            Ty = float(fields[2])
        else:
            Ty = fields[2]
            
        thresholdsX[snp] = Tx
        thresholdsY[snp] = Ty


# Array for keeping track of counts for (GenCall, zCall)
counts = {("AA","AA"):0,("AA","AB"):0,("AA","BB"):0,("AA","NC"):0,
          ("AB","AA"):0,("AB","AB"):0,("AB","BB"):0,("AB","NC"):0,
          ("BB","AA"):0,("BB","AB"):0,("BB","BB"):0,("BB","NC"):0,
          ("NC","AA"):0,("NC","AB"):0,("NC","BB"):0,("NC","NC"):0}

# Iterate through Illumina GenomeStudio Report
for line in open(options.report, 'r'):
    line = line.replace("\n", "")
    line = line.replace("\r", "")

    fields = line.split("\t")

    if line.find("Name") != -1: # header row -- skip
        continue
    else:
        snp = fields[0]
        Tx = thresholdsX[snp]
        Ty = thresholdsY[snp]

        if Tx == "NA" or Ty == "NA": # skip site
            continue
        
        if snp not in clusterCounts: # Don't include sites that don't meet thresholding criteria
            continue
        nAA = clusterCounts[snp][0]
        nBB = clusterCounts[snp][1]

        for i in range(3, len(fields), 3):
            oldGT = fields[i]

            if nBB > nAA: # Make everything such that AA is the common homozygote cluster
                if oldGT == "AA":
                    oldGT = "BB"
                elif oldGT == "BB":
                    oldGT = "AA"
            
            x = float(fields[i + 1])
            y = float(fields[i + 2])

            if x >= Tx and y <= Ty:
                newGT = "AA"
                if nBB > nAA: # Make everything such that AA is the common homozygote cluster          
                    newGT = "BB"
            elif x <= Tx and y >= Ty:
                newGT = "BB"
                if nBB > nAA: # Make everything such that AA is the common homozygote cluster
                    newGT = "AA"
            elif x > Tx and y > Ty:
                newGT = "AB"
            elif x < Tx and y < Ty:
                newGT = "NC"

            counts[(oldGT,newGT)] += 1


## Output Concordance Statistics
total = sum([counts[i] for i in counts if i[0] != "NC"])

x = [("AA","AA"),("AA","AB"),("AA","BB"),("AA","NC"),
          ("AB","AA"),("AB","AB"),("AB","BB"),("AB","NC"),
          ("BB","AA"),("BB","AB"),("BB","BB"),("BB","NC"),
          ("NC","AA"),("NC","AB"),("NC","BB"),("NC","NC")]

c = [counts[i] for i in x]

globalConcordance = ((c[0] + c[5] + c[10]) / float(total)) * 100
specificity = (c[0] / float(c[0] + c[1] + c[2] + c[3])) * 100
sensitivityAB = (c[5] / float(c[4] + c[5] + c[6] + c[7])) * 100
sensitivityBB = (c[10] / float(c[8] + c[9] + c[10] + c[11])) * 100
npv = (c[0] / float(c[0] + c[4] + c[8] + c[12])) * 100
ppvAB = (c[5] / float(c[1] + c[5] + c[9] + c[13])) * 100
ppvBB = (c[10] / float(c[2] + c[6] + c[10] + c[14])) * 100

print "Thresholds used:", options.thresholds
print
print "Concordance Stats (%)"
print "Global Concordance:", globalConcordance
print "Specificity:", specificity
print "SensitivityAB:",sensitivityAB
print "SensitivityBB:",sensitivityBB
print "Negative Predictive Value:", npv
print "Positive Predictive Value AB:", ppvAB
print "Positive Predictive Value BB:", ppvBB
print
print "nAA:", (c[0] + c[1] + c[2] + c[3])
print "nAB:", (c[4] + c[5] + c[6] + c[7])
print "nBB:", (c[8] + c[9] + c[10] + c[11])
print "nNC:", (c[12] + c[13] + c[14] + c[15])
