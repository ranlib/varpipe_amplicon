#! /usr/bin/env python

""" 
Author

Matthew Ezewudo

DTBE LB CDC Atlanta

"""

import sys

""" The scripts accepts the stats, target coverage and final annotation file """
""" It parses and merge those to the summary output file """

input1 = sys.argv[1]
input2 = sys.argv[2]
input3 = sys.argv[3]
matrix = []
flag = ""

annot  = ['760314','760315','760316','761277','761278','761279','4247730','4247731','2715342',
          '1674046', '1674047','1674048','4247429','4247430','4247431','6575','6576','6577','2715369',
          '6620','6621','6622','1473246','1473247','1473329','1475927','1476471','2715344','2715346',
          '4247573','4247574','4247575','4247729','4247468', '4247469', '4247470']

print "Sample Summary:"

fh1 = open(input1,'r')
for lines in fh1:
    lined = lines.rstrip("\r\n").split("\t")
    matrix.append(lined)
print matrix[0][0] + ":" + "\t" + matrix[1][0]
print matrix[0][1] + ":" + "\t" + matrix[1][1]
print matrix[0][5] + ":" + "\t" + matrix[1][5]
print matrix[0][6] + ":" + "\t" + matrix[1][6]
print matrix[0][7] + ":" + "\t" + matrix[1][7]
fh1.close()
    
print "\n"
print "Amplicon Target Summary:"

fh2 = open(input2,'r')
for lines in fh2:
    lined = lines.rstrip("\r\n").split("\t")
    print lined[5] + "\t" + lined[2] + "\t" + lined[3] + "\t" + lined[6] + "\t" + lined[7] + "\t" + lined[8]
fh2.close()

print "\n"
print "Variant Summary:"

fh3 = open(input3,'r')
print "POS" + "\t" + "Gene Name" + "\t" + "Nucleotide Change" + "\t" + "Amino acid Change" + "\t" + "Read Depth" + "\t" + "Percent Alt Allele" + "\t" + "Annotation" + "\t" + "Flag"

for lines in fh3:
    if lines.startswith("Sample ID"):
       continue
    lined = lines.rstrip("\r\n").split("\t")
    if lined[2] in annot:
       flag = "Y"
    elif 761081 < int(lined[2]) and int(lined[2]) < 761163:
         flag = "Y"
    elif 6733 < int(lined[2]) and int(lined[2]) < 6743:
         flag = "Y"
    elif 7562 < int(lined[2]) and int(lined[2]) < 7584:
         flag = "Y"
    elif 2153878 < int(lined[2]) and int(lined[2]) < 2156149:
         flag = "Y"
    elif 2288681 < int(lined[2]) and int(lined[2]) < 2289282:
         flag = "Y"
    elif 778905 < int(lined[2]) and int(lined[2]) < 779498:
         flag = "Y"
    elif 2859290 < int(lined[2]) and int(lined[2]) < 2860452:
         flag = "Y"
    elif 1460996 < int(lined[2]) and int(lined[2]) < 1461301:    
         flag = "Y"
    elif 800792 < int(lined[2]) and int(lined[2]) < 801462:
         flag = "Y"
    elif 1673422 < int(lined[2]) and int(lined[2]) < 1673434:
         flag = "Y"
    else:
         flag = ""
    print lined[2] + "\t" + lined[15] + "\t" + lined[9] + "\t" + lined[11] + "\t" + lined[5] + "\t" + lined[6]  + "\t" + lined[7] + "\t" + flag
fh3.close()

