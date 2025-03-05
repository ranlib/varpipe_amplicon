#! /usr/bin/env python

""" Aceepts unmapped_reads,mapped_reads,target_coverage text files and sample name, generates stats output file """

"""
Author:
Matthew Ezewudo
DTBE LB CDC Atlanta

"""


import sys
import os
from string import join
from datetime import datetime

input1 = sys.argv[1]
input2 = sys.argv[2]
input3 = sys.argv[3]
input4 = sys.argv[4]
input5 = input4.split("-")[0]

(low_cov_count,unmapped,mapped,percent_mapped) = (0,0,0,0)

fh1 = open(input1)
for lines in fh1:
    unmapped = int(lines.rstrip("\r\n"))
fh1.close()

fh2 = open(input2)
for lines in fh2:
    mapped = int(lines.rstrip("\r\n"))
fh2.close()

fh3 = open(input3)
for lines in fh3:  
  if lines.startswith("SAMPLE_ID"):
       continue
  lined = lines.rstrip("\r\n").split("\t")
  if float(lined[6]) < 500.0:
     low_cov_count += 1
fh3.close()

percent_mapped = (float(mapped)/float(unmapped))*100.00
str_percent_mapped = "{0:.2f}".format(percent_mapped)
i = datetime.now()

print("Sample ID" + "\t" + "Sample Name" + "\t" + "Total reads from sequence file" + "\t" + "Number of mapped reads" + "\t" + "Percent mapped reads" + "\t" + "Low Coverage" + "\t" + "Pipeline Version" + "\t" + "Date\r")
print(input4 + "\t" + input5 + "\t" + str(unmapped) + "\t" +  str(mapped) + "\t" +  str_percent_mapped + "\t" + str(low_cov_count) + "\t" + "Varpipeline: Varpipe_amplicons_1.0.1" + "\t" + i.strftime('%Y/%m/%d %H:%M:%S') + "\r")

