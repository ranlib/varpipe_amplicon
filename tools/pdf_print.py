#! /usr/bin/env python

""" Accepts summary text file and prints out PDF report """

"""
Author:
Matthew Ezewudo
DTBE LB CDC Atlanta

"""


from fpdf import FPDF
import csv
#from string import join
import sys

class PDF(FPDF):
    def header(self):
        input1 = sys.argv[1] 
        f = open(input1, "r")
        f_reader = csv.reader(f, delimiter='\t')
        k = []
        for z in f_reader:
            x = '\t'.join(z) 
            if len(z) > 0 and 'Sample Summary' not in x and 'Date' not in x:
               k.append(z)
            if 'Date' in x:
               k.append(z)
               break
        self.set_font('Arial', '', 9)
        self.ln(0.5)
        for r in k:
            for d in r:
                self.cell(30, 5, str(d), border=0)
            self.ln(5)
        f.close()
        # Line break
        self.ln(10)
        b = [['POS', 'Gene Name', 'Nucleotide Change', 'Amino acid Change', 'Read Depth', 'Percent Alt allele', 'Annotation', 'Flag'],[]]
        if self.page_no() > 2:
           self.set_font('Arial', '', 9)
           self.ln(0.5)
           for c in b:
               if len(c) > 1:
                  for d in c:
                      self.cell(38, 5, str(d), border=0)
                  self.ln(5)
    # Page footer
    def footer(self):
        # Position at 1.5 cm from bottom
        self.set_y(-15)
        # Arial italic 8
        self.set_font('Arial', 'I', 8)
        # Page number
        self.cell(0, 10, 'Page ' + str(self.page_no()) + ' of {nb}', 0, 0, 'C')


pdf = PDF('L', 'mm', 'A4')
pdf.alias_nb_pages()
pdf.add_page()
pdf.set_font('Arial', '', 10)

input1 = sys.argv[1]
input2 = sys.argv[2]

f2 = open(input1, 'r')
f2_reader = csv.reader(f2, delimiter='\t')
p = []
for x in f2_reader:
    y = '\t'.join(x)
    if 'Sample' in y or 'SAMPLE ID' in y or  'Low' in y or 'Pipeline' in y or 'Date' in y:
       continue
    if 'Variant' in y:
       break
    if 'Amplicon' in y:
       pdf.set_font('Arial', 'B', 9)
       pdf.cell(0, 5, y, 0, 1)
       pdf.set_font('Arial', '', 9)
    elif len(x) > 0:
       p.append(x)
f2.close()

pdf.ln(0.5)
for row in p:
    for datum in row:
       pdf.cell(40, 5, str(datum), border=0)
    pdf.ln(5)
pdf.ln(160)
f3 = open(input1, 'r')
f3_reader = csv.reader(f3, delimiter='\t')
m = []
pdf.set_font('Arial', 'B', 9)
pdf.cell(0, 5, 'Variant Summary:', 0, 1)
lines = False
for y in f3_reader:
    n = '\t'.join(y)
    if 'Variant' in n:
       lines = True
       continue
    if lines == True and len(n) > 0:
       m.append(y)
f3.close()

pdf.set_font('Arial', '', 9)
pdf.ln(0.5)
for row in m:
    for datum in row:
        if len(datum) > 24 and 'c.' in datum:
           datum = datum[:25] + '*'
        if len(datum) > 18:
           pdf.set_font('Arial', '', 8)
           pdf.cell(38, 5, str(datum), border=0)
        else:
           pdf.set_font('Arial', '', 9)
           pdf.cell(38, 5, str(datum), border=0)
    pdf.ln(5)
pdf.output(input2, 'F')

