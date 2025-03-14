#! /usr/bin/env python

"""
Author:
Matthew Ezewudo
DTBE LB CDC Atlanta

"""


if __name__ == '__main__':
    import sys
    import os
    import re
    import snp
    import argparse as ap
    import gzip
    import subprocess

    parser = ap.ArgumentParser(prog='Varpipeline', conflict_handler='resolve', 
                               description="Varpipeline - Call SNPs and InDels")
    group1 = parser.add_argument_group('Input', '')
    group1.add_argument('-q', '--fastq', required=True, help='Input FASTQ file', metavar="STRING")
    group1.add_argument('-r', '--reference', required=True, metavar="STRING", help='Reference genome in FASTA format.')
    group1.add_argument('-n', '--name', required=True, metavar="STRING", help='Sample name to be used as a prefix.')
    group1.add_argument('-q2', '--fastq2', default='', metavar="STRING", help='Second paired-end FASTQ file.')
    group2 = parser.add_argument_group('Output', '')
    group2.add_argument('-o', '--outdir', help='Output directory', metavar="STRING")
    group2.add_argument('--keepfiles', action='store_true', help='Keep intermediate files.')
    group3 = parser.add_argument_group('Aligners', 'Select a specific aligner.')
    group3.add_argument('--bwa', action='store_true', help='Align Illumina reads using bwa. (Default)')
    group4 = parser.add_argument_group('Callers', 'Choose program(s) to call SNPs/InDels with.')
    group4.add_argument('--all', action='store_true', help='Run all SNP / InDel calling programs.')
    group4.add_argument('--gatk', action='store_true', help='Run GATK SNP / InDel calling. (Default)')
    group4.add_argument('--samtools', action='store_true', help='Run SamTools SNP / InDel calling.')
    group5 = parser.add_argument_group('Annotation', 'Use snpEff to annotate VCF file')
    group5.add_argument('-a', action='store_true', help='Run snpEff functional annotation.')
    group6 = parser.add_argument_group('Optional', '')
    group6.add_argument('-v', '--verbose', action='store_true', help='Produce status updates of the run.') 
    group6.add_argument('-h', '--help', action='help', help='Show this help message and exit')
    group6.add_argument('--version', action='version', version='%(prog)s_amplicons v1.01', 
                        help='Show program\'s version number and exit')
   

    if len(sys.argv)==1:
        parser.print_usage()
        sys.exit(1)

    args = parser.parse_args()
    
    error = 0
    """ Validate input arguements. """
    # Verify input FASTQ file exists
    if not os.path.isfile(args.fastq):
        error += 1
        print("Please check if input '%s' exists, then try again." % args.fastq)
    
    # Verify reference file exists
    if not os.path.isfile(args.reference):
        error += 1
        print("Please check if reference '%s' exists, then try again." % args.reference)
    
    # Check paired end, and verify exists
    if args.fastq2:
       if not os.path.isfile(args.fastq2):
          error += 1
          print("Please check that '%s' exists, then try again." % args.fastq2)
       else:
            paired = True
    else:
        paired = False

    # Choose an aligner 
    aligners = 0
    if args.bwa: aligners += 1
    if not aligners: args.bwa = True 
    
    # Check for errors if so, print usage
    if error:
        print("")
        print("Use --help for more information.")
        parser.print_usage()
        sys.exit(2)
    
    # Choose a Caller
    if args.all:
        args.gatk = True
        args.samtools = True
    else:
        callers = 0
        if args.gatk: callers += 1
        if args.samtools: callers += 1   
        if not callers: args.gatk = True
    
    # If no outdir, use the given name as the outdir
    if not args.outdir:
        args.outdir = args.name

    if args.verbose:
        print(' '.join(sys.argv))
        print("")
        
    # All is well let's get started!
    s = snp.snp(args.fastq, args.outdir, args.reference, args.name, paired, args.fastq2, 
                        args.verbose, ' '.join(sys.argv))
    # run Clockwork decontamination
    #s.runClockwork()

    # run Trimmomatic timmer
    s.runTrimmomatic()

 
    # Run the aligner
    s.runBWA(args.bwa)

    # If asked, run SNP callers
    if args.gatk: s.runGATK()
    if args.samtools: s.runSamTools()
    
    
    # Annotate Final VCF
    s.annotateVCF()
   
    # Perform Coverage statistics Analysis
    s.runCoverage()

    # Print final reports
    s.runPrint()

    # By default clean up intermediate files
    if not args.keepfiles: s.cleanUp()
