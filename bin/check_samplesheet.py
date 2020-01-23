#!/usr/bin/env python

#######################################################################
#######################################################################
## Modified from nf-core/nanoseq script on 27th December 2019
#######################################################################
#######################################################################

import os
import sys
import argparse

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Reformat nf-core/nanornabam design file and check its contents.'
Epilog = """Example usage: python check_samplesheet.py <DESIGN_FILE_IN> <DESIGN_FILE_OUT>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('DESIGN_FILE_IN', help="Input design file.")
argParser.add_argument('DESIGN_FILE_OUT', help="Output design file.")

args = argParser.parse_args()

############################################
############################################
## MAIN SCRIPT
############################################
############################################

ERROR_STR = 'ERROR: Please check samplesheet'
HEADER = ['sample', 'bam', 'annotation']

## CHECK HEADER
fin = open(args.DESIGN_FILE_IN,'r')
header = fin.readline().strip().split(',')
if header != HEADER:
    print("{} header: {} != {}".format(ERROR_STR,','.join(header),','.join(HEADER)))
    sys.exit(1)

outLines = []
while True:
    line = fin.readline()
    if line:
        lspl = [x.strip() for x in line.strip().split(',')]
        sample,bam, annotation = lspl

        ## CHECK VALID NUMBER OF COLUMNS PER SAMPLE
        numCols = len([x for x in lspl if x])
        if numCols < 2:
            print("{}: Invalid number of columns (minimum of 2)!\nLine: '{}'".format(ERROR_STR,line.strip()))
            sys.exit(1)

        ## CHECK SAMPLE ID ENTRIES
        if sample:
            if sample.find(' ') != -1:
                print("{}: Sample ID contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)
        else:
            print("{}: Sample ID not specified!\nLine: '{}'".format(ERROR_STR,line.strip()))
            sys.exit(1)

        ## CHECK BAM ENTRIES
        if bam:
            if bam[-4:] != '.bam':
                print("{}: BAM file has incorrect extension (has to be '.bam')!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)
        else:
            print("{}: BAM file not specified!\nLine: '{}'".format(ERROR_STR,line.strip()))
            sys.exit(1)

        ## CHECK ANNOTATION ENTRIES
        if annotation:
            if annotation.find(' ') != -1:
                print("{}: Annotation field contains spaces!\nLine: '{}'".format(ERROR_STR,line.strip()))
                sys.exit(1)

            if len(annotation.split('.')) > 1:
                if annotation[-4:] != '.gtf' and annotation[-4:] != '.gff' and annotation[-7:] != '.gtf.gz' and annotation[-7:] != '.gff.gz':
                    print("{}: Annotation field has incorrect extension (has to be '.gtf', '.gff', '.gtf.gz' or '.gff.gz')!\nLine: '{}'".format(ERROR_STR,line.strip()))
                    sys.exit(1)

        else:
            print("{}: Annotation file not specified!\nLine: '{}'".format(ERROR_STR,line.strip()))
            sys.exit(1)

        outLines.append([sample,bam,annotation])
    else:
        fin.close()
        break

## WRITE TO FILE
fout = open(args.DESIGN_FILE_OUT,'w')
fout.write(','.join(['sample', 'bam', 'annotation']) + '\n')
for line in outLines:
    fout.write(','.join(line) + '\n')
fout.close()
