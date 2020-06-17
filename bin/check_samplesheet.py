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
HEADER = ['sample', 'bam', 'genome', 'transcriptome', 'condition']

## CHECK HEADER
fin = open(args.DESIGN_FILE_IN,'r')
header = fin.readline().strip().split(',')
if header != HEADER:
    print("{} header: {} != {}".format(ERROR_STR,','.join(header),','.join(HEADER)))
    sys.exit(1)

outLines,condition_list = [],[]
while True:
    line = fin.readline()
    if line:
        lspl = [x.strip() for x in line.strip().split(',')]
        sample, bam, genome, transcriptome, condition = lspl

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

        ## CHECK GENOME ENTRIES
        if genome:
            if genome.find(' ') != -1:
                print("Genome entry contains spaces!",line)
                sys.exit(1)

            if len(genome.split('.')) > 1:
                if genome[-6:] != '.fasta' and genome[-3:] != '.fa' and genome[-9:] != '.fasta.gz' and genome[-6:] != '.fa.gz':
                    print("Genome entry does not have extension '.fasta', '.fa', '.fasta.gz' or '.fa.gz'!",line)
                    sys.exit(1)

        ## CHECK TRANSCRIPTOME ENTRIES
        gtf = ''
        if transcriptome:
            if transcriptome.find(' ') != -1:
                print("Transcriptome entry contains spaces!",line)
                sys.exit(1)

            if transcriptome[-6:] != '.fasta' and transcriptome[-3:] != '.fa' and transcriptome[-9:] != '.fasta.gz' and transcriptome[-6:] != '.fa.gz' and transcriptome[-4:] != '.gtf' and transcriptome[-7:] != '.gtf.gz':
                print("Transcriptome entry does not have extension '.fasta', '.fa', '.fasta.gz', '.fa.gz', '.gtf' or '.gtf.gz'!",line)
                sys.exit(1)

            if transcriptome[-4:] == '.gtf' or transcriptome[-7:] == '.gtf.gz':
                gtf = transcriptome
                if not genome:
                    print("If genome isn't provided, transcriptome must be in fasta format for mapping!",line)
                    sys.exit(1)
            
        ## CHECK CONDITION ENTRIES
        if condition:
            if condition not in condition_list:
                condition_list.append(condition)

        outLines.append([sample,bam,genome,transcriptome])
    else:
        fin.close()
        break

## WRITE TO FILE
fout = open(args.DESIGN_FILE_OUT,'w')
fout.write(','.join(['sample', 'bam', 'genome', 'transcriptome', 'is_transcripts']) + '\n')
for line in outLines:
    print(line)
    fout.write(','.join(line) + '\n')
fout.close()

outfile=open("total_conditions.csv","w")
outfile.write(",".join(condition_list)+"\n")
