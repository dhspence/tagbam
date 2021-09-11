from __future__ import division
import sys
import csv
import argparse
import pysam
import string
import math
import pyranges as pr
import pandas as pd

pd.options.mode.chained_assignment = None

parser = argparse.ArgumentParser(description='Add amplicon data to VCF file with variants identified from HaloplexHS data.')
parser.add_argument('bamin',help='Consensus BAM file')
parser.add_argument('ampliconBed',help='Consensus BAM file')
parser.add_argument('bamout',help='Consensus BAM file')

parser.add_argument('-d',"--maxDist",default=5,help='maximum sum of the distances between amplicon and fragment ends')

# parse arguments
args = parser.parse_args()

inbam = args.bamin
outbam = args.bamout
ampfile = args.ampliconBed # must be in standard bed format and have strand info

# maximum sum of the distances between amplicon and read fragment ends to assign a read pair to an amplicon
maxDist = args.maxDist

####################################
#
# Main script
#
####################################

# open bam file
bamin = pysam.AlignmentFile(inbam,"rb")
# create out bam file
bamout = pysam.AlignmentFile(outbam,"wb",template=bamin)

# read in bed file and make a pd object--seems to be faster than pyranges
amppr = pr.read_bed(ampfile).df

# get the first chromosome and stratify the search by chrom to try and speed things up
ampdf = amppr[(amppr["Chromosome"]==amppr["Chromosome"].iloc[0])]
# current chrom
chrom = ampdf["Chromosome"].iloc[0]

# keep track of read names--if amplicon has already been assigned to a read pair then no need to do the lookup for the second read 
reads = {}

# progress counter
readnum = 0

# get all reads
for read in bamin.fetch():
    readnum += 1

    if readnum % 10000 == 0:
        print("number of reads processed: " + str(readnum),file=sys.stderr)

    # if one read has been assigned then assign the same to the other one
    if read.query_name in reads.keys():
        read.set_tag("XN",reads[read.query_name],replace=True)

    # has to be a proper pair with both reads mapped
    elif  read.is_proper_pair and not read.is_unmapped and not read.mate_is_unmapped and read.template_length is not None:
        mymatch = '';
        refseq = read.reference_name
        refstart = 0
        refend = 0
        probestrand = ''

        # speed things (maybe?) up by subsetting the pd object by chromosome and resetting the reads dict
        if refseq != chrom:
            ampdf = amppr[(amppr["Chromosome"]==refseq)]
            chrom = refseq
            reads = {}
            print("Now on chromosome " + str(chrom),file=sys.stderr)
            
        # get probestrand--if read 1 is not reverse then its +, otherwise -. Probably doesnt matter as long as these are disting. from each other
        if read.is_read1 and not read.is_reverse:
            probestrand = '+'
        else:
            probestrand = '-'

        # calc fragment start and end, depending on orientation
        if not read.is_reverse:
            refstart = read.reference_start
            refend = read.reference_start + read.template_length - 1
        else:
            refstart = read.next_reference_start
            refend = read.next_reference_start + abs(read.template_length) - 1

        # lookup the amplicon by chrom, start, end, and strand
        ov = ampdf[(ampdf["End"]<=(refend+maxDist)) & (ampdf["Start"]>=(refstart-maxDist)) & (ampdf["Strand"]==probestrand)]

        # if there are hits...
        if ov.shape[0] > 0:

            # get sum of distance between amplicon and read pair ends
            ov['distance'] = abs(ov['End'] - refend) + abs(ov['Start'] - refstart)
            ov = ov.sort_values(by=['distance'])

            # only assign if its less than maxDist
            if ov["distance"].iloc[0] < maxDist:
                read.set_tag("XN",ov["Name"].iloc[0],replace=True)
                reads[read.query_name] = ov["Name"].iloc[0]
        else:
            reads[read.query_name] = ''

    else:
        reads[read.query_name] = ''

    # write bam file
    bamout.write(read)

bamin.close()
bamout.close()

