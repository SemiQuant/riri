#!/bin/python3
import argparse
from plastid import *
parser = argparse.ArgumentParser(description="Make 5' UTRs")
parser.add_argument('--gtf_in', dest='gtf_in', type=str, help='GTF input')
args = parser.parse_args()

gtf_out = args.gtf_in + "_5p.gtf"

new_gtf = open(gtf_out,"w")
upstream_utr_length = 30

ORFS = GTF2_TranscriptAssembler(open(args.gtf_in))
for my_orf in ORFS:
    span = my_orf.spanning_segment
    if my_orf.strand == "+":
        new_region = GenomicSegment(span.chrom,span.start - upstream_utr_length,span.end,span.strand)
    else:
        new_region = GenomicSegment(span.chrom,span.start,span.end+upstream_utr_length,span.strand)  
    new_transcript = Transcript(new_region,**my_orf.attr) # copy metadata attributes from old ORF
    new_gtf.write(new_transcript.as_gtf())

new_gtf.close()