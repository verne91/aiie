#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from aiie.main import aiie

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    input_parser = parser.add_argument_group("Input BAM file")
    input_parser.add_argument('-i', '--inputBam', nargs='+', required=True, metavar='INPUT_BAM', dest='input_bam_list', help='input bam file(s)')

    ref_parser = parser.add_argument_group("Reference files")
    ref_parser.add_argument('-r', '--ref', required=True, metavar='REFERENCE', dest='ref_file', help='referenece genome fasta file')
    ref_parser.add_argument('-t', '-targetRegion', required=True, metavar='TARGET', dest='target_file', help='Exome target region (bed file)')

    repeat_parser = parser.add_argument_group("Repeat files")
    repeat_parser.add_argument('--aluSeq', required=True, metavar='ALU_SEQ', dest='alu_seq_file', help='Alu sequence fasta file')
    repeat_parser.add_argument('--allRepeatSeq',required=True, metavar='ALL_REPEAT_SEQ', dest='all_repeat_file', help='All repeat sequence fasta file')
    repeat_parser.add_argument('--aluCoord', required=True, metavar='ALU_COORD', dest='alu_coord_file', help='ALU coordinates in genome (3-column bed file)')
    repeat_parser.add_argument('--simpleRepeatCoord', required=True, metavar='SIMPLE_REPEAT_COORD', dest='simple_repeat_coord_file', help='simple repeat coordinates in genome (3-column bed file)')

    cutoff_parser = parser.add_argument_group("Cutoff parameters")
    cutoff_parser.add_argument('--aluKmer', default=13, type=int, metavar='ALU_KMER', dest='alu_k', help='k-mer of Alu side (default: 13)')
    cutoff_parser.add_argument('--refKmer', default=13, type=int, metavar='REF_KMER', dest='ref_k', help='k-mer of reference side (default: 13)')
    cutoff_parser.add_argument('--extendLen', default=10, type=int, metavar='EXTENDED_LENGTH', dest='extend_len', help='extended length to check the match (default: 10)')
    cutoff_parser.add_argument('--allowMis', default=2, type=int, metavar='ALLOW_MISMATCH', dest='allow_mis', help='maximum allowing mismatch base pairs (default: 2)')
    cutoff_parser.add_argument('--minSupportRead', default=2, type=int, metavar='MIN_SUPPORT_READ', dest='min_support_read', help="minimum supported reads (default: 2")

    output_parser = parser.add_argument_group("Output options")
    output_parser.add_argument('-o', '--out', required=True, metavar='OUT', dest='out_file', help='output file')

    args = parser.parse_args()

    # print(args.allow_mis)
    sys.exit(aiie(args))
