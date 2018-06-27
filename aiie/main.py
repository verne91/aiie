#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
from Bio import SeqIO
from Bio.Seq import Seq
import pysam
import gc

from aiie._version import __version__
from aiie.parser import *
from aiie.extend import *
from aiie.mask import *



FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

def aiie(args):
    logger.info("Running aiie %s" % __version__)
    logger.info("Command-line %s" % (" ".join(sys.argv)))
    logger.info("Arguments are " + str(args))

    # k-mer
    alu_k = args.alu_k
    ref_k = args.ref_k
    km = alu_k + ref_k
    extend_len = args.extend_len
    allow_mis = args.allow_mis

    #reference alu coordinates
    logger.info("alu coordinates")
    alu_coords = bed2dict(args.alu_coord_file)

    #reference simple repeat coordinates
    logger.info("simple repeat coordinates")
    sr_coords = bed2dict(args.simple_repeat_coord_file)
   
    #all repeat k-mer index
    logger.info("all repeat k-mer index")
    all_repeat_dict = kmer_index(args.all_repeat_file, ref_k)

    #build k-mer index for alu
    logger.info("build k-mer index for alu")
    alu_dict = {}
    revcomp_alu_dict = {}
    long_alu_dict = {}
    alu_len = 0

    with open(args.alu_seq_file) as fin:
        for rec in SeqIO.parse(fin, "fasta"):
            myseq = str(rec.seq).upper()
            myseq_rev_comp = str(Seq(myseq).reverse_complement())
            alu_len = len(myseq)
            alu_seq = myseq
            alu_seq_revcomp = myseq_rev_comp
            for i in range(len(myseq)-alu_k):
                try:
                    alu_dict[myseq[i:i+alu_k]].append(("plus",i))
                except KeyError:
                    alu_dict[myseq[i:i+alu_k]] = [("plus",i)]
                revcomp_kmer = str(Seq(myseq[i:i+alu_k]).reverse_complement())
                try:
                    revcomp_alu_dict[revcomp_kmer].append(("neg",i))
                except KeyError:
                    revcomp_alu_dict[revcomp_kmer] = [("neg",i)]
            for i in range(len(myseq)-km):
                long_alu_dict[myseq[i:i+km]] = 1
                long_alu_dict[myseq_rev_comp[i:i+km]] = 1
    
    
    #build k-mer index for ref
    logger.info("build k-mer index for ref")
    ref_dict = {}
    ref_seq_dict = {}
    long_ref_dict = {}
    whole_ref_dict = SeqIO.to_dict(SeqIO.parse(args.ref_file, "fasta"))
    with open(args.target_file) as fin:
        for line in fin:
            tmp = line.strip().split()
            start = int(tmp[1])
            end = int(tmp[2])
            myseq = str(whole_ref_dict[tmp[0]][start:end].seq).upper()
            seq_id = tmp[0]+":"+tmp[1]
            ref_seq_dict[seq_id] = myseq
            alu_splitted_region = mask_alu_region(tmp[0], start, end, alu_coords, ref_k)
            sr_splitted_region = mask_SR_region(tmp[0], start, end, sr_coords)
            if alu_splitted_region == 0 or sr_splitted_region == 0:
                continue
            else:
                splitted_region = intersect_regions(alu_splitted_region, sr_splitted_region)
                for coords in splitted_region:
                    if coords[1] - coords[0] >= ref_k:
                        for i in range(coords[0],coords[1]-ref_k):
                            if (myseq[i:i+ref_k] not in alu_dict) and (myseq[i:i+ref_k] not in revcomp_alu_dict):
                                ##test
                                ##print(myseq[i:i+ref_k])
                                if filter_SR(myseq[i:i+ref_k]) and filter_all_repeat(myseq, i, ref_k, all_repeat_dict):
                                    try:
                                        ref_dict[myseq[i:i+ref_k]].append((seq_id,i))
                                    except KeyError:
                                        ref_dict[myseq[i:i+ref_k]] = [(seq_id,i)]
            for i in range(len(myseq)-km):
                long_ref_dict[myseq[i:i+km]] = 1
    del whole_ref_dict
    gc.collect()
    #extend k-mer
    logger.info("extend k-mer")
    
    extend_results = []
    for bam_file in args.input_bam_list:
        samfile = pysam.AlignmentFile(bam_file, "rb")
        logger.info("begin "+bam_file)
        for read in samfile.fetch():
            if read.get_cigar_stats()[0][0]/read.query_length < 0.98 and read.get_cigar_stats()[0][10] <= 5:
                myseq = read.query_sequence.upper()
                omit_region = filter_consecutiveAT(myseq)

                for i in range(len(myseq)-km):
                    if i+km >= omit_region[0] and i+km < omit_region[1]:
                        continue
                    elif i >= omit_region[0] and i < omit_region[1]:
                        continue

                    if (myseq[i:i+km] not in long_alu_dict) and (myseq[i:i+km] not in long_ref_dict):
                        if (myseq[i:i+ref_k] in ref_dict) and (myseq[i+ref_k:i+km] in alu_dict):
                            ref_match_list = ref_dict[myseq[i:i+ref_k]]
                            alu_match_list = alu_dict[myseq[i+ref_k:i+km]]
                            for alu_match in alu_match_list:
                                for ref_match in ref_match_list:
                                    left = left_extend(ref_seq_dict[ref_match[0]], myseq, ref_match[1], i, extend_len, allow_mis)
                                    right = right_extend(alu_seq, myseq, alu_match[1], i, alu_k, km, extend_len, allow_mis)
                                    if left and right:
                                        if checkmate(read, ref_match[0]):
                                            # out_list = [ref_match[0], str(ref_match[1]+ref_k), alu_match[0], str(alu_match[1]), myseq[i:i+km], read.query_name, str(int(left)), str(int(right)), str(int(left and right)), "1"]
                                            # fout.write("\t".join(out_list)+"\n")
                                            out_list = [ref_match[0], ref_match[1]+ref_k, alu_match[0], alu_match[1], read.query_name]
                                            extend_results.append(out_list)
                        elif (myseq[i:i+alu_k] in revcomp_alu_dict) and (myseq[i+alu_k:i+km] in ref_dict):
                            rev_alu_match_list = revcomp_alu_dict[myseq[i:i+alu_k]]
                            ref_match_list = ref_dict[myseq[i+alu_k:i+km]]
                            for alu_match in rev_alu_match_list:
                                for ref_match in ref_match_list:
                                    left = rev_comp_left_extend(alu_seq_revcomp, myseq, alu_len-alu_match[1], i, alu_k, extend_len, allow_mis)
                                    right = right_extend(ref_seq_dict[ref_match[0]], myseq, ref_match[1], i, ref_k, km, extend_len, allow_mis)
                                    if left and right:
                                        if checkmate(read, ref_match[0]):
                                            # out_list = [ref_match[0], str(ref_match[1]), alu_match[0], str(alu_match[1]), myseq[i:i+km], read.query_name, str(int(left)), str(int(right)), str(int(left and right)), "2"]
                                            # fout.write("\t".join(out_list)+"\n")
                                            out_list = [ref_match[0], ref_match[1], alu_match[0], alu_match[1], read.query_name]
                                            extend_results.append(out_list)

        samfile.close()
    
    #cutoff by supported reads number
    logger.info("process the result")
    summary_dict = {}
    for record in extend_results:
        if record[0] in summary_dict:
            flag = 0
            for pos in summary_dict[record[0]]:
                if abs(int(record[1])-pos[0]) == abs(int(record[3])-pos[1]):
                    if record[4] not in summary_dict[record[0]][pos]:
                        summary_dict[record[0]][pos].append(record[4])
                    flag = 1
                    break
            if flag == 0:
                summary_dict[record[0]][(int(record[1]), int(record[3]))] = [record[4]]
        else:
            summary_dict[record[0]] = {}
            summary_dict[record[0]][(int(record[1]), int(record[3]))] = [record[4]]

    cutoff = args.min_support_read

    fout = open(args.out_file, "w")
    for k in summary_dict:
        for p in summary_dict[k]:
            if len(summary_dict[k][p]) >= cutoff:
                tmp = k.split(":")
                fout.write(tmp[0]+"\t"+str(int(tmp[1])+p[0])+"\t"+str(p[1])+"\t"+str(len(summary_dict[k][p]))+"\n")

    fout.close()

    logger.info("FINISH!")