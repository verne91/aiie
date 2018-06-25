#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging
from Bio import SeqIO
from Bio.Seq import Seq

from aiie._version import __version__
from aiie.parser import *



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

    with open(args.alu_seq_file) as fin:
        for rec in SeqIO.parse(fin, "fasta"):
            myseq = str(rec.seq).upper()
            myseq_rev_comp = str(Seq(myseq).reverse_complement())
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
    ref_dict = []
    ref_seq_dict = {}
    whole_ref_dict = SeqIO.to_dict(SeqIO.parse(args.ref_file), "fasta")
    with open(args.target_file) as fin:
        for line in fin:
            tmp = line.strip().split()
            start = int(tmp[1])
            end = int(tmp[2])
            alu_splitted_region = mask_alu_region(tmp[0], start, end, alu_coords, ref_k)
            sr_splitted_region = mask_SR_region(tmp[0], start, end, sr_coords)

            myseq = str(whole_ref_dict[tmp[0]][start:end]).upper()
            
                