#!/usr/bin/env python
# -*- coding: utf-8 -*-
from Bio import SeqIO

def bed2dict(bed_file):
    coords = {}
    with open(bed_file) as fin:
        for line in fin:
            tmp = line.split()
            try:
                coords[tmp[0]].append((int(tmp[1]), int(tmp[2])))
            except:
                coords[tmp[0]] = [(int(tmp[1]), int(tmp[2]))]
    for ch in coords:
        coords[ch].sort(key = lambda x:x[0])
    return coords

def kmer_index(fasta_file, k):
    a_dict = {}
    with open(fasta_file) as fin:
        for rec in SeqIO.parse(fin, "fasta"):
            myseq = str(rec.seq).upper()
            for i in range(len(myseq)-k):
                a_dict[myseq[i:i+k]] = 1