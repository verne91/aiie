# **aiie**


----------

Description
-----------
This software is designed to identify the active ALU insertion in the exon regions.

Required Resources
------------------
    Python:		https://www.python.org
    pysam:		http://pysam.readthedocs.io/en/latest/
    Biopython:  https://biopython.org/

Quick Start
-----------
Download and install:
    
    git clone git@github.com:verne91/aiie.git

Index alignments:

    samtools index <in.bam>

Run aiie with default parameters:

    python run_aiie.py -i INPUT_BAM [INPUT_BAM ...] -r REFERENCE -t TARGET
                   --aluSeq ALU_SEQ --allRepeatSeq ALL_REPEAT_SEQ --aluCoord
                   ALU_COORD --simpleRepeatCoord SIMPLE_REPEAT_COORD
                   [--aluKmer ALU_KMER] [--refKmer REF_KMER]
                   [--extendLen EXTENDED_LENGTH] [--allowMis ALLOW_MISMATCH]
                   [--minSupportRead MIN_SUPPORT_READ] -o OUT_FILE

Usage
-----
    python run_aiie.py [parameters]
    

**Parameters:**
optional arguments:

    -h, --help                              show this help message and exit

Input BAM file:

    -i INPUT_BAM [INPUT_BAM ...],
    --inputBam INPUT_BAM [INPUT_BAM ...]    input bam file(s)

Reference files:

    -r REFERENCE, --ref REFERENCE           referenece genome fasta file
    -t TARGET, -targetRegion TARGET         Exome target region (bed file)

Repeat files:

    --aluSeq ALU_SEQ                        Alu sequence fasta file
    --allRepeatSeq ALL_REPEAT_SEQ           All repeat sequence fasta file
    --aluCoord ALU_COORD                    ALU coordinates in genome (3-column bed file)
    --simpleRepeatCoord SIMPLE_REPEAT_COORD simple repeat coordinates in genome (3-column bed file)

Cutoff parameters:

    --aluKmer ALU_KMER                      k-mer of Alu side (default: 13)
    --refKmer REF_KMER                      k-mer of reference side (default: 13)
    --extendLen EXTENDED_LENGTH             extended length to check the match (default: 10)
    --allowMis ALLOW_MISMATCH               maximum allowing mismatch base pairs (default: 2)
    --minSupportRead MIN_SUPPORT_READ       minimum supported reads (default: 2

Output options:

    -o OUT, --out OUT                       output file
