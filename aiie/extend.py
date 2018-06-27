#!/usr/bin/env python
# -*- coding: utf-8 -*-

def right_extend(ref_seq, read_seq, ref_pos, read_pos, ref_k, km, l=10, nm=3):
    num_mismatch = 0
    if read_pos+km+l > len(read_seq):
        residual_len = len(read_seq) - read_pos - km
        if residual_len > 4 and read_seq[read_pos+km: read_pos+km+residual_len] == ref_seq[ref_pos+ref_k:ref_pos+ref_k+residual_len]:
            return True
        else:
            return False
        # for i in range(residual_len):
        #     if read_seq[read_pos+km+i] != ref_seq[ref_pos+ref_k+i]:
        #         return False
        # return True
    elif ref_pos+ref_k+l > len(ref_seq):
        residual_len = len(ref_seq) - ref_pos - ref_k
        if residual_len > 4 and read_seq[read_pos+km: read_pos+km+residual_len] == ref_seq[ref_pos+ref_k:ref_pos+ref_k+residual_len]:
            return True
        else:
            return False
        # for i in range(residual_len):
        #     if read_seq[read_pos+km+i] != ref_seq[ref_pos+ref_k+i]:
        #         return False
        # return True
    else:
        for i in range(l):
            if read_seq[read_pos+km+i] != ref_seq[ref_pos+ref_k+i]:
                num_mismatch += 1
                if num_mismatch > nm:
                    return False
        return True

def left_extend(ref_seq, read_seq, ref_pos, read_pos, l=10, nm=3):
    num_mismatch = 0
    if read_pos-l < 0:
        residual_len = read_pos
        if residual_len > 4 and read_seq[read_pos-residual_len: read_pos] == ref_seq[ref_pos-residual_len:ref_pos]:
            return True
        else:
            return False
        # for i in range(residual_len):
        #     if read_seq[read_pos-i] != ref_seq[ref_pos-i]:
        #         return False
        # return True
    elif ref_pos-l < 0:
        residual_len = ref_pos
        if residual_len > 4 and read_seq[read_pos-residual_len: read_pos] == ref_seq[ref_pos-residual_len:ref_pos]:
            return True
        else:
            return False
        # for i in range(residual_len):
        #     if read_seq[read_pos-i] != ref_seq[ref_pos-i]:
        #         return False
        # return True
    else:
        for i in range(1, l+1):
            if read_seq[read_pos-i] != ref_seq[ref_pos-i]:
                num_mismatch += 1
                if num_mismatch > nm:
                    return False
        return True

def rev_comp_left_extend(rev_ref_seq, read_seq, rev_ref_pos, read_pos, k, l=10, nm=3):
    num_mismatch = 0
    if read_pos-l < 0:
        residual_len = read_pos
        if residual_len > 4 and read_seq[read_pos-residual_len: read_pos] == rev_ref_seq[rev_ref_pos-k-residual_len:rev_ref_pos-k]:
            return True
        else:
            return False
    else:
        for i in range(1, l+1):
            if read_seq[read_pos-i] != rev_ref_seq[rev_ref_pos-k-i]:
                num_mismatch += 1
                if num_mismatch > nm:
                    return False
        return True

def rev_comp_right_extend(rev_ref_seq, read_seq, rev_ref_pos, read_pos, l=10, nm=3):
    num_mismatch = 0
    read_len = len(read_seq)
    if read_pos+l > read_len:
        residual_len = read_len - read_pos
        if residual_len > 4 and read_seq[read_pos: read_pos+residual_len] == rev_ref_seq[rev_ref_pos:rev_ref_pos+residual_len]:
            return True
        else:
            return False
    else:
        for i in range(1, l+1):
            if read_seq[read_pos+i] != rev_ref_seq[rev_ref_pos+i]:
                num_mismatch += 1
                if num_mismatch > nm:
                    return False
        return True

def checkmate(read, ref_match):
    tmp = ref_match.split(":")
    ref_chr = tmp[0]
    ref_pos = int(tmp[1])
    mate_chr = read.next_reference_name
    mate_pos = read.next_reference_start
    if mate_chr == ref_chr and abs(mate_pos-ref_pos) < 1000:
        return True
    else:
        return False
