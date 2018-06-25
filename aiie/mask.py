#!/usr/bin/env python
# -*- coding: utf-8 -*-

def mask_alu_region(chr_name, start, end, region_coords, ref_k):
    out_list = [start]
    flag = 0

    try:
        for coord in region_coords[chr_name]:
            if coord[0] < start:
                if coord[1] > start:
                    if coord[1] >= end:
                        return 0
                    else:
                        out_list[0] = coord[1] + ref_k
            else:
                if coord[0] < end:
                    if coord[1] < end:
                        out_list.append(coord[0] - ref_k)
                        out_list.append(coord[1] + ref_k)
                    else:
                        out_list.append(coord[0] - ref_k)
                        flag = 1
            if coord[0] > end:
                break
    except KeyError:
        flag = 0
    if flag == 0:
        out_list.append(end)
    b_list = []
    for i in range(int(len(out_list)/2)):
        b_list.append((out_list[2*i]-start, out_list[2*i+1]-start))
    return b_list


def mask_SR_region(chr_name, start, end, region_coords):
    out_list = [start]
    flag = 0

    try:
        for coord in region_coords[chr_name]:
            if coord[0] < start:
                if coord[1] > start:
                    if coord[1] >= end:
                        return 0
                    else:
                        out_list[0] = coord[1]
            else:
                if coord[0] < end:
                    if coord[1] < end:
                        out_list.append(coord[0])
                        out_list.append(coord[1])
                    else:
                        out_list.append(coord[0])
                        flag = 1
            if coord[0] > end:
                break
    except KeyError:
        flag = 0
    if flag == 0:
        out_list.append(end)
    b_list = []
    for i in range(int(len(out_list)/2)):
        b_list.append((out_list[2*i]-start, out_list[2*i+1]-start))
    return b_list

def intersect_regions(region1, region2):
    i = j = 0
    out_list = []
    while i < len(region1) and j < len(region2):
        if region1[i][0] > region2[j][1]:
            j += 1
        elif region1[i][1] < region2[j][0]:
            i += 1
        elif region1[i][0] <= region2[j][0]:
            if region1[i][1] >= region2[j][1]:
                out_list.append(region2[j][0])
                out_list.append(region2[j][1])
                j += 1
            else:
                out_list.append(region2[j][0])
                out_list.append(region1[i][1])
                i += 1
        elif region1[i][0] > region2[j][0]:
            if region1[i][1] <= region2[j][1]:
                out_list.append(region1[i][0])
                out_list.append(region1[i][1])
                i += 1
            else:
                out_list.append(region1[i][0])
                out_list.append(region2[j][1])
                j += 1
    b_list = []
    for i in range(int(len(out_list)/2)):
        b_list.append((out_list[2*i], out_list[2*i+1]))
    return b_list

def filter_SR(seq):
    num_A = num_T = num_C = num_G = 0
    for i in range(len(seq)):
        if seq[i] == "A":
            num_A += 1
        elif seq[i] == "T":
            num_T += 1
        elif seq[i] == "G":
            num_G += 1
        elif seq[i] == "C":
            num_C += 1
    if (num_A/len(seq) > 0.7) or (num_T/len(seq) > 0.7) or (num_G/len(seq) > 0.7) or (num_C/len(seq) > 0.7):
        return False
    else:
        return True

def filter_consecutiveAT(read):
    aa = "AAAAAAAAAAAAA"
    tt = "TTTTTTTTTTTTT"
    if aa in read:
        start = read.index(aa)
        after = read[start+13:]
        i = 0
        end = 0
        for base in after:
            if base != "A":
                end = i
                break
            i += 1
        if end == 0:
            end = len(read) - 1
        return (start, end+13+start)
    elif tt in read:
        start = read.index(tt)
        after = read[start+13:]
        i = 0
        end = 0
        for base in after:
            if base != "T":
                end = i
                break
            i += 1
        if end == 0:
            end = len(read) - 1
        return (start, end+13+start)
    else:
        return (0, 0)