# __author__ = 'tianfuzneg'
# !/usr/bin/python
# -*- coding:utf-8 -*-

########################################################################################
# 20230331 add clipped ref loc; reverse loc  debug (ref_info_dic need to decide +-).
# 20241109 debug: remove flanking without human; remove unmmaped reads
########################################################################################
import os
import pysam
import itertools
from argparse import ArgumentParser


def get_names(names):
    with open(names, 'r') as infile:
        n = infile.read().splitlines()
    if '' in n:
        n.remove('')
    n_uniq = list(set(n))
    return n_uniq


def reverse_complement(dna):
    revc = ""
    basepair = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    for c in dna:
        revc = basepair[c] + revc
    return revc

def range_diff(r1, r2):
    s1, e1 = r1
    s2, e2 = r2
    endpoints = sorted((s1, s2, e1, e2))
    result = []
    if endpoints[0] == s1 and endpoints[1] != s1:
        result.append((endpoints[0], endpoints[1]))
    if endpoints[3] == e1 and endpoints[2] != e1:
        result.append((endpoints[2], endpoints[3]))
    return result

def multirange_diff(r1_list, r2_list):
    for r2 in r2_list:
        r1_list = list(itertools.chain(*[range_diff(r1, r2) for r1 in r1_list]))
    return r1_list

def extract_clipped_seq(options):
    jas_somatic_rnames = options.names
    jas_somatic_rnames_bam = options.bam
    n = get_names(jas_somatic_rnames)
    bamfile = pysam.AlignmentFile(jas_somatic_rnames_bam, 'rb')
    name_indexed = pysam.IndexedReads(bamfile)
    name_indexed.build()
    out_fa = options.out
    with open(out_fa, 'w') as out:
        for name in n:
            try:
                name_indexed.find(name)
            except KeyError:
                print('%s: extract read false' % name)
            else:
                iterator = name_indexed.find(name)
                # get alignment coordinate
                aligned_region_list = []
                ref_info_dic = {}
                query_length = 0
                query_sequence = 0
                for x in iterator:
                    if not (x.is_supplementary or x.is_secondary or x.is_unmapped):  # select primary, primary has full SEQ
                        ref_n = x.reference_name
                        cigartuples2 = x.cigartuples
                        aligned_pairsv2 = x.get_aligned_pairs(matches_only=False, with_seq=False)
                        reference_positions2 = x.get_reference_positions()
                        if x.is_reverse:
                            if cigartuples2[-1][0] in [4, 5]:  # may not exist clip
                                left_clip_len = cigartuples2[-1][1]
                            else:
                                left_clip_len = 0
                            if cigartuples2[0][0] in [4, 5]:
                                right_clip_len = cigartuples2[0][1]
                            else:
                                right_clip_len = 0
                            query_sequence = reverse_complement(
                                x.query_sequence)  # reverse SEQ in minimap2 was the same as fq
                        else:
                            if cigartuples2[0][0] in [4, 5]:
                                left_clip_len = cigartuples2[0][1]
                            else:
                                left_clip_len = 0
                            if cigartuples2[-1][0] in [4, 5]:
                                right_clip_len = cigartuples2[-1][1]
                            else:
                                right_clip_len = 0
                            query_sequence = x.query_sequence
                        clip_len_info = [left_clip_len, right_clip_len]
                        query_length = x.query_length
                        aligned_region_info = (left_clip_len, query_length - right_clip_len)
                        aligned_region_list.append(aligned_region_info)
                        if not x.is_reverse:
                            ref_info_dic[left_clip_len] = "%s_%s" % (
                                ref_n, str(reference_positions2[0]))  # clipped loc on ref
                            ref_info_dic[query_length - right_clip_len] = "%s_%s" % (
                            ref_n, str(reference_positions2[-1]))
                        else:
                            ref_info_dic[left_clip_len] = "%s_%s" % (
                                ref_n, str(reference_positions2[-1]))  # clipped loc on ref
                            ref_info_dic[query_length - right_clip_len] = "%s_%s" % (
                            ref_n, str(reference_positions2[0]))
                    elif x.is_supplementary:
                        ref_n = x.reference_name
                        cigartuples1 = x.cigartuples
                        aligned_pairsv1 = x.get_aligned_pairs(matches_only=False, with_seq=False)
                        reference_positions1 = x.get_reference_positions()
                        if x.is_reverse:  # left clip in reverse is right!
                            if cigartuples1[-1][0] in [4, 5]:
                                left_clip_len = cigartuples1[-1][1]
                            else:
                                left_clip_len = 0
                            if cigartuples1[0][0] in [4, 5]:
                                right_clip_len = cigartuples1[0][1]
                            else:
                                right_clip_len = 0
                        else:
                            if cigartuples1[0][0] in [4, 5]:
                                left_clip_len = cigartuples1[0][1]
                            else:
                                left_clip_len = 0
                            if cigartuples1[-1][0] in [4, 5]:
                                right_clip_len = cigartuples1[-1][1]
                            else:
                                right_clip_len = 0
                        clip_len_info = [left_clip_len, right_clip_len]
                        query_length = x.query_length + sum(clip_len_info)  # supplymentary was hard clip
                        aligned_region_info = (left_clip_len, query_length - right_clip_len)
                        aligned_region_list.append(aligned_region_info)
                        if not x.is_reverse:
                            ref_info_dic[left_clip_len] = "%s_%s" % (
                                ref_n, str(reference_positions1[0]))  # clipped loc on ref
                            ref_info_dic[query_length - right_clip_len] = "%s_%s" % (
                            ref_n, str(reference_positions1[-1]))
                        else:
                            ref_info_dic[left_clip_len] = "%s_%s" % (
                                ref_n, str(reference_positions1[-1]))  # clipped loc on ref
                            ref_info_dic[query_length - right_clip_len] = "%s_%s" % (
                            ref_n, str(reference_positions1[0]))
                    else:
                        continue
                        # print("alignment is secondary")
                # get clipped region
                if query_length == '0' or query_sequence == '0':
                    print("%s: extract BUG" % name)
                else:
                    full_region = [(0, query_length)]
                    clipped_region_list = multirange_diff(full_region, aligned_region_list)
                    # get clipped sequence
                    for clip_region in clipped_region_list:
                        if clip_region[1] - clip_region[0] >= 30:  # only output clipped seq >= 30
                            clip_seq = query_sequence[clip_region[0]:clip_region[1]]
                            try:
                                ref_info_dic[clip_region[0]]
                            except KeyError:
                                clipped_left_loc = 'N'
                            else:
                                clipped_left_loc = ref_info_dic[clip_region[0]]
                            try:
                                ref_info_dic[clip_region[1]]
                            except KeyError:
                                clipped_right_loc = 'N'
                            else:
                                clipped_right_loc = ref_info_dic[clip_region[1]]
                            clipped_ref_loc = '-'.join([clipped_left_loc, clipped_right_loc])
                            clip_seq_id = ">%s:%s_%s_%s_%s" % (
                                name, str(clip_region[0]), str(clip_region[1]), str(len(clip_seq)), clipped_ref_loc)
                            if 'chr' in clipped_ref_loc:  # 20241109 debug
                                out.write(clip_seq_id + '\n')
                                out.write(clip_seq + '\n')

if __name__ == "__main__":
    parser = ArgumentParser(description='Extract clipped seq by read name from bam file')
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-n', '--names', help='list of read names to extract', required=True)
    parser.add_argument('-o', '--out', help='file name for extracted sequence', required=True)
    options = parser.parse_args()
    extract_clipped_seq(options)
