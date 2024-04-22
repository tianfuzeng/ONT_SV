# __author__ = 'tianfuzneg'
# !/usr/bin/python
# -*- coding:utf-8 -*-

########################################################################################
# 20230403
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


def extract_ins_seq(options):
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
                for x in iterator:
                    if not x.is_secondary:  # secondary aligment didn't have SEQ info
                        query_n = x.query_name
                        ref_n = x.reference_name
                        # get query sequence
                        query_sequence = x.query_sequence
                        # get INS loc on query
                        ins_loc_list = []
                        cigartuples = x.cigartuples
                        aligned_pairsv = x.get_aligned_pairs(matches_only=False, with_seq=False)
                        cigar_loc = 0  # coordinate in cigar
                        query_loc = 0  # coordinate in query
                        for cigar in cigartuples:
                            if cigar[0] != 2 and cigar[0] != 5:
                                # DEL, # supplementary is hard clip in minimap2, query_loc should exclude
                                query_loc += cigar[1]
                            if cigar[0] != 5:  # supplementary is hard clip
                                cigar_loc += cigar[1]
                            if cigar[0] == 1 and cigar[1] >= 30:  # INS length >= 30bp
                                query_sv_end = query_loc
                                query_sv_start = query_loc - cigar[1]
                                align_sv_end = cigar_loc
                                align_sv_start = cigar_loc - cigar[1]
                                sv_len = cigar[1]
                                sv_ref_start = aligned_pairsv[align_sv_start - 1][1] + 1
                                sv_ref_end = aligned_pairsv[align_sv_end][1] + 1
                                out_list = [query_n, query_sv_start, query_sv_end, sv_len, ref_n,
                                            sv_ref_start, sv_ref_end, align_sv_start, align_sv_end]
                                ins_loc_list.append(out_list)
                        # get INS seq
                        if len(ins_loc_list) != '0':
                            for ins_loc in ins_loc_list:
                                # print(ins_loc_list)
                                seq_ins = query_sequence[ins_loc[1]:ins_loc[2]]
                                ins_seq_id = ">%s:%s" % (name, '_'.join(list(map(str, ins_loc[1:7]))))
                                out.write(ins_seq_id + '\n')
                                out.write(seq_ins + '\n')
                                out.flush()
                        else:
                            print('%s: none INS' % name)
                    else:
                        continue
                        # print("alignment is secondary")

if __name__ == "__main__":
    parser = ArgumentParser(description='Extract INS seq by read name from bam file')
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-n', '--names', help='list of read names to extract', required=True)
    parser.add_argument('-o', '--out', help='file name for extracted sequence', required=True)
    options = parser.parse_args()
    extract_ins_seq(options)