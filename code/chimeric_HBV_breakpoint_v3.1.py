# __author__ = 'tianfuzneg'
# !/usr/bin/python
# -*- coding:utf-8 -*-

import os
import pysam
import multiprocessing.pool
import datetime
# from multiprocessing import Pool
# from multiprocessing.pool import ThreadPool as Pool

########################################################################################
# 20230926
########################################################################################
samtools = "/data/fs01/biosoft/samtools-1.9/samtools"
extract_py = "/data/fs09/wangzf/nanopore/ztf/HCC/ONT/program/extract_reads_from_bam_by_ids_v2.py"


### extract HBV bam
# extract reads on HBV genome, get readsid
# if reads contain HBV seq, extract all reads alignments results (bam)
def extract_hbv_readsid(hbv_id_txt, bam_sort, bam_sort_hbv, hbv_reads_id_txt, bam_hbv_all):
    hbv_newid_list = os.popen("cat %s" % hbv_id_txt).read().strip().split()
    os.system("{samtools} view -@ 5 -hb {bam_in} {str} > {bam_out}".format(
        samtools=samtools, bam_in=bam_sort, str=' '.join(hbv_newid_list), bam_out=bam_sort_hbv))
    os.system("{samtools} index {bam}".format(samtools=samtools, bam=bam_sort_hbv))
    os.system("{samtools} view {bam}|cut -f 1|sort -u > {out}".format(samtools=samtools, bam=bam_sort_hbv, out=hbv_reads_id_txt))
    os.system("python %s -b %s -n %s -o %s" % (extract_py, bam_sort, hbv_reads_id_txt, bam_hbv_all) + '\n')


### get chimeric reads info
def get_names(names):
    with open(names, 'r') as infile:
        n = infile.read().splitlines()
    if '' in n:
        n.remove('')
    return n

def pysam_rnames_info(name_indexed, name):
    try:
        name_indexed.find(name)
    except KeyError:
        pass
    else:
        # get query_length
        query_length = 0
        iterator = name_indexed.find(name)
        for AlignedSegmennt in iterator:
            if not AlignedSegmennt.is_secondary and not AlignedSegmennt.is_supplementary:
                query_length = AlignedSegmennt.query_length
        iterator = name_indexed.find(name)
        out_list = []
        out_line = ""
        for AlignedSegmennt in iterator:
            mapq = AlignedSegmennt.mapping_quality
            # primary, supplementary, secondary
            if AlignedSegmennt.is_supplementary:
                supplementary = "supplementary"
            else:
                if AlignedSegmennt.is_secondary:
                    supplementary = "secondary"
                else:
                    supplementary = "primary"
            # coordinate
            query_n = AlignedSegmennt.query_name
            ref_n = AlignedSegmennt.reference_name
            reference_start = AlignedSegmennt.reference_start
            reference_end = AlignedSegmennt.reference_end
            # is_reverse
            if AlignedSegmennt.is_reverse:
                strand = 'reverse'
            else:
                strand = 'forward'
            # clip info
            cigartuples = AlignedSegmennt.cigartuples
            if cigartuples[0][0] in [4, 5]:
                left_clip_in_bam = cigartuples[0][1]
            else:
                left_clip_in_bam = 0
            if cigartuples[-1][0] in [4, 5]:
                right_clip_in_bam = cigartuples[-1][1]
            else:
                right_clip_in_bam = 0
            # get query clip info in forward strand
            if strand == "forward":
                left_clip_in_forawrd = left_clip_in_bam
                right_clip_in_forawrd = right_clip_in_bam
            else:
                left_clip_in_forawrd = right_clip_in_bam
                right_clip_in_forawrd = left_clip_in_bam
            # adjust alignment coordinate in query(all change to forward)
            query_alignment_start = left_clip_in_forawrd
            query_alignment_end = query_length - right_clip_in_forawrd
            # output
            ref_loc = "%s:%s-%s" % (ref_n, str(reference_start), str(reference_end))
            alignment_info = "|".join([supplementary, str(mapq), strand,
                                       str(query_length), str(left_clip_in_forawrd),
                                       str(right_clip_in_forawrd)])
            out_line += '\t'.join([query_n, str(query_alignment_start), str(query_alignment_end), ref_loc, alignment_info]) + '\n'
            out_list.append('\t'.join([query_n, str(query_alignment_start), str(query_alignment_end), ref_loc, alignment_info]) + '\n')
        return out_list

def get_chimeric_reads(bam_hbv_all_sort, hbv_reads_id_txt, hbv_reads_info_bed, threads):
    read_names = get_names(hbv_reads_id_txt)
    hbv_pysam = pysam.AlignmentFile(bam_hbv_all_sort, 'rb')
    name_indexed = pysam.IndexedReads(hbv_pysam)
    name_indexed.build()
    # query_loc, ref_loc, primary|strand|query_length|left_clip|right_clip
    result = []
    pools = multiprocessing.pool.ThreadPool(threads)
    for name in read_names:
        result.append(pools.apply_async(pysam_rnames_info, args=(name_indexed, name)).get())
    pools.close()
    pools.join()
    del pools
    with open(hbv_reads_info_bed, 'w') as out:
        for newline_info in result:
            for newline in newline_info:
                out.write(newline)
                out.flush()

### get breakpoint
def get_breakpoint(alingment_info, lr):
    c_x = alingment_info.split()
    ref_loc_x = c_x[3]
    chrom_x = ref_loc_x.split(':')[0]
    start_x = ref_loc_x.split(':')[1].split('-')[0]
    end_x = ref_loc_x.split(':')[1].split('-')[1]
    if lr == "left":
        if "reverse" in alingment_info:
            loc_bk = "%s:%s-%s" % (chrom_x, start_x, str(int(start_x) + 1))
        else:
            loc_bk = "%s:%s-%s" % (chrom_x, str(int(end_x) - 1), end_x)
    else:
        if "reverse" in alingment_info:
            loc_bk = "%s:%s-%s" % (chrom_x, str(int(end_x) - 1), end_x)
        else:
            loc_bk = "%s:%s-%s" % (chrom_x, start_x, str(int(start_x) + 1))
    return loc_bk

def get_breakpoint_on_hbv(alingment_info, lr):
    c_x = alingment_info.split()
    ref_loc_x = c_x[3]
    chrom_x = ref_loc_x.split(':')[0]
    start_x = ref_loc_x.split(':')[1].split('-')[0]
    end_x = ref_loc_x.split(':')[1].split('-')[1]
    if lr == "left":
        if "reverse" in alingment_info:
            loc_bk = "%s:%s-%s" % (chrom_x, str(int(end_x) - 1), end_x)
        else:
            loc_bk = "%s:%s-%s" % (chrom_x, start_x, str(int(start_x) + 1))
    else:
        if "reverse" in alingment_info:
            loc_bk = "%s:%s-%s" % (chrom_x, start_x, str(int(start_x) + 1))
        else:
            loc_bk = "%s:%s-%s" % (chrom_x, str(int(end_x) - 1), end_x)
    return loc_bk

def get_breakpoint_read(name, hbv_reads_info_bed, sampleid):
    ### just use primary alignment
    hbv_nr = os.popen(''' grep %s %s|grep -v secondary|sort -k1,1 -k2,2n|grep -n HBV|cut -f 1 -d ":" ''' % (
        name, hbv_reads_info_bed)).read().strip().split()
    # print("%s: %s" % (name, hbv_nr))
    bk_loc_list = []
    hbv_seq_list = []
    for nr in hbv_nr:
        left_nr = int(nr) - 1
        right_nr = int(nr) + 1
        # consecutive HBV alignment
        if str(left_nr) not in hbv_nr:
            left_bk = os.popen(
                ''' grep %s %s|grep -v secondary|sort -k1,1 -k2,2n|awk '{if(NR==%s) {print $0}}' ''' % (
                    name, hbv_reads_info_bed, left_nr)).read().strip()
            if left_bk:
                # print(left_bk)
                # get breakpoint on HM
                left_bk_loc = get_breakpoint(left_bk, 'left')
                # get breakpoint on HBV
                hbv_nr_left = int(left_nr) + 1
                left_bk_hbv = os.popen(
                    ''' grep %s %s|grep -v secondary|sort -k1,1 -k2,2n|awk '{if(NR==%s) {print $0}}' ''' % (
                        name, hbv_reads_info_bed, hbv_nr_left)).read().strip()
                hbv_seq_loc = left_bk_hbv.split()[3]
                if hbv_seq_loc not in hbv_seq_list:
                    hbv_seq_list.append(hbv_seq_loc)
                left_bk_loc_hbv = get_breakpoint_on_hbv(left_bk_hbv, 'left')
                bk_loc_list.append(left_bk_loc + '|' + left_bk_loc_hbv)
            else:
                # print("left empty")
                bk_loc_list.append("N")
        if str(right_nr) not in hbv_nr:
            right_bk = os.popen(
                ''' grep %s %s|grep -v secondary|sort -k1,1 -k2,2n|awk '{if(NR==%s) {print $0}}' ''' % (
                    name, hbv_reads_info_bed, right_nr)).read().strip()
            if right_bk:
                # print(right_bk)
                # get breakpoint on HM
                right_bk_loc = get_breakpoint(right_bk, 'right')
                # get breakpoint on HBV
                hbv_nr_right = int(right_nr) - 1
                right_bk_hbv = os.popen(
                    ''' grep %s %s|grep -v secondary|sort -k1,1 -k2,2n|awk '{if(NR==%s) {print $0}}' ''' % (
                        name, hbv_reads_info_bed, hbv_nr_right)).read().strip()
                hbv_seq_loc = right_bk_hbv.split()[3]
                if hbv_seq_loc not in hbv_seq_list:
                    hbv_seq_list.append(hbv_seq_loc)
                right_bk_loc_hbv = get_breakpoint_on_hbv(right_bk_hbv, 'right')
                bk_loc_list.append(right_bk_loc + '|' + right_bk_loc_hbv)
            else:
                # print("right empty")
                bk_loc_list.append("N")
    # out txt:
    out_list1 = [name, ','.join(bk_loc_list), ','.join(hbv_seq_list)]
    # out bed
    out_list2 = []
    for bk_loc in bk_loc_list:
        if bk_loc != 'N':
            bk_loc_info = bk_loc.split('|')
            chrom = bk_loc_info[0].split(':')[0]
            start = bk_loc_info[0].split(':')[1].split('-')[0]
            end = bk_loc_info[0].split(':')[1].split('-')[1]
            out_list2.append('|'.join([chrom, start, end, bk_loc_info[1], sampleid, name]))
    return [out_list1, out_list2]


def get_breakpoint_out(hbv_reads_info_bed, hbv_reads_id_txt, hbv_bk_hm_txt, hbv_bk_hm_bed, sampleid):
    read_names = get_names(hbv_reads_id_txt)
    result = []
    pools = multiprocessing.Pool(10)
    for name in read_names:
        result.append(pools.apply_async(get_breakpoint_read, args=(name, hbv_reads_info_bed, sampleid)).get())
    pools.close()
    pools.join()
    del pools
    with open(hbv_bk_hm_txt, 'w') as out:
        with open(hbv_bk_hm_bed, 'w') as out_bed:
            for bk_info in result:
                out.write('\t'.join(bk_info[0]) + '\n')
                out.flush()
                if len(bk_info[1]) != 0:
                    for bk_info2 in bk_info[1]:
                        out_bed.write('\t'.join(bk_info2.split('|')) + '\n')
                        out_bed.flush()

def run(options):
    hbv_id_txt = options.names
    bam_sort = options.bam
    sampleid = options.sampleid
    out_dir = options.out_dir
    threads = options.threads
    # 1 extract HBV bam
    print(f"{datetime.datetime.now()} *** extract HBV bam start.")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    bam_sort_hbv = os.path.join(out_dir, "%s_minimap2_sorted_hbv.bam" % sampleid)
    hbv_reads_id_txt = os.path.join(out_dir, "%s_minimap2_hbv_readsid.txt" % sampleid)
    bam_hbv_all = os.path.join(out_dir, "%s_minimap2_hbv_all.bam" % sampleid)
    extract_hbv_readsid(hbv_id_txt, bam_sort, bam_sort_hbv, hbv_reads_id_txt, bam_hbv_all)
    print(f"{datetime.datetime.now()} *** extract HBV bam finished.")
    # 2
    print(f"{datetime.datetime.now()} *** breakponts start: chimeric")
    bam_hbv_all_sort = os.path.join(out_dir, "%s_minimap2_hbv_all_sorted.bam" % sampleid)
    hbv_reads_info_bed = os.path.join(out_dir, "%s_minimap2_hbv_chimeric_reads.bed" % sampleid)
    get_chimeric_reads(bam_hbv_all_sort, hbv_reads_id_txt, hbv_reads_info_bed, threads)
    # 3
    sample_bk_dir = os.path.join(out_dir, "breakpoint")
    if not os.path.exists(sample_bk_dir):
        os.makedirs(sample_bk_dir)
    hbv_bk_hm_txt = os.path.join(sample_bk_dir, "%s_chimeric_breakpoint_hm_hbv.txt" % sampleid)
    hbv_bk_hm_bed = os.path.join(sample_bk_dir, "%s_chimeric_breakpoint_hm_hbv.bed" % sampleid)
    get_breakpoint_out(hbv_reads_info_bed, hbv_reads_id_txt, hbv_bk_hm_txt, hbv_bk_hm_bed, sampleid)
    print(f"{datetime.datetime.now()} *** breakponts finished: chimeric")

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Identify HBV integration events by chrmeric reads')
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-n', '--names', help='virus ids per line', required=True)
    parser.add_argument('-s', '--sampleid', help='output file prefix', required=True)
    parser.add_argument('-o', '--out_dir', help='output directory', required=True)
    parser.add_argument('-t', '--threads', help='number of threads (default: 30)', type=int, default=30)
    options = parser.parse_args()
    run(options)
