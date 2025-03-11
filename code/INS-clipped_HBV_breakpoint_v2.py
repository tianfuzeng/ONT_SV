# __author__ = 'tianfuzneg'
# !/usr/bin/python
# -*- coding:utf-8 -*-

########################################################################################
# 20241129
# use minimap2 aligned to main-HBV+hg38
########################################################################################
import os
import datetime
import pysam
import multiprocessing.pool
import pandas as pd
import subprocess
import pandas as pd

samtools = "/data/fs01/biosoft/samtools-1.9/samtools"
minimap2 = "/data/fs01/wangzf/software/minimap2-2.24_x64-linux/minimap2"
bedtools = "/data/fs01/wangzf/software/bedtools-2.30.0"

clip_py = '/data/fs09/wangzf/nanopore/ztf/HCC_xm/program/extract_clipped_seq_from_bam_by_ids_v4.py'
ins_py = "/data/fs09/wangzf/nanopore/ztf/HCC_xm/program/extract_INS_seq_v2.py"

########################################################################################
# INS-clipped speed up: minimap2 to main-HBV + human
########################################################################################
# extract ins and clipped seq
# clip-INS seq re-alignment
def get_rnames(bam, txt):
    os.system("%s view %s | awk '{print $1}' | sort | uniq > %s" % (samtools, bam, txt))

def extract_ins_clipped_seq(rnames_bam, rnames_txt, clipped_fa, ins_fa):
    os.system("python %s -b %s -n %s -o %s" % (clip_py, rnames_bam, rnames_txt, clipped_fa))
    os.system("python %s -b %s -n %s -o %s" % (ins_py, rnames_bam, rnames_txt, ins_fa))

# re-alignment
def m2_sr(ref, fa, sam_x):
    bam_x = sam_x.replace('.sam', '.bam')
    bam_sort_x = sam_x.replace('.sam', '_sorted.bam')
    os.system(f"{minimap2} -ax sr --MD -t 10 {ref} {fa} -o {sam_x}")  # -N 10 -p 0.1
    os.system(f"{samtools} view -b -S -F 4 %s > %s" % (sam_x, bam_x))
    os.system(f"{samtools} sort -@ 10 %s > %s" % (bam_x, bam_sort_x))
    os.system(f"{samtools} index -@ 10 %s" % bam_sort_x)
    os.system("rm %s" % sam_x)
    os.system("rm %s" % bam_x)

def extract_hbv_bam(hbv_id_txt, bam_sort, bam_sort_hbv):
    hbv_newid_list = os.popen("cat %s" % hbv_id_txt).read().strip().split()
    os.system("{samtools} view -@ 5 -hb {bam_in} {str} > {bam_out}".format(
        samtools=samtools, bam_in=bam_sort, str=' '.join(hbv_newid_list), bam_out=bam_sort_hbv))
    os.system("{samtools} index {bam}".format(samtools=samtools, bam=bam_sort_hbv))

########################################################################################
# extract breakpoints
########################################################################################
### get re-alignment reads info: queryid, query_start, query_end, ref_loc(hbv_loc), infos
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
            mapq = str(AlignedSegmennt.mapping_quality)
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
            alignment_info = "|".join([supplementary, mapq, strand,
                                       str(query_length), str(left_clip_in_forawrd),
                                       str(right_clip_in_forawrd)])
            out_line += '\t'.join([query_n, str(query_alignment_start), str(query_alignment_end), ref_loc, alignment_info]) + '\n'
            out_list.append('\t'.join([query_n, str(query_alignment_start), str(query_alignment_end), ref_loc, alignment_info]) + '\n')
        return out_list

def get_chimeric_reads(bam_hbv_all_sort, hbv_reads_id_txt, hbv_reads_info_bed, threads):  # from chimeric pip
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
            if newline_info is not None:
                for newline in newline_info:
                    out.write(newline)
                    out.flush()

### breakpoints
def get_loc_info(loc_info_x):
    # chr6_70300829-chr2_74273672
    c_x = loc_info_x.split('_')
    chr_x = c_x[0]
    start_x = c_x[1]
    end_x = str(int(c_x[1]) + 1)
    return [chr_x, start_x, end_x]

def get_bp(hbv_reads_info_bed, hm_loc_bed, hm_hbv_txt, type1, sampleid, rnames_txt):
    with open(hm_loc_bed, 'w') as out, open(hm_hbv_txt, 'w') as out1:
        if os.path.getsize(hbv_reads_info_bed):
            read_names = get_names(rnames_txt)
            column_names = ['queryid', 'query_start', 'query_end', 'ref_loc', 'infos']
            df = pd.read_csv(hbv_reads_info_bed, sep='\t', header=None, names=column_names)
            for readid in read_names:
                filtered_df = df[df['queryid'] == readid]
                filtered_df = filtered_df[filtered_df['ref_loc'].str.contains('HBV', na=False)]
                remaining_rows = len(filtered_df)
                if remaining_rows == 0:
                    print(f"{readid} not on HBV")
                else:
                    # hm breakpoint
                    c1 = readid.split(':')
                    read_name = c1[0]
                    c2 = c1[1].split('_')
                    if type1 == "INS":
                        hm_loc_list = [c2[3:]]
                    else:
                        hm_loc_list = []
                        c3 = '_'.join(c2[3:]).split('-')
                        for loc_info in c3:
                            if 'HBV' not in loc_info and loc_info != 'N':
                                hm_loc_list.append(get_loc_info(loc_info))
                            elif loc_info == 'N':
                                hm_loc_list.append('N')
                    # hbv loc
                    for index, row in filtered_df.iterrows():  # if 1 query aligned to multiple HBV, output separately
                        ref_loc_value = row['ref_loc']
                        ref_loc_info = ref_loc_value.split(':')
                        hbv_name = ref_loc_info[0]
                        ref_loc_loc = ref_loc_info[1].split('-')
                        hbv_loc = hbv_name + ':' + '_'.join(ref_loc_loc)
                        hbv_loc1 = hbv_name + ':' + '_'.join([ref_loc_loc[0], str(int(ref_loc_loc[0]) + 1)])
                        hbv_loc2 = hbv_name + ':' + '_'.join([ref_loc_loc[1], str(int(ref_loc_loc[1]) + 1)])
                        if 'reverse' in filtered_df.iloc[0]['infos']:
                            hbv_loc_list = [hbv_loc2, hbv_loc1]
                        else:
                            hbv_loc_list = [hbv_loc1, hbv_loc2]
                        # hm-hbv
                        out_list2_info = []
                        for hm_loc in hm_loc_list:
                            if hm_loc != 'N':
                                hbv_loc_pair = hbv_loc_list[hm_loc_list.index(hm_loc)]  # breakpoints paired info
                                out_list = hm_loc + [hbv_loc_pair, sampleid, read_name]
                                print(out_list)
                                out.write('\t'.join(out_list) + '\n')
                                out.flush()
                                out_list2_info.append(
                                    '%s|%s' % ("%s:%s-%s" % (hm_loc[0], hm_loc[1], hm_loc[2]), hbv_loc_pair))
                            else:
                                out_list2_info.append('N')
                        # piared txt, 1 event per line
                        out_list2 = [read_name, ','.join(out_list2_info), hbv_loc]
                        out1.write('\t'.join(out_list2) + '\n')

########################################################################################
# run configuration
########################################################################################
def run(options):
    sampleid = options.sampleid
    out_dir = options.out_dir
    bam_sort = options.bam  # '%s_minimap2_sorted.bam'
    threads = options.threads
    hg38_main_hbv_fa = options.main_fa
    main_hbv_txt = options.main_hbv
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # 1. get ins and clipped fasta
    print(f"{datetime.datetime.now()} *** extract_ins_clipped_seq start.")
    rnames_txt = f"{out_dir}/{sampleid}_readsid.txt"
    get_rnames(bam_sort, rnames_txt)
    clipped_fa = f"{out_dir}/{sampleid}_clipped.fasta"
    ins_fa = f"{out_dir}/{sampleid}_INS.fasta"
    extract_ins_clipped_seq(bam_sort, rnames_txt, clipped_fa, ins_fa)
    print(f"{datetime.datetime.now()} *** extract_ins_clipped_seq finished.")
    for type1 in ['INS', 'clipped']:
        # 2. re-alignment
        print(f"{datetime.datetime.now()} *** minimap2 start: {type1}")
        ic_fa = f"{out_dir}/{sampleid}_{type1}.fasta"
        ic_sam = f"{out_dir}/{sampleid}_{type1}_minimap2.sam"
        ic_bam = ic_sam.replace('.sam', '_sorted.bam')
        ic_bam_hbv = ic_sam.replace('.sam', '_hbv.bam')  # f"{out_dir}/{sampleid}_{type1}_minimap2_hbv.bam"
        m2_sr(hg38_main_hbv_fa, ic_fa, ic_sam)
        extract_hbv_bam(main_hbv_txt, ic_bam, ic_bam_hbv)
        print(f"{datetime.datetime.now()} *** minimap2 finished: {type1}")
        # 3. brekpoints
        print(f"{datetime.datetime.now()} *** brekpoints start: {type1}")
        # get hbv bam and read names
        hbv_rnames_txt = f"{out_dir}/{sampleid}_{type1}_minimap2_hbv_readsid.txt"
        get_rnames(ic_bam_hbv, hbv_rnames_txt)
        # get re-alignment reads info
        hbv_reads_info_bed = f"{out_dir}/{sampleid}_minimap2_hbv_{type1}_reads.bed"
        get_chimeric_reads(ic_bam, hbv_rnames_txt, hbv_reads_info_bed, threads)
        # breakpoints
        sample_bk_dir = os.path.join(out_dir, "breakpoint")
        if not os.path.exists(sample_bk_dir):
            os.makedirs(sample_bk_dir)
        hm_loc_bed_clip = f"{sample_bk_dir}/{sampleid}_{type1}_breakpoint_hm_hbv.bed"
        hm_hbv_txt_clip = f"{sample_bk_dir}/{sampleid}_{type1}_breakpoint_hm_hbv.txt"
        get_bp(hbv_reads_info_bed, hm_loc_bed_clip, hm_hbv_txt_clip, type1, sampleid, hbv_rnames_txt)
        print(f"{datetime.datetime.now()} *** brekpoints finished: {type1}")

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Identify HBV integration events by INS-clipped reads')
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-s', '--sampleid', help='output file prefix', required=True)
    parser.add_argument('-o', '--out_dir', help='output directory', required=True)
    parser.add_argument('-fa', '--main_fa', help='hg38 and main hbv (fa)', required=True)
    parser.add_argument('-hbv', '--main_hbv', help='main hbv genotype (txt)', required=True)
    parser.add_argument('-t', '--threads', help='number of threads (default: 30)', type=int, default=30)
    options = parser.parse_args()
    run(options)




