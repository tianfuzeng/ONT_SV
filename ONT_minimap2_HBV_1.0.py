# __author__ = 'tianfuzneg'
# !/usr/bin/python
# -*- coding:utf-8 -*-

########################################################################################
# 20221107
# new 21 HBV reference
########################################################################################
import os
import pysam

ref_dir = '/data/fs09/wangzf/nanopore/ztf/HCC/ref'
hg38_fa = os.path.join(ref_dir, 'hg38_mainChr.fa')
hg38_hbv_fa = os.path.join(ref_dir, 'HBV_ztf', 'hg38_mainChr_and_21_HBV_genome.fa')
# gene5kb_cn_bed = "/data/fs08/wangzf/hg38_ztf/gencode_v33_annotation_gtf_gene5kb_cosmic-ncg6.bed"

# minimap2 = "/data/fs01/wangzf/software/minimap2/minimap2"
minimap2 = "/data/fs01/wangzf/software/minimap2-2.24_x64-linux/minimap2"
samtools = "/data/fs01/biosoft/samtools-1.9/samtools"
sniffles2 = "/data/fs01/wangzf/software/anaconda3/envs/nanoplot/bin/sniffles"
tr_bed = "/data/fs09/wangzf/nanopore/ztf/HCC/ref/sniffles/human_GRCh38_no_alt_analysis_set.trf.bed"

program_dir = '/data/fs09/wangzf/nanopore/ztf/HCC/ONT/program/'
work_dir = '/data/fs08/wangzf/nanopore/ztf/HCC/ONT/HBV_minimap2'
ns_dir = '/data/fs09/wangzf/nanopore/ztf/HCC/ONT/ngmlr_sniffles'

sample_list = ['HCC8_WBC', 'HCC8_N1', 'HCC8_N3', 'HCC8_T1', 'HCC8_T2', 'HCC8_T3', 'HCC8_T4', 'HCC8_T5',
               'HCC9_WBC', 'HCC9_N1', 'HCC9_N3', 'HCC9_T1', 'HCC9_T2', 'HCC9_T3', 'HCC9_T4', 'HCC9_T5',
               'HCC10_WBC', 'HCC10_N1', 'HCC10_N3', 'HCC10_N5', 'HCC10_T1', 'HCC10_T2', 'HCC10_T3', 'HCC10_T4', 'HCC10_T5', 'HCC10_T6', 'HCC10_T7',
               'HCC12_WBC', 'HCC12_N1', 'HCC12_N3', 'HCC12_T1', 'HCC12_T2', 'HCC12_T3', 'HCC12_T4', 'HCC12_T5',
               'HCC13_WBC', 'HCC13_N1', 'HCC13_T1', 'HCC13_T2', 'HCC13_T3', 'HCC13_T4', 'HCC13_T5']
sample_n50_list = ['HCC8_T4', 'HCC9_N3', 'HCC9_T2', 'HCC10_N5', 'HCC12_T3', 'HCC13_T2']
sample_list_3rd = ['HCC13_T2', 'HCC13_T5', 'HCC13_N1', 'HCC13_T4', 'HCC13_T1', 'HCC13_WBC', 'HCC13_T3']
sample_list_16 = ['HCC10_T5', 'HCC10_T3', 'HCC10_T6', 'HCC10_T7', 'HCC10_N3', 'HCC12_T1', 'HCC12_T2', 'HCC12_T4',
                  'HCC10_T4', 'HCC12_T5', 'HCC12_N1', 'HCC12_N3', 'HCC10_N1', 'HCC10_N5', 'HCC12_T3', 'HCC10_T2']

########################################################################################
# minimap2
########################################################################################
script = os.path.join(program_dir, 'ONT_minimap2_HBV_1.0.sh')
IDs = [31, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50][10:]
i = 0
for sampleid in sample_list:
    # fastq configure
    if sampleid in sample_n50_list:
        fq = os.path.join('/data/fs09/wangzf/nanopore/ztf/HCC/ONT/data/fastq_5th', sampleid, "%s_clean.fq.gz" % sampleid)
    elif sampleid in ['HCC10_WBC']:
        fq = os.path.join("/data/fs09/wangzf/nanopore/ztf/HCC/ONT/data/20220215", sampleid, "%s_clean.fq.gz" % sampleid)
    elif sampleid in sample_list_3rd:
        fq = os.path.join('/data/fs09/wangzf/nanopore/ztf/HCC/ONT/data/fastq_3rd', sampleid, "%s_clean.fq.gz" % sampleid)
    elif sampleid in sample_list_16:
        fq = os.path.join("/data/fs09/wangzf/nanopore/ztf/HCC/ONT/data/20210929", sampleid, "%s_clean.fq.gz" % sampleid)
    else:
        fq = os.path.join("/data/fs09/wangzf/nanopore/ztf/HCC/ONT/fastq_pass", sampleid, "%s_clean.fq.gz" % sampleid)
    # run
    output_dir = os.path.join(work_dir, 'output')
    sample_dir = os.path.join(output_dir, sampleid)
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)
    stdout = os.path.join(sample_dir, '%s_pipeline.o' % sampleid)
    stderr = os.path.join(sample_dir, '%s_pipeline.e' % sampleid)
    for std in [stdout, stderr]:
        if os.path.exists(std):
            os.system("rm %s" % std)
    if i <= 5:
        os.system(
            'qsub -l hostname=PMC-{server} -S /bin/bash -o {out} -e {err} -N {name} -cwd {script} {a} {b} {c}'.format(
                server=IDs[i], out=stdout, err=stderr, name="m2_%s" % sampleid[4:], script=script, a=sampleid, b=fq, c=output_dir))
        i = i + 1
    else:
        i = 0
        os.system(
            'qsub -l hostname=PMC-{server} -S /bin/bash -o {out} -e {err} -N {name} -cwd {script} {a} {b} {c}'.format(
                server=IDs[i], out=stdout, err=stderr, name="m2_%s" % sampleid[4:], script=script, a=sampleid, b=fq,
                c=output_dir))
        i = i + 1

# check
for sampleid in sample_list:
    output_dir = os.path.join(work_dir, 'output')
    sample_dir = os.path.join(output_dir, sampleid)
    stderr = os.path.join(sample_dir, '%s_pipeline.e' % sampleid)
    if 'HCC13' in sampleid:
        print(sampleid)
        print(os.popen("tail %s" % stderr).read())

########################################################################################
# extract reads on HBV genome
########################################################################################
# grep ">" HBV_genome_21_newid.fasta|cut -d ">" -f 2 > HBV_genome_21_newid.txt
hbv_newid_txt = os.path.join(ref_dir, 'HBV_ztf', 'HBV_genome_21_newid.txt')
hbv_newid_list = os.popen("cat %s" % hbv_newid_txt).read().strip().split()

for sampleid in sample_list:
    bam_sort = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_sorted.bam" % sampleid)
    bam_sort_hbv = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_sorted_hbv.bam" % sampleid)
    os.system("{samtools} view -@ 5 -hb {bam_in} {str} > {bam_out}".format(
        samtools=samtools, bam_in=bam_sort, str=' '.join(hbv_newid_list), bam_out=bam_sort_hbv))
    os.system("{samtools} index {bam}".format(samtools=samtools, bam=bam_sort_hbv))

########################################################################################
# if reads contain HBV seq, extract all reads alignments results
# https://timoast.github.io/blog/2015-10-12-extractreads/
########################################################################################
# get reads id
for sampleid in sample_list:
    bam_sort = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_sorted.bam" % sampleid)
    bam_sort_hbv = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_sorted_hbv.bam" % sampleid)
    hbv_reads_id_txt = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_readsid.txt" % sampleid)
    os.system("{samtools} view {bam}|cut -f 1|sort -u > {out}".format(samtools=samtools, bam=bam_sort_hbv, out=hbv_reads_id_txt))

# extract bam
extract_py = os.path.join(program_dir, 'extract_reads_from_bam_by_ids.py')
for sampleid in sample_list:
    bam_sort = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_sorted.bam" % sampleid)
    hbv_reads_id_txt = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_readsid.txt" % sampleid)
    bam_sort_hbv_all = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_sorted_hbv_all.bam" % sampleid)
    script_pip = os.path.join(work_dir, 'output', sampleid, "%s_extract_reads.sh" % sampleid)
    with open(script_pip, 'w') as out:
        out.write("#! /bin/bash" + '\n')
        out.write('''echo "$(date) 1. Start : %s" ''' % sampleid + '\n')
        out.write("python %s -b %s -n %s -o %s" % (extract_py, bam_sort, hbv_reads_id_txt, bam_sort_hbv_all) + '\n')
        out.write('''echo "$(date) 1. Finish : %s" ''' % sampleid + '\n')

i = 0
IDs = [32, 35, 49, 50]
for sampleid in sample_list:
    script_pip = os.path.join(work_dir, 'output', sampleid, "%s_extract_reads.sh" % sampleid)
    stdout = script_pip.replace(".sh", ".o")
    stderr = script_pip.replace(".sh", ".e")
    for std in [stdout, stderr]:
        if os.path.exists(std):
            os.system("rm %s" % std)
    if i <= 3:
        os.system(
            'qsub -l hostname=PMC-{server} -S /bin/bash -o {out} -e {err} -N {name} -cwd {script}'.format(
                server=IDs[i], out=stdout, err=stderr, name="e_%s" % sampleid.replace('HCC', ''), script=script_pip))
        i = i + 1
    else:
        i = 0
        os.system(
            'qsub -l hostname=PMC-{server} -S /bin/bash -o {out} -e {err} -N {name} -cwd {script}'.format(
                server=IDs[i], out=stdout, err=stderr, name="e_%s" % sampleid.replace('HCC', ''), script=script_pip))
        i = i + 1

for sampleid in sample_list:
    if "HCC8" not in sampleid:
        bam_sort_hbv_all = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_sorted_hbv_all.bam" % sampleid)
        bam_hbv_all_sort = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_all_sorted.bam" % sampleid)
        os.system("{samtools} sort {bam_in} > {bam_out}".format(
            samtools=samtools, bam_in=bam_sort_hbv_all, bam_out=bam_hbv_all_sort))
        os.system("{samtools} index {bam}".format(samtools=samtools, bam=bam_hbv_all_sort))
        os.system("rm %s" % bam_sort_hbv_all)

########################################################################################
# get breakpoints using sniffles
########################################################################################
for sampleid in ['HCC8_T2']:
    bam_hbv_all_sort = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_all_sorted.bam" % sampleid)
    bam_hbv_vcf = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_all_sniffles.vcf" % sampleid)
    os.system("{sniffles} -i {bam_sort} -v {vcf} --tandem-repeats {tr} -t 4 --minsupport 1 --minsvlen 50 "
              "--mapq 0 --min-alignment-length 100 --output-rnames --allow-overwrite "
              "--long-ins-length 100000 ".format(sniffles=sniffles2, bam_sort=bam_hbv_all_sort, vcf=bam_hbv_vcf, tr=tr_bed))

########################################################################################
# get chimeric reads info
########################################################################################
def get_names(names):
    with open(names, 'r') as infile:
        n = infile.read().splitlines()
    if '' in n:
        n.remove('')
    return n


for sampleid in sample_list:
    bam_hbv_all_sort = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_all_sorted.bam" % sampleid)
    hbv_reads_id_txt = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_readsid.txt" % sampleid)
    read_names = get_names(hbv_reads_id_txt)
    hbv_pysam = pysam.AlignmentFile(bam_hbv_all_sort, 'rb')
    name_indexed = pysam.IndexedReads(hbv_pysam)
    name_indexed.build()
    hbv_reads_info_txt = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_chimeric_reads.bed" % sampleid)
    # query_loc, ref_loc, primary|strand|query_length|left_clip|right_clip
    with open(hbv_reads_info_txt, 'w') as out:
        for name in read_names:
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
                for AlignedSegmennt in iterator:
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
                    query_loc = "\t".join([query_n, str(query_alignment_start), str(query_alignment_end)])
                    ref_loc = "%s:%s-%s" % (ref_n, str(reference_start), str(reference_end))
                    alignment_info = "|".join([supplementary, strand,
                                               str(query_length), str(left_clip_in_forawrd),
                                               str(right_clip_in_forawrd)])
                    out_list = [query_loc, ref_loc, alignment_info]
                    out.write('\t'.join(out_list) + '\n')
                    out.flush()

########################################################################################
# get breakpoint
########################################################################################
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


for sampleid in sample_list:
    hbv_reads_info_txt = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_chimeric_reads.bed" % sampleid)
    hbv_reads_id_txt = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_readsid.txt" % sampleid)
    read_names = get_names(hbv_reads_id_txt)
    sample_bk_dir = os.path.join(work_dir, 'output', sampleid, "breakpoint")
    if not os.path.exists(sample_bk_dir):
        os.makedirs(sample_bk_dir)
    hbv_bk_hm_txt = os.path.join(sample_bk_dir, "%s_minimap2_hbv_breakpoint_hm.txt" % sampleid)
    hbv_bk_hm_bed = os.path.join(sample_bk_dir, "%s_minimap2_hbv_breakpoint_hm.bed" % sampleid)
    out_bed = open(hbv_bk_hm_bed, 'w')
    with open(hbv_bk_hm_txt, 'w') as out:
        for name in read_names:
            ### just use primary alignment
            hbv_nr = os.popen(''' grep %s %s|grep -v secondary|sort -k1,1 -k2,2n|grep -n HBV|cut -f 1 -d ":" ''' % (
                name, hbv_reads_info_txt)).read().strip().split()
            print("%s: %s" % (name, hbv_nr))
            bk_loc_list = []
            for nr in hbv_nr:
                left_nr = int(nr) - 1
                right_nr = int(nr) + 1
                # consecutive HBV alignment
                if str(left_nr) not in hbv_nr:
                    left_bk = os.popen(
                        ''' grep %s %s|grep -v secondary|sort -k1,1 -k2,2n|awk '{if(NR==%s) {print $0}}' ''' % (
                            name, hbv_reads_info_txt, left_nr)).read().strip()
                    if left_bk:
                        print(left_bk)
                        left_bk_loc = get_breakpoint(left_bk, 'left')
                        if left_bk_loc != 'HBV':
                            bk_loc_list.append(left_bk_loc)
                    else:
                        print("left empty")
                        bk_loc_list.append("N")
                else:
                    print("consecutive HBV")
                if str(right_nr) not in hbv_nr:
                    right_bk = os.popen(
                        ''' grep %s %s|grep -v secondary|sort -k1,1 -k2,2n|awk '{if(NR==%s) {print $0}}' ''' % (
                            name, hbv_reads_info_txt, right_nr)).read().strip()
                    if right_bk:
                        print(right_bk)
                        right_bk_loc = get_breakpoint(right_bk, 'right')
                        if right_bk_loc != 'HBV':
                            bk_loc_list.append(right_bk_loc)
                    else:
                        print("right empty")
                        bk_loc_list.append("N")
                else:
                    print("consecutive HBV")
            out_list = [name, ','.join(bk_loc_list)]
            out.write('\t'.join(out_list) + '\n')
            out.flush()
            # out bed
            for bk_loc in bk_loc_list:
                if bk_loc != 'N':
                    chrom = bk_loc.split(':')[0]
                    start = bk_loc.split(':')[1].split('-')[0]
                    end = bk_loc.split(':')[1].split('-')[1]
                    out_list2 = [chrom, start, end, sampleid, name]
                    out_bed.write('\t'.join(out_list2) + '\n')
    out_bed.close()


########################################################################################
# get breakpoint 2.0: add HBV breakpoint; add secondary alinment
########################################################################################
def extract_breakoint(name_x, hbv_reads_info_txt_x, type1):
    hbv_nr = os.popen(''' grep %s %s|grep %s|sort -k1,1 -k2,2n|grep -n HBV|cut -f 1 -d ":" ''' % (
        name_x, hbv_reads_info_txt_x, type1)).read().strip().split()
    print("%s: %s" % (name_x, hbv_nr))
    bk_loc_list = []
    hbv_seq_list = []
    for nr in hbv_nr:
        left_nr = int(nr) - 1
        right_nr = int(nr) + 1
        # consecutive HBV alignment
        if str(left_nr) not in hbv_nr:
            left_bk = os.popen(
                ''' grep %s %s|grep %s|sort -k1,1 -k2,2n|awk '{if(NR==%s) {print $0}}' ''' % (
                    name_x, hbv_reads_info_txt_x, type1, left_nr)).read().strip()
            if left_bk:
                print(left_bk)
                # get breakpoint on HM
                left_bk_loc = get_breakpoint(left_bk, 'left')
                # get breakpoint on HBV
                hbv_nr_left = int(left_nr) + 1
                left_bk_hbv = os.popen(
                    ''' grep %s %s|grep %s|sort -k1,1 -k2,2n|awk '{if(NR==%s) {print $0}}' ''' % (
                        name_x, hbv_reads_info_txt_x, type1, hbv_nr_left)).read().strip()
                hbv_seq_loc = left_bk_hbv.split()[3]
                if hbv_seq_loc not in hbv_seq_list:
                    hbv_seq_list.append(hbv_seq_loc)
                left_bk_loc_hbv = get_breakpoint_on_hbv(left_bk_hbv, 'left')
                bk_loc_list.append(left_bk_loc + '|' + left_bk_loc_hbv)
            else:
                print("left empty")
                bk_loc_list.append("N")
        else:
            print("consecutive HBV")
        if str(right_nr) not in hbv_nr:
            right_bk = os.popen(
                ''' grep %s %s|grep %s|sort -k1,1 -k2,2n|awk '{if(NR==%s) {print $0}}' ''' % (
                    name_x, hbv_reads_info_txt_x, type1, right_nr)).read().strip()
            if right_bk:
                print(right_bk)
                # get breakpoint on HM
                right_bk_loc = get_breakpoint(right_bk, 'right')
                # get breakpoint on HBV
                hbv_nr_right = int(right_nr) - 1
                right_bk_hbv = os.popen(
                    ''' grep %s %s|grep %s|sort -k1,1 -k2,2n|awk '{if(NR==%s) {print $0}}' ''' % (
                        name_x, hbv_reads_info_txt_x, type1, hbv_nr_right)).read().strip()
                hbv_seq_loc = right_bk_hbv.split()[3]
                if hbv_seq_loc not in hbv_seq_list:
                    hbv_seq_list.append(hbv_seq_loc)
                right_bk_loc_hbv = get_breakpoint_on_hbv(right_bk_hbv, 'right')
                bk_loc_list.append(right_bk_loc + '|' + right_bk_loc_hbv)
            else:
                print("right empty")
                bk_loc_list.append("N")
        else:
            print("consecutive HBV")
    return [bk_loc_list, hbv_seq_list]


for sampleid in sample_list:
    hbv_reads_info_txt = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_chimeric_reads.bed" % sampleid)
    hbv_reads_id_txt = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_readsid.txt" % sampleid)
    read_names = get_names(hbv_reads_id_txt)
    sample_bk_dir = os.path.join(work_dir, 'output', sampleid, "breakpoint2.0")
    if not os.path.exists(sample_bk_dir):
        os.makedirs(sample_bk_dir)
    # colnames: read_name HM_breakpoint_left|HBV_breakpoint_left,HM_breakpoint_right|HBV_breakpoint_right HBV_seq_loc
    hbv_bk_hm_txt = os.path.join(sample_bk_dir, "%s_minimap2_breakpoint_hm_hbv_secondary.txt" % sampleid)
    # colnames: chr start end(human)  HBV_breakpoint sampleid read_name
    hbv_bk_hm_bed = os.path.join(sample_bk_dir, "%s_minimap2_breakpoint_hm_hbv_secondary.bed" % sampleid)
    out_bed = open(hbv_bk_hm_bed, 'w')
    with open(hbv_bk_hm_txt, 'w') as out:
        for name in read_names:
            ### 1:use primary alignment, get NR of HBV
            hbv_nr = os.popen(''' grep %s %s|grep -v secondary|sort -k1,1 -k2,2n|grep -n HBV|cut -f 1 -d ":" ''' % (
                name, hbv_reads_info_txt)).read().strip().split()
            print("%s: %s" % (name, hbv_nr))
            bk_loc_list = extract_breakoint(name, hbv_reads_info_txt, '-v secondary')[0]
            hbv_seq_list = extract_breakoint(name, hbv_reads_info_txt, '-v secondary')[1]
            ### 2： if primary without breakpoint, then use secondary
            if len(bk_loc_list) == 0:
                bk_loc_list = extract_breakoint(name, hbv_reads_info_txt, 'secondary')[0]
                hbv_seq_list = extract_breakoint(name, hbv_reads_info_txt, '-v secondary')[1]
            out_list = [name, ','.join(bk_loc_list), ','.join(hbv_seq_list)]
            out.write('\t'.join(out_list) + '\n')
            out.flush()
            # out bed
            for bk_loc in bk_loc_list:
                if bk_loc != 'N':
                    bk_loc_info = bk_loc.split('|')
                    chrom = bk_loc_info[0].split(':')[0]
                    start = bk_loc_info[0].split(':')[1].split('-')[0]
                    end = bk_loc_info[0].split(':')[1].split('-')[1]
                    out_list2 = [chrom, start, end, bk_loc_info[1], sampleid, name]
                    out_bed.write('\t'.join(out_list2) + '\n')
    out_bed.close()

########################################################################################
# 20231101 get breakpoint 3.0: add HBV breakpoint; delete secondary alinment
########################################################################################
for sampleid in sample_list:
    hbv_reads_info_txt = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_chimeric_reads.bed" % sampleid)
    hbv_reads_id_txt = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_hbv_readsid.txt" % sampleid)
    read_names = get_names(hbv_reads_id_txt)
    sample_bk_dir = os.path.join(work_dir, 'output', sampleid, "breakpoint3.0")
    if not os.path.exists(sample_bk_dir):
        os.makedirs(sample_bk_dir)
    # colnames: read_name HM_breakpoint_left|HBV_breakpoint_left,HM_breakpoint_right|HBV_breakpoint_right HBV_seq_loc
    hbv_bk_hm_txt = os.path.join(sample_bk_dir, "%s_minimap2_breakpoint_hm_hbv.txt" % sampleid)
    # colnames: chr start end(human)  HBV_breakpoint sampleid read_name
    hbv_bk_hm_bed = os.path.join(sample_bk_dir, "%s_minimap2_breakpoint_hm_hbv.bed" % sampleid)
    out_bed = open(hbv_bk_hm_bed, 'w')
    with open(hbv_bk_hm_txt, 'w') as out:
        for name in read_names:
            ### 1:use primary alignment, get NR of HBV
            hbv_nr = os.popen(''' grep %s %s|grep -v secondary|sort -k1,1 -k2,2n|grep -n HBV|cut -f 1 -d ":" ''' % (
                name, hbv_reads_info_txt)).read().strip().split()
            print("%s: %s" % (name, hbv_nr))
            bk_loc_list = extract_breakoint(name, hbv_reads_info_txt, '-v secondary')[0]
            hbv_seq_list = extract_breakoint(name, hbv_reads_info_txt, '-v secondary')[1]
            out_list = [name, ','.join(bk_loc_list), ','.join(hbv_seq_list)]
            out.write('\t'.join(out_list) + '\n')
            out.flush()
            # out bed
            for bk_loc in bk_loc_list:
                if bk_loc != 'N':
                    bk_loc_info = bk_loc.split('|')
                    chrom = bk_loc_info[0].split(':')[0]
                    start = bk_loc_info[0].split(':')[1].split('-')[0]
                    end = bk_loc_info[0].split(':')[1].split('-')[1]
                    out_list2 = [chrom, start, end, bk_loc_info[1], sampleid, name]
                    out_bed.write('\t'.join(out_list2) + '\n')
    out_bed.close()

########################################################################################
# SV type test
########################################################################################
# cat ../../output/*/breakpoint2.0/*minimap2_breakpoint_hm_hbv_secondary.txt|grep -v "N,N" > all_samples_minimap2_breakpoint_hm_hbv_secondary.txt
rd_bed = os.path.join(work_dir, 'breakpoint2.0/SV_class', 'all_samples_minimap2_breakpoint_hm_hbv_secondary.txt')
rd_bed2 = os.path.join(work_dir, 'breakpoint2.0/SV_class', 'all_samples_minimap2_breakpoint_hm_hbv_secondary2.txt')
with open(rd_bed2, 'w') as out:
    with open(rd_bed, 'r') as f:
        for line in f:
            c = line.strip().split()
            al_list = c[1].split(',')
            if 'N' in al_list:
                out_list = c + ['Single']
            else:
                chrom1 = al_list[0].split(':')[0]
                chrom2 = al_list[1].split(':')[0]
                if chrom1 == chrom2:
                    out_list = c + ['SV']
                else:
                    out_list = c + ['TRA']
            out.write('\t'.join(out_list) + '\n')

########################################################################################
# check with SV, HCC9_T1, HCC10_T4
########################################################################################
# /data/fs08/wangzf/nanopore/ztf/HCC/ONT/Somatic_2.4/HCC10_T4
# sniffles/HCC10_T4_merge_sniffles_v1.vcf
# sniffles/jasmine/HCC10_T4_merge_sniffles_v1_somatic_s3.vcf
# sniffles/jasmine/HCC10_T4_merge_sniffles_v1_somatic_s3_jasmine_filter1.vcf

# sniffles2/HCC10_T4_merge_sniffles_v2.vcf
# sniffles2/jasmine/HCC10_T4_merge_sniffles_v2_somatic_s3.vcf
# sniffles2/jasmine/HCC10_T4_merge_sniffles_v2_somatic_s3_jasmine_filter1.vcf

# /data/fs08/wangzf/nanopore/ztf/HCC/ONT/HBV_minimap2/output
# HCC8_T3/breakpoint2.0/HCC8_T3_minimap2_breakpoint_hm_hbv_secondary_merge.bed
# HCC8_T3/breakpoint2.0/HCC8_T3_minimap2_breakpoint_hm_hbv_secondary_sort.bed
# HCC8_T3/breakpoint2.0/HCC8_T3_minimap2_breakpoint_hm_hbv_secondary.txt
# /data/fs08/wangzf/nanopore/ztf/HCC/ONT/HBV_minimap2/breakpoint2.0/HBV_integration_sequence_length_per_read.txt

# /data/fs08/wangzf/nanopore/ztf/HCC/ONT/Somatic_3.0
# HCC10_T4/sniffles2.1/HCC10_T4_merge_sniffles2.vcf
# HCC10_T4/sniffles2.1/jasmine/HCC10_T4_merge_sniffles2_somatic_s3_jasmine.vcf

def bam_extarct(bam_in, bam_out, region_x):
    # os.system("{samtools} view -hb {bam_in} {region} > {bam_out}".format(
    #     samtools=samtools, bam_in=bam_in, bam_out=bam_out, region=region_x))
    os.system("{samtools} view -hb  {bam_in} {region} > {bam_out}".format(
        samtools=samtools, bam_in=bam_in, bam_out=bam_out, region=region_x))
    os.system("%s index %s" % (samtools, bam_out))

for sampleid1 in ['HCC8_T3']:
    sample_bk_dir = os.path.join(work_dir, 'output', sampleid1, "breakpoint2.0")
    hbv_bk_hm_bed = os.path.join(sample_bk_dir, "%s_minimap2_breakpoint_hm_hbv_secondary_merge.bed" % sampleid1)
    with open(hbv_bk_hm_bed, 'r') as f:
        for line in f:
            c = line.strip().split('\t')
            for sampleid in ['HCC8_WBC', 'HCC8_N1', 'HCC8_N3', 'HCC8_T1', 'HCC8_T2', 'HCC8_T3', 'HCC8_T4', 'HCC8_T5']:
                if sampleid in sample_n50_list:
                    bam_sort = os.path.join(ns_dir, 'output/6_QC', sampleid, '%s_ngmlr_sorted.bam' % sampleid)
                else:
                    bam_sort = os.path.join(ns_dir, 'output', sampleid, '%s_ngmlr_sorted.bam' % sampleid)
                region1 = c[0] + ':' + str(int(c[1]) - 20000) + '-' + str(int(c[2]) + 20000)
                # for region in ["chr2:74263683-74283684", "chr6:70290829-70310830"]:
                for region in [region1]:
                    # ngmlr to hg38
                    region_bam = os.path.join(work_dir, 'IGV', sampleid.split('_')[0],
                                              '%s_ngmlr_HBV.all.%s.bam' % (sampleid, region.replace(':', '-')))
                    print(region)
                    # bam_extarct(bam_sort, region_bam, region)
                    # minimap2 to HBV+hg38
                    bam_m2_newref = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_sorted.bam" % sampleid)
                    region_bam2 = os.path.join(work_dir, 'IGV', sampleid.split('_')[0],
                                               '%s_%s_minimap2_newref_all.bam' % (region.replace(':', '-'), sampleid))
                    print(region_bam2)
                    bam_extarct(bam_m2_newref, region_bam2, region)
                    # ngmlr to HBV+hg38
                    bam_ng_newref = os.path.join('/data/fs08/wangzf/nanopore/ztf/HCC/ONT/HBV_ngmlr', 'output', sampleid, "%s_ngmlr_sorted.bam" % sampleid)
                    region_bam3 = os.path.join(work_dir, 'IGV', sampleid.split('_')[0],
                                           '%s_ngmlr_newref.all.%s.bam' % (sampleid, region.replace(':', '-')))
                    # bam_extarct(bam_ng_newref, region_bam3, region)

for sampleid in ['HCC9_T1']:
    if sampleid in sample_n50_list:
        bam_sort = os.path.join(ns_dir, 'output/6_QC', sampleid, '%s_ngmlr_sorted.bam' % sampleid)
    else:
        bam_sort = os.path.join(ns_dir, 'output', sampleid, '%s_ngmlr_sorted.bam' % sampleid)
    for region in ["chr2:74263683-74283684", "chr6:70290829-70310830"]:
        # ngmlr to hg38
        region_bam = os.path.join(work_dir, 'IGV', 'HCC9',
                                  '%s_ngmlr_HBV.all.%s.bam' % (sampleid, region.replace(':', '-')))
        print(region)
        # bam_extarct(bam_sort, region_bam, region)
        # minimap2 to HBV+hg38
        bam_m2_newref = os.path.join(work_dir, 'output', sampleid, "%s_minimap2_sorted.bam" % sampleid)
        region_bam2 = os.path.join(work_dir, 'IGV', 'HCC9',
                                   '%s_minimap2_newref.all.%s.bam' % (sampleid, region.replace(':', '-')))
        bam_extarct(bam_m2_newref, region_bam2, region)
        # ngmlr to HBV+hg38
        bam_ng_newref = os.path.join('/data/fs08/wangzf/nanopore/ztf/HCC/ONT/HBV_ngmlr', 'output', sampleid,
                                     "%s_ngmlr_sorted.bam" % sampleid)
        region_bam3 = os.path.join(work_dir, 'IGV', 'HCC9',
                                   '%s_ngmlr_newref.all.%s.bam' % (sampleid, region.replace(':', '-')))
        # bam_extarct(bam_ng_newref, region_bam3, region)

# clliped reads BLAST
bam1 = "/data/fs08/wangzf/nanopore/ztf/HCC/ONT/HBV_minimap2/IGV/HCC9/HCC9_T1_ngmlr_newref.all.chr2-74263683-74283684.bam"
seq1 = os.popen("%s view %s|grep e4bb1096-7d88-4629-9b79-3c52efb22d32|cut -f 10" % (samtools, bam1)).read().strip()
