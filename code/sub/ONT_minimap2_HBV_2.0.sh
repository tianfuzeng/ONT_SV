#! /bin/bash
sampleID=$1
fastq=$2
Workdir=$3
ref=/data/fs09/wangzf/nanopore/ztf/HCC/ref/HBV_ztf/hg38_mainChr_and_JQ688404.1_C2_HBV_genome.fa
mappingDir=${Workdir}/${sampleID}
SVDir=${Workdir}/${sampleID}
#files configure
sam=$mappingDir/${sampleID}_minimap2.sam
bam=$mappingDir/${sampleID}_minimap2.bam
bam_sort=$mappingDir/${sampleID}_minimap2_sorted.bam
bam_sort_flagstat=$mappingDir/${sampleID}_minimap2_sorted.bam.flagstat
bam_sort_depth=$mappingDir/${sampleID}_minimap2_sorted_depth.bed
vcf1=$SVDir/${sampleID}_minimap2_sorted_sniffles1.vcf
vcf2=$SVDir/${sampleID}_minimap2_sorted_sniffles2.vcf
bedpe=$SVDir/${sampleID}_minimap2_sorted_sniffles.bedpe
#sam_paf=$mappingDir/ngmlr_$sampleID.paf
####ngmlr mapping####
echo "$(date) 1. Start to mapping on sample: $sampleID"
/data/fs01/wangzf/software/minimap2-2.24_x64-linux/minimap2 -N 10 -p 0.3 -ax map-ont --MD -t 30 $ref "$fastq" -o "$sam"
  echo "$(date) 1. Finish mapping on sample: $sampleID"
####sam to bam####
echo "$(date) 2. sam to bam: $sampleID"
/data/fs01/biosoft/samtools-1.9/samtools view -@ 10 -hb -S "$sam" >"$bam"
  echo "$(date) 2. sam to bam finish: $sampleID"
####sort bam####
echo "$(date) 3. sort bam: $sampleID"
/data/fs01/biosoft/samtools-1.9/samtools sort -@ 10 "$bam" >"$bam_sort"
/data/fs01/biosoft/samtools-1.9/samtools index -@ 10 "$bam_sort"
echo "$(date) 3. sort bam finish: $sampleID"
####rmsam###
rm "$sam"
rm "$bam"
