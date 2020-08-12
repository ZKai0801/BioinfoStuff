#!/usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
#      Noted: You need to purchase sentieon's license for use of this pipeline.         #
# ------------------------------------------------------------------------------------- #
# Perform:                                                                              #
#   Trimming (fastp) +                                                                  #
#   Alignment (bwa) +                                                                   #
#   statistic analysis +                                                                #
#   Deduplication (optional) +                                                          #
#   SNV/INDEL calling (TNscope)                                                         #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash snv_calling_tumor_only_TNscope.sh [input_folder] [output_folder]   #
#                                                                                       #
# input_folder should contain pair-end sequenced Fastq files.                           #
# Each sample is expected to have two Fastq files, with naming convention of:           #
# ${sampleID}_R[1|2].fastq.gz                                                           #
# ------------------------------------------------------------------------------------- #
# 1. Depends on DNA capture methods (i.e. targeted amplicon based or hybrid capture     #
#    based), user can choose to either perform or not perform deduplication step.       #
# 2. PON is NOT used by default, user can add it themselves if they have one            #
# 3. BED file is also NOT used by default,                                              #
#    as BED file will only specify regions of variant calling.                          #
#    This aim can be achieve equally by filtering obtained VCF file based on            #
#    genomic coordinate.                                                                #
# ------------------------------------------------------------------------------------- #


input_folder=$1
output_folder=$2

if [[ ! -d $output_folder ]]; then
    mkdir $output_folder
fi
################# set parameters ##################
sentieon_license="xxxxxxx"
thread=16
# whether perform deduplicate step (true || false)
dedup=false
# path to software 
fastp="/data/ngs/softs/fastp/fastp"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-201808.08/bin/sentieon"
# additional files
# reference fasta file
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
# sentieon provided three additional vcf file that can be used in Base Recalibration
k1="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
k2="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
k3="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
# dbSNP database
dbsnp="/data/ngs/database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
####################################################


export SENTIEON_LICENSE=${sentieon_license};


for ifile in $input_folder/*_R1.fastq.gz
    do
    sampleID=`basename ${ifile%%"_R1"*}`

    # step1 - trim reads
    trim_dir=$output_folder/trim/;
    if [[ ! -d $trim_dir ]]; then
        mkdir $trim_dir
    fi
    $fastp --in1 $input_folder/${sampleID}_R1.fastq.gz \
    --in2 $input_folder/${sampleID}_R2.fastq.gz \
    --out1 $trim_dir/${sampleID}_trim_R1.fastq.gz \
    --out2 $trim_dir/${sampleID}_trim_R2.fastq.gz \
    -c --length_required 3 --detect_adapter_for_pe -p \
    --thread ${thread} \
    --html $trim_dir/${sampleID}.trim.html \
    --json $trim_dir/${sampleID}.trim.json;

    # step2 - align & sort
    align_dir=$output_folder/align/;
    if [[ ! -d $align_dir ]]; then
        mkdir $align_dir
    fi
    ($sentieon bwa mem -M -R "@RG\tID:${sampleID}\tSM:${sampleID}\tPL:illumina" \
    -t ${thread} -K 10000000 ${ref} \
    $trim_dir/${sampleID}_trim_R1.fastq.gz \
    $trim_dir/${sampleID}_trim_R2.fastq.gz \
    || echo -n 'error' ) \
    | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.sorted.bam \
    -t ${thread} --sam2bam -i -;

    # step3 - statistic analysis
    ${sentieon} driver -t ${thread} -r ${ref} -i ${align_dir}/${sampleID}.sorted.bam \
    --algo GCBias --summary ${align_dir}/${sampleID}.gc_summary.txt \
    ${align_dir}/${sampleID}.gc_metric.txt \
    --algo MeanQualityByCycle ${align_dir}/${sampleID}.mq_metric.txt \
    --algo QualDistribution ${align_dir}/${sampleID}.qd_metric.txt \
    --algo InsertSizeMetricAlgo ${align_dir}/${sampleID}.is_metric.txt \
    --algo AlignmentStat ${align_dir}/${sampleID}.aln_metric.txt;

    # step4 (optional) - remove duplicates
    if [ $dedup == true ]; then
        ${sentieon} driver -t ${thread} \
        -i ${align_dir}/${sampleID}.sorted.bam \
        --algo LocusCollector \
        --fun score_info ${align_dir}/${sampleID}.score.txt;

        ${sentieon} driver -t ${thread} \
        -i ${align_dir}/${sampleID}.sorted.bam \
        --algo Dedup --score_info ${align_dir}/${sampleID}.score.txt \
        --metrics ${align_dir}/${sampleID}.dedup_metrics.txt \
        ${align_dir}/${sampleID}.sorted.dedup.bam;
    fi

    # set path of the bam file
    if [ $dedup == true ]; then
        bam=${align_dir}/${sampleID}.sorted.dedup.bam
    else
        bam=${align_dir}/${sampleID}.sorted.bam
    fi

    # step5 - BQSR
    ${sentieon} driver -t ${thread} -r ${ref} \
    -i ${bam} --algo QualCal -k ${k1} -k ${k2} -k ${k3} \
    ${align_dir}/${sampleID}.recal.table;

    # step6 - variant calling (TNscope)
    snv_dir=$output_folder/snv/;
    if [[ ! -d $snv_dir ]]; then
        mkdir $snv_dir
    fi
    ${sentieon} driver -t ${thread} -r ${ref} \
    -i ${bam} -q ${align_dir}/${sampleID}.recal.table \
    --algo TNscope --tumor_sample ${sampleID} \
    --dbsnp ${dbsnp} ${snv_dir}/${sampleID}.tumor.raw.vcf;

done