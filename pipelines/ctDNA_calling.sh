#! /usr/bin/bash

version="v1.0"

# ----------------------------------- Description ------------------------------------- #
# Perform:                                                                              #
#   Trimming (fastp) +                                                                  #
#   Alignment (Sentieon-bwa + realignment)                                              #
#   Variant Calling (Sentieon TNscope + Varscan)                                        #
#   VCF Normalisation +                                                                 #
#   Remove low-quality variants +                                                       #
#   VEP annotation                                                                      #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash ctDNA_calling.sh [input_folder] [output_folder]                    #
#                                                                                       #
# 常规ctDNA实验必然会添加UMI以保证特异性；因此在参数中若不提供 umi_template，则默认无法去重   #
# ctDNA 往往为小panel，故不进行 BQSR
# ------------------------------------------------------------------------------------- #


# --------------------------- set parameters --------------------------- #
# ---------------------------------------------------------------------- #

mode="single"
umi_template="6M+T,6M+T"

# path to software & scripts
fastp="/data/ngs/softs/fastp/fastp"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-201808.08/bin/sentieon"
bcftools="/public/software/bcftools-1.9/bcftools"
samtools="/public/software/samtools-1.14/samtools"
bamdst="/public/software/bamdst/bamdst"
vep="/public/software/98vep/ensembl-vep/vep"
bgzip="/public/software/htslib-1.9/bin/bgzip"
tabix="/public/software/htslib-1.9/bin/tabix"
VIC="/public/software/VIC"
annovar="/public/software/annovar/"

merge_mnv="/public/home/kai/BioinfoStuff/tertiary_analysis/merge_mnv.py"
anno_hgvs="/public/home/kai/BioinfoStuff/tertiary_analysis/anno_hgvs.py"
rm_common_variant="/public/home/kai/BioinfoStuff/rm_common_variant.py"

sentieon_license="192.168.1.186:8990"
thread=8

# additional files
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
dbsnp="/public/database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf"
k1="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
k2="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
k3="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
refflat="/public/database/GATK_Resource_Bundle/refFlat.txt"
vep_dir="/public/software/vep_98/"
cache_version="98"
clinic_transcripts="/public/home/kai/database/LRG/parsed_LRG.tsv"



# switch (on||off)
do_trim="on"
do_align="on"
do_realign="on"
do_qc="on"
do_tnscope="on"
do_mpileup="on"
do_varscan="on"
do_filt="on"
do_anno="on"


# ------------------------------ argparser ----------------------------- #
# ---------------------------------------------------------------------- #
if [[  $1 == '-h'  ]]; then
    echo "Usage: ./ctDNA_calling.sh [input_folder] [output_folder] [BED]"
    echo "-------------------------------------------------------------------------"
    echo "[input_folder] should contain fastq files with following naming system:"
    echo "  single -- \${sampleID}_R[1|2].fastq.gz"
    echo "  matched -- \${sampleID}_normal/tumor_R[1|2].fastq.gz"
    exit 0
fi

input_folder=`realpath $1`
output_folder=`realpath $2`

if [[  -f $bed  ]];
then
    bed=`realpath $3`
else
    bed="empty"
fi


if [[ ! -d $input_folder  ]];
then
    echo "Error: input_folder does not Found!"
    exit 1
fi



# ---------------  orgnise output directory structure  ----------------- #
# ---------------------------------------------------------------------- #
if [[ ! -d $output_folder ]]; 
then
    mkdir $output_folder
fi

trim_dir=$output_folder/trim/;
if [[ ! -d $trim_dir ]]; then
    mkdir $trim_dir
fi

align_dir=$output_folder/align/;
if [[ ! -d $align_dir ]]; then
    mkdir $align_dir
fi

snv_dir=${output_folder}/SNV/;
if [[  ! -d $snv_dir  ]];
then
    tnscope_dir=$snv_dir/tnscope
    varscan_dir=$snv_dir/varscan

    mkdir -p $tnscope_dir
    mkdir -p $varscan_dir
fi

qc_dir=$output_folder/qc/;
if [[ ! -d $qc_dir ]]; then
    mkdir $qc_dir
fi

# ---------------------------  LOGGING  -------------------------------- #
# ---------------------------------------------------------------------- #
echo "LOGGING: `date --rfc-3339=seconds` -- Analysis started"
echo "LOGGING: This is the ctDNA_calling.sh pipeline"
echo "========================================================"
echo "LOGGING: -- settings -- input folder -- ${input_folder}"
echo "LOGGING: -- settings -- output folder -- ${output_folder}"
echo "LOGGING: -- settings -- BED file -- ${bed}"
echo "========================================================"


echo "sampleID,fastq_size,raw_reads,raw_bases,clean_reads,clean_bases,\
qc30(%),trim_percent(%),mapping_rate(%),on-target_percent(%),\
mean_depth,mean_dedup_depth,dup_rate(%),\
average_insert_size,std_insert_size,\
Uniformity_0.1X(%),Uniformity_0.2X(%),\
Uniformity_0.5X(%),Uniformity_1X(%),\
50x_depth_percent(%),100x_depth_percent(%),\
150x_depth_percent(%),200x_depth_percent(%),\
300x_depth_percent(%),400x_depth_percent(%),\
500x_depth_percent(%)" \
> $qc_dir/QC_summary.csv


# ---------------------------------------------------------------------- #
# ---------------------------  Pipeline  ------------------------------- #
# ---------------------------------------------------------------------- #

export SENTIEON_LICENSE=${sentieon_license};

if [[  $mode == 'matched' ]]; then
    for ifile in $input_folder/*_normal_R1.fastq.gz;
    do 
        sampleID=`basename ${ifile%%"_normal"*}`

        # step1 - trim reads
        # remove reads with length < 50bp;
        # remove reads with >30% of bases below Q30;
        # remove low-complexity reads;
        # remove polyG tail; polyX tails;
        # remove read pairs that contain mismatches in overlapped regions

        if [[  $do_trim == "on" ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- trimming reads";

            $fastp --in1 $input_folder/${sampleID}_tumor_R1.fastq.gz \
            --in2 $input_folder/${sampleID}_tumor_R2.fastq.gz \
            --out1 $trim_dir/${sampleID}_tumor_trim_R1.fastq.gz \
            --out2 $trim_dir/${sampleID}_tumor_trim_R2.fastq.gz \
            --length_required 50 -y \
            -u 30 -q 30 -x -g \
            --overlap_diff_limit 0 \
            --detect_adapter_for_pe -p \
            --thread ${thread} \
            --html $trim_dir/${sampleID}_tumor.trim.html \
            --json $trim_dir/${sampleID}_tumor.trim.json;

            $fastp --in1 $input_folder/${sampleID}_normal_R1.fastq.gz \
            --in2 $input_folder/${sampleID}_normal_R2.fastq.gz \
            --out1 $trim_dir/${sampleID}_normal_trim_R1.fastq.gz \
            --out2 $trim_dir/${sampleID}_normal_trim_R2.fastq.gz \
            --length_required 50 -y \
            -u 30 -q 30 -x -g \
            --overlap_diff_limit 0 \
            --detect_adapter_for_pe -p \
            --thread ${thread} \
            --html $trim_dir/${sampleID}_normal.trim.html \
            --json $trim_dir/${sampleID}_normal.trim.json;
        fi

        # step2 - align & sort
        if [[  $do_align == "on"  ]] && [[ ! -n $umi_template  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- alignment & sorting";

            ($sentieon bwa mem -M -R "@RG\tID:${sampleID}_tumor\tSM:${sampleID}_tumor\tPL:illumina" \
            -t ${thread} -K 10000000 ${ref} \
            $trim_dir/${sampleID}_tumor_trim_R1.fastq.gz \
            $trim_dir/${sampleID}_tumor_trim_R2.fastq.gz \
            || echo -n 'error' ) \
            | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}_tumor.sorted.bam \
            -t ${thread} --sam2bam -i -;

            ($sentieon bwa mem -M -R "@RG\tID:${sampleID}_normal\tSM:${sampleID}_normal\tPL:illumina" \
            -t ${thread} -K 10000000 ${ref} \
            $trim_dir/${sampleID}_normal_trim_R1.fastq.gz \
            $trim_dir/${sampleID}_normal_trim_R2.fastq.gz \
            || echo -n 'error' ) \
            | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}_normal.sorted.bam \
            -t ${thread} --sam2bam -i -;
        fi

        tumor_bam=${align_dir}/${sampleID}_tumor.sorted.bam
        normal_bam=${align_dir}/${sampleID}_normal.sorted.bam


        # step2 - umi consensus & align
        
        if [[  $do_align == "on"  ]] && [[ -n $umi_template ]]; 
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- deduplication/umi consensus";

            $sentieon umi $umi_template $trim_dir/${sampleID}_tumor_trim_R1.fastq.gz $trim_dir/${sampleID}_tumor_trim_R2.fastq.gz | \
            ($sentieon bwa mem -t ${thread} -p -C -M -R "@RG\\tID:${sampleID}_tumor\\tSM:${sampleID}_tumor\\tPL:ILLUMINA\\tLB:lib" -K 10000000 ${ref} - || echo -n 'error') | \
            $sentieon util sort -t ${thread} --sam2bam -i - -o ${align_dir}/${sampleID}_tumor.sorted.bam;

            $sentieon umi consensus -i ${align_dir}/${sampleID}_tumor.sorted.bam --input_format BAM -o ${trim_dir}/${sampleID}_tumor.umi_consensus.fq.gz;

            ($sentieon bwa mem -t ${thread} -p -C -M -R "@RG\\tID:${sampleID}_tumor\\tSM:${sampleID}_tumor\\tPL:ILLUMINA\\tLB:lib" -K 10000000 ${ref} ${trim_dir}/${sampleID}_tumor.umi_consensus.fq.gz || echo -n 'error' ) | \
            $sentieon util sort --umi_post_process --sam2bam -i - -o ${align_dir}/${sampleID}_tumor.sorted.dedup.bam


            $sentieon umi $umi_template $trim_dir/${sampleID}_normal_trim_R1.fastq.gz $trim_dir/${sampleID}_normal_trim_R2.fastq.gz | \
            ($sentieon bwa mem -t ${thread} -p -C -M -R "@RG\\tID:${sampleID}_normal\\tSM:${sampleID}_normal\\tPL:ILLUMINA\\tLB:lib" -K 10000000 ${ref} - || echo -n 'error') | \
            $sentieon util sort -t ${thread} --sam2bam -i - -o ${align_dir}/${sampleID}_normal.sorted.bam;

            $sentieon umi consensus -i ${align_dir}/${sampleID}_normal.sorted.bam --input_format BAM -o ${trim_dir}/${sampleID}_normal.umi_consensus.fq.gz;

            ($sentieon bwa mem -t ${thread} -p -C -M -R "@RG\\tID:${sampleID}_normal\\tSM:${sampleID}_normal\\tPL:ILLUMINA\\tLB:lib" -K 10000000 ${ref} ${trim_dir}/${sampleID}_normal.umi_consensus.fq.gz || echo -n 'error' ) | \
            $sentieon util sort --umi_post_process --sam2bam -i - -o ${align_dir}/${sampleID}_normal.sorted.dedup.bam

        fi

        if [[ -n $umi_template ]];
        then
            tumor_bam=${align_dir}/${sampleID}_tumor.sorted.dedup.bam
            normal_bam=${align_dir}/${sampleID}_normal.sorted.dedup.bam
        fi

        # step4 - quality control
        if [[  $do_qc == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- qc";

            if [[ ! -d $qc_dir/${sampleID}_normal ]]; then
                mkdir $qc_dir/${sampleID}_normal
            fi

            if [[ ! -d $qc_dir/${sampleID}_tumor ]]; then
                mkdir $qc_dir/${sampleID}_tumor
            fi

            if [[  $bed != "empty"  ]];
            then
                $bamdst -p $bed -o $qc_dir/${sampleID}_normal \
                ${align_dir}/${sampleID}_normal.sorted.dedup.bam;

                $bamdst -p $bed -o $qc_dir/${sampleID}_tumor \
                ${align_dir}/${sampleID}_tumor.sorted.dedup.bam;
            fi

            $samtools stats -@ ${thread} ${align_dir}/${sampleID}_normal.sorted.dedup.bam > ${qc_dir}/${sampleID}_normal.stats.txt;
            $samtools stats -@ ${thread} ${align_dir}/${sampleID}_tumor.sorted.dedup.bam > ${qc_dir}/${sampleID}_tumor.stats.txt;

            tumor_r1=$(du $input_folder/${sampleID}_tumor_R1.fastq.gz -shL |awk '{print $1}');
            tumor_r2=$(du $input_folder/${sampleID}_tumor_R2.fastq.gz -shL |awk '{print $1}');
            normal_r1=$(du $input_folder/${sampleID}_normal_R1.fastq.gz -shL |awk '{print $1}');
            normal_r2=$(du $input_folder/${sampleID}_normal_R2.fastq.gz -shL |awk '{print $1}');

            normal_raw_reads=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}_normal.trim.json', 'r')); \
            print(fh['summary']['before_filtering']['total_reads'])"`

            normal_clean_reads=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}_normal.trim.json', 'r')); \
            print(fh['summary']['after_filtering']['total_reads'])"`

            normal_raw_bases=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}_normal.trim.json', 'r')); \
            print(fh['summary']['before_filtering']['total_bases'])"`

            normal_clean_bases=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}_normal.trim.json', 'r')); \
            print(fh['summary']['after_filtering']['total_bases'])"`

            tumor_raw_reads=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}_tumor.trim.json', 'r')); \
            print(fh['summary']['before_filtering']['total_reads'])"`

            tumor_clean_reads=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}_tumor.trim.json', 'r')); \
            print(fh['summary']['after_filtering']['total_reads'])"`

            tumor_raw_bases=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}_tumor.trim.json', 'r')); \
            print(fh['summary']['before_filtering']['total_bases'])"`

            tumor_clean_bases=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}_tumor.trim.json', 'r')); \
            print(fh['summary']['after_filtering']['total_bases'])"`

            normal_qc_rate=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}_normal.trim.json', 'r')); \
            print(fh['summary']['before_filtering']['q30_rate'])"`

            tumor_qc_rate=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}_tumor.trim.json', 'r')); \
            print(fh['summary']['before_filtering']['q30_rate'])"`

            normal_mapping_rate=$(grep "Fraction of Mapped Reads" $qc_dir/${sampleID}_normal/coverage.report | awk -F"\t" '{print $2}');
            tumor_mapping_rate=$(grep "Fraction of Mapped Reads" $qc_dir/${sampleID}_tumor/coverage.report | awk -F"\t" '{print $2}');

            normal_mean_depth=$(grep "Average depth" $qc_dir/${sampleID}_normal/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
            normal_mean_dedup_depth=$(grep "Average depth(rmdup)" $qc_dir/${sampleID}_normal/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
            normal_dup_rate=$(grep "Fraction of PCR duplicate reads" $qc_dir/${sampleID}_normal/coverage.report |awk -F"\t" '{print $2}');
            
            tumor_mean_depth=$(grep "Average depth" $qc_dir/${sampleID}_tumor/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
            tumor_mean_dedup_depth=$(grep "Average depth(rmdup)" $qc_dir/${sampleID}_tumor/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
            tumor_dup_rate=$(grep "Fraction of PCR duplicate reads" $qc_dir/${sampleID}_tumor/coverage.report |awk -F"\t" '{print $2}');

            normal_on_target=$(grep "Fraction of Target Reads in all reads" $qc_dir/${sampleID}_normal/coverage.report |awk -F"\t" '{print $2}');
            tumor_on_target=$(grep "Fraction of Target Reads in all reads" $qc_dir/${sampleID}_tumor/coverage.report |awk -F"\t" '{print $2}');

            normal_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}_normal.stats.txt);
            normal_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}_normal.stats.txt);

            tumor_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}_tumor.stats.txt);
            tumor_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}_tumor.stats.txt);

            normal_50x=$(less -S $qc_dir/${sampleID}_normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 50) count+=1} END {print count/NR*100}');
            normal_100x=$(less -S $qc_dir/${sampleID}_normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 100) count+=1} END {print count/NR*100}');
            normal_150x=$(less -S $qc_dir/${sampleID}_normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 150) count+=1} END {print count/NR*100}');
            normal_200x=$(less -S $qc_dir/${sampleID}_normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 200) count+=1} END {print count/NR*100}');
            normal_300x=$(less -S $qc_dir/${sampleID}_normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 300) count+=1} END {print count/NR*100}');
            normal_400x=$(less -S $qc_dir/${sampleID}_normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 400) count+=1} END {print count/NR*100}');
            normal_500x=$(less -S $qc_dir/${sampleID}_normal/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 500) count+=1} END {print count/NR*100}');

            tumor_50x=$(less -S $qc_dir/${sampleID}_tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 50) count+=1} END {print count/NR*100}');
            tumor_100x=$(less -S $qc_dir/${sampleID}_tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 100) count+=1} END {print count/NR*100}');
            tumor_150x=$(less -S $qc_dir/${sampleID}_tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 150) count+=1} END {print count/NR*100}');
            tumor_200x=$(less -S $qc_dir/${sampleID}_tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 200) count+=1} END {print count/NR*100}');
            tumor_300x=$(less -S $qc_dir/${sampleID}_tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 300) count+=1} END {print count/NR*100}');
            tumor_400x=$(less -S $qc_dir/${sampleID}_tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 400) count+=1} END {print count/NR*100}');
            tumor_500x=$(less -S $qc_dir/${sampleID}_tumor/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 500) count+=1} END {print count/NR*100}');

            normal_01x=$(less -S $qc_dir/${sampleID}_normal/depth.tsv.gz | awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.1) count+=1} END {print count/NR*100}');
            normal_02x=$(less -S $qc_dir/${sampleID}_normal/depth.tsv.gz | awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.2) count+=1} END {print count/NR*100}');
            normal_05x=$(less -S $qc_dir/${sampleID}_normal/depth.tsv.gz | awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.5) count+=1} END {print count/NR*100}');
            normal_1x=$(less -S $qc_dir/${sampleID}_normal/depth.tsv.gz | awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth) count+=1} END {print count/NR*100}');

            tumor_01x=$(less -S $qc_dir/${sampleID}_tumor/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.1) count+=1} END {print count/NR*100}');
            tumor_02x=$(less -S $qc_dir/${sampleID}_tumor/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.2) count+=1} END {print count/NR*100}');
            tumor_05x=$(less -S $qc_dir/${sampleID}_tumor/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.5) count+=1} END {print count/NR*100}');
            tumor_1x=$(less -S $qc_dir/${sampleID}_tumor/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth) count+=1} END {print count/NR*100}');

            echo "${sampleID}.normal,${normal_r1}/${normal_r2},${normal_raw_reads},${normal_raw_bases},${normal_clean_reads},${normal_clean_bases},\
            ${normal_qc_rate},${normal_mapping_rate},${normal_on_target},${normal_mean_depth},${normal_mean_dedup_depth},\
            ${normal_dup_rate},${normal_insert_size},${normal_insert_std},${normal_01x},${normal_02x},${normal_05x},${normal_1x},\
            ${normal_50x},${normal_100x},${normal_150x},${normal_200x},${normal_300x},${normal_400x},${normal_500x}" \
            >> ${qc_dir}/QC_summary.csv

            echo "${sampleID}.tumor,${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},\
            ${tumor_qc_rate},${tumor_mapping_rate},${tumor_on_target},${tumor_mean_depth},${tumor_mean_dedup_depth},${tumor_dup_rate},\
            ${tumor_insert_size},${tumor_insert_std},${tumor_01x},${tumor_02x},${tumor_05x},${tumor_1x},\
            ${tumor_50x},${tumor_100x},${tumor_150x},${tumor_200x},${tumor_300x},${tumor_400x},${tumor_500x}" \
            >> ${qc_dir}/QC_summary.csv
        fi

        # step5 - Indel Realignment
        if [[  $do_realign == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- Indel Realignment";

            $sentieon driver -t ${thread} -r ${ref} \
            -i $tumor_bam --algo Realigner -k ${k1} -k ${k2} -k ${k3} \
            ${align_dir}/${sampleID}_tumor.realign.bam

            $sentieon driver -t ${thread} -r ${ref} \
            -i $normal_bam --algo Realigner -k ${k1} -k ${k2} -k ${k3} \
            ${align_dir}/${sampleID}_normal.realign.bam
            
        fi

        tumor_bam=${align_dir}/${sampleID}_tumor.realign.bam
        normal_bam=${align_dir}/${sampleID}_normal.realign.bam

        # step6 - variant calling
        if [[  $do_tnscope == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- variant calling (TNscope)";

            $sentieon driver -t ${thread} -r ${ref} \
            -i $tumor_bam -i $normal_bam \
            --algo TNscope \
            --tumor_sample ${sampleID}_tumor \
            --normal_sample ${sampleID}_normal \
            --dbsnp ${dbsnp} \
            --trim_soft_clip \
            --min_tumor_allele_frac 0.0002 \
            --filter_t_alt_frac 0.0002 \
            --clip_by_minbq 1 \
            --min_base_qual 30 \
            --max_fisher_pv_active 0.05 \
            --assemble_mode 4 \
            --min_init_tumor_lod 1.0 \
            --min_tumor_lod 2.0 \
            --resample_depth 1000000000 \
            ${tnscope_dir}/${sampleID}_tumor.raw.vcf;
        fi

        if [[ $do_mpileup == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- generate mpileup";

            if [[  $bed != "empty"  ]];
            then
                $samtools mpileup -B -Q 20 -C 50 -q 30 -d 500000 \
                -f $ref $tumor_bam > ${align_dir}/${sampleID}_tumor.realign.mpileup;
            else 
                $samtools mpileup -B -Q 20 -C 50 -q 30 -d 500000 \
                -f $ref $tumor_bam > ${align_dir}/${sampleID}_tumor.realign.mpileup;
            fi
        fi

        if [[ $do_varscan == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- variant calling (Varscan)";

            java -Xmx3550m -Xms3550m -Xmn2g -jar $varscan \
            mpileup2snp ${align_dir}/${sampleID}_tumor.realign.mpileup \
            --output-vcf 1 \
            --min-reads2 2 \
            --min-avg-qual 25 \
            --min-var-freq 0.0002 \
            --min-strandedness 0.0001 > \
            ${varscan_dir}/${sampleID}_tumor.raw.snp.vcf;

            java -Xmx3550m -Xms3550m -Xmn2g -jar $varscan \
            mpileup2indel ${align_dir}/${sampleID}_tumor.realign.mpileup \
            --output-vcf 1 \
            --min-reads2 2 \
            --min-avg-qual 25 \
            --min-var-freq 0.0002 \
            --min-strandedness 0.0001 > \
            ${varscan_dir}/${sampleID}_tumor.raw.indel.vcf;
        fi

        # step7 - normalisation + remove low quality variants
        if [[  $do_filt == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- normalise VCF + filter low-support";
            
            # remove SV
            grep -v "SVTYPE=BND" ${snv_dir}/${sampleID}.raw.vcf \
            >  ${snv_dir}/${sampleID}.step1_snv.vcf;

            # split multiallelic sites + left-alignment
            $bcftools norm -m -both -f ${ref} \
            ${snv_dir}/${sampleID}.step1_snv.vcf \
            -o ${snv_dir}/${sampleID}.step2_norm.vcf;

            # bgzip and make index file
            $bgzip ${snv_dir}/${sampleID}.step2_norm.vcf;
            $tabix -p vcf ${snv_dir}/${sampleID}.step2_norm.vcf.gz;

            # filter off-target variants
            $bcftools view -R $bed \
            ${snv_dir}/${sampleID}.step2_norm.vcf.gz \
            > ${snv_dir}/${sampleID}.step3_on_target.vcf;

            # filter low-support variants
            awk '($1 ~ /^#/) || ($7 ~ /PASS/)' ${snv_dir}/${sampleID}.step3_on_target.vcf > \
            ${snv_dir}/${sampleID}.step4_filter.vcf

            $bcftools filter -i "(FORMAT/AF[:0]) >= 0.0002" \
            ${snv_dir}/${sampleID}.step4_filter.vcf > \
            ${snv_dir}/${sampleID}.step5_filter.vcf;

            $bcftools filter -i "(FORMAT/AD[:0]+AD[:1]) >= 1000" \
            ${snv_dir}/${sampleID}.step5_filter.vcf > \
            ${snv_dir}/${sampleID}.step6_filter.vcf; 

            # merge MNV
            python3 $merge_mnv ${snv_dir}/${sampleID}.step6_filter.vcf ${ref} \
            -o $snv_dir/${sampleID}.step7_MNV_merged.vcf;
        fi


        # step8 - annotation
        if [[  $do_anno == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- VEP Annotation";

            $vep --dir ${vep_dir} --cache --offline --cache_version ${cache_version} \
            --assembly GRCh37 --format vcf --fa ${ref} --force_overwrite --vcf \
            --gene_phenotype --use_given_ref --refseq --check_existing \
            --hgvs --hgvsg --transcript_version --max_af \
            --vcf_info_field ANN -i ${snv_dir}/${sampleID}.step7_MNV_merged.vcf \
            -o ${snv_dir}/${sampleID}.step8_anno.vcf;

            python3 $anno_hgvs ${snv_dir}/${sampleID}.step8_anno.vcf \
            $clinic_transcripts $refflat -o $snv_dir/$sampleID.step9_anno.vcf;
        fi
    done

elif [[  $mode == 'single'  ]]; then
    for ifile in $input_folder/*_R1.fastq.gz;
    do
        sampleID=`basename ${ifile%%"_R1"*}`

        # step1 - trim reads
        if [[ $do_trim == "on" ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- trimming reads";

            $fastp --in1 $input_folder/${sampleID}_R1.fastq.gz \
            --in2 $input_folder/${sampleID}_R2.fastq.gz \
            --out1 $trim_dir/${sampleID}_trim_R1.fastq.gz \
            --out2 $trim_dir/${sampleID}_trim_R2.fastq.gz \
            --length_required 50 -y \
            -u 30 -q 30 -x -g \
            --overlap_diff_limit 0 \
            --detect_adapter_for_pe -p \
            --thread ${thread} \
            --html $trim_dir/${sampleID}.trim.html \
            --json $trim_dir/${sampleID}.trim.json;
        fi


        # step2 - align & sort
        if [[ $do_align == "on" ]] && [[ ! -n $umi_template  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- alignment & sorting";

            ($sentieon bwa mem -M -R "@RG\tID:${sampleID}\tSM:${sampleID}\tPL:illumina" \
            -t ${thread} -K 10000000 ${ref} \
            ${trim_dir}/${sampleID}_trim_R1.fastq.gz \
            ${trim_dir}/${sampleID}_trim_R2.fastq.gz \
            || echo -n 'error' ) \
            | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.sorted.bam \
            -t ${thread} --sam2bam -i -;
        fi


        # step2 - align & umi consensus
        if [[  $do_align == "on"  ]] && [[ -n $umi_template ]]; 
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- deduplication";

            $sentieon umi $umi_template $trim_dir/${sampleID}_tumor_trim_R1.fastq.gz $trim_dir/${sampleID}_trim_R2.fastq.gz | \
            ($sentieon bwa mem -t ${thread} -p -C -M -R "@RG\\tID:${sampleID}_tumor\\tSM:${sampleID}_tumor\\tPL:ILLUMINA\\tLB:lib" -K 10000000 ${ref} - || echo -n 'error') | \
            $sentieon util sort -t ${thread} --sam2bam -i - -o ${align_dir}/${sampleID}.sorted.bam;

            $sentieon umi consensus -i ${align_dir}/${sampleID}.sorted.bam --input_format BAM -o ${trim_dir}/${sampleID}.umi_consensus.fq.gz;

            ($sentieon bwa mem -t ${thread} -p -C -M -R "@RG\\tID:${sampleID}_tumor\\tSM:${sampleID}_tumor\\tPL:ILLUMINA\\tLB:lib" -K 10000000 ${ref} ${trim_dir}/${sampleID}.umi_consensus.fq.gz || echo -n 'error' ) | \
            $sentieon util sort --umi_post_process --sam2bam -i - -o ${align_dir}/${sampleID}.sorted.dedup.bam
        fi

        if [[ -n $umi_template ]];
        then
            bam=${align_dir}/${sampleID}.sorted.dedup.bam
        else
            bam=${align_dir}/${sampleID}.sorted.bam
        fi

        # step5 - quality control
        if [[  $do_qc == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- qc";

            if [[ ! -d $qc_dir/${sampleID} ]]; then
                mkdir $qc_dir/${sampleID};
            fi;
            
            $bamdst -p $bed -o $qc_dir/${sampleID} \
            ${align_dir}/${sampleID}.sorted.dedup.bam;

            $samtools stats -@ ${thread} ${align_dir}/${sampleID}.sorted.dedup.bam > ${qc_dir}/${sampleID}.stats.txt;

            tumor_r1=$(du $input_folder/${sampleID}_R1.fastq.gz -shL |awk '{print $1}');
            tumor_r2=$(du $input_folder/${sampleID}_R2.fastq.gz -shL |awk '{print $1}');

            tumor_raw_reads=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
            print(fh['summary']['before_filtering']['total_reads'])"`

            tumor_clean_reads=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
            print(fh['summary']['after_filtering']['total_reads'])"`

            tumor_raw_bases=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
            print(fh['summary']['before_filtering']['total_bases'])"`

            tumor_clean_bases=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
            print(fh['summary']['after_filtering']['total_bases'])"`

            tumor_qc_rate=`python3 -c "import json; \
            fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
            print(fh['summary']['before_filtering']['q30_rate'])"`

            tumor_mapping_rate=$(grep "Fraction of Mapped Reads" $qc_dir/${sampleID}/coverage.report | awk -F"\t" '{print $2}');
            
            tumor_mean_depth=$(grep "Average depth" $qc_dir/${sampleID}/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
            tumor_mean_dedup_depth=$(grep "Average depth(rmdup)" $qc_dir/${sampleID}/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
            tumor_dup_rate=$(grep "Fraction of PCR duplicate reads" $qc_dir/${sampleID}/coverage.report |awk -F"\t" '{print $2}');

            tumor_on_target=$(grep "Fraction of Target Reads in all reads" $qc_dir/${sampleID}/coverage.report |awk -F"\t" '{print $2}');

            tumor_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);
            tumor_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);

            tumor_50x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 50) count+=1} END {print count/NR*100}');
            tumor_100x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 100) count+=1} END {print count/NR*100}');
            tumor_150x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 150) count+=1} END {print count/NR*100}');
            tumor_200x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 200) count+=1} END {print count/NR*100}');
            tumor_300x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 300) count+=1} END {print count/NR*100}');
            tumor_400x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 400) count+=1} END {print count/NR*100}');
            tumor_500x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 500) count+=1} END {print count/NR*100}');

            tumor_01x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.1) count+=1} END {print count/NR*100}');
            tumor_02x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.2) count+=1} END {print count/NR*100}');
            tumor_05x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.5) count+=1} END {print count/NR*100}');
            tumor_1x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${tumor_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth) count+=1} END {print count/NR*100}');

            echo "${sampleID},${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},\
            ${tumor_qc_rate},${tumor_mapping_rate},${tumor_on_target},${tumor_mean_depth},${tumor_mean_dedup_depth},${tumor_dup_rate},\
            ${tumor_insert_size},${tumor_insert_std},${tumor_01x},${tumor_02x},${tumor_05x},${tumor_1x},\
            ${tumor_50x},${tumor_100x},${tumor_150x},${tumor_200x},${tumor_300x},${tumor_400x},${tumor_500x}" \
            >> $qc_dir/QC_summary.csv
        fi


        # step6 - Indel Realignment
        if [[  $do_realign == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- Indel Realignment";

            $sentieon driver -t ${thread} -r ${ref} \
            -i $bam --algo Realigner -k ${k1} -k ${k2} -k ${k3} \
            ${align_dir}/${sampleID}.realign.bam            
        fi

        bam=${align_dir}/${sampleID}.realign.bam 

        # step7 - variant calling
        if [[ $do_tnscope == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- variant calling";

            $sentieon driver -t ${thread} -r ${ref} \
            -i $tumor_bam -i $normal_bam \
            --algo TNscope \
            --tumor_sample ${sampleID}_tumor \
            --normal_sample ${sampleID}_normal \
            --dbsnp ${dbsnp} \
            --trim_soft_clip \
            --min_tumor_allele_frac 0.0002 \
            --filter_t_alt_frac 0.0002 \
            --clip_by_minbq 1 \
            --min_base_qual 30 \
            --max_fisher_pv_active 0.05 \
            --assemble_mode 4 \
            --min_init_tumor_lod 1.0 \
            --min_tumor_lod 2.0 \
            --resample_depth 1000000000 \
            ${tnscope_dir}/${sampleID}.raw.vcf;

        fi

        if [[ $do_mpileup == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- generate mpileup";

            if [[  $bed != "empty"  ]];
            then
                $samtools mpileup -B -Q 20 -C 50 -q 30 -d 500000 \
                -f $ref $tumor_bam > ${align_dir}/${sampleID}_tumor.realign.mpileup;
            else 
                $samtools mpileup -B -Q 20 -C 50 -q 30 -d 500000 -l $bed\
                -f $ref $tumor_bam > ${align_dir}/${sampleID}_tumor.realign.mpileup;
            fi
        fi

        if [[ $do_varscan == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- variant calling (Varscan)";

            java -Xmx3550m -Xms3550m -Xmn2g -jar $varscan \
            mpileup2snp ${align_dir}/${sampleID}_tumor.realign.mpileup \
            --output-vcf 1 \
            --min-reads2 2 \
            --min-avg-qual 25 \
            --min-var-freq 0.0002 \
            --min-strandedness 0.0001 > \
            ${varscan_dir}/${sampleID}.raw.snp.vcf;

            java -Xmx3550m -Xms3550m -Xmn2g -jar $varscan \
            mpileup2indel ${align_dir}/${sampleID}.realign.mpileup \
            --output-vcf 1 \
            --min-reads2 2 \
            --min-avg-qual 25 \
            --min-var-freq 0.0002 \
            --min-strandedness 0.0001 > \
            ${varscan_dir}/${sampleID}.raw.indel.vcf;
        fi

        # step8 - normalise & remove low-quality variants
        if [[ $do_filt == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- normalise VCF + filter low-support";

            # remove 'germline_risk' annotation
            $bcftools annotate -x FILTER/germline_risk \
            ${snv_dir}/${sampleID}.tumor.raw.vcf > \
            ${snv_dir}/${sampleID}.step1_deanno.vcf;

            # remove SV
            grep -v "SVTYPE=BND" ${snv_dir}/${sampleID}.step1_deanno.vcf \
            >  ${snv_dir}/${sampleID}.step2_snv.vcf;

            # split multiallelic sites + left-alignment
            $bcftools norm -m -both -f ${ref} \
            ${snv_dir}/${sampleID}.step2_snv.vcf \
            -o ${snv_dir}/${sampleID}.step3_norm.vcf;

            # bgzip and make index file
            $bgzip ${snv_dir}/${sampleID}.step3_norm.vcf;
            $tabix -p vcf ${snv_dir}/${sampleID}.step3_norm.vcf.gz;

            # filter off-target variants
            $bcftools view -R $bed \
            ${snv_dir}/${sampleID}.step3_norm.vcf.gz \
            > ${snv_dir}/${sampleID}.step4_on_target.vcf;

            # filter low-support variants
            grep "#\|PASS" ${snv_dir}/${sampleID}.step4_on_target.vcf > \
            ${snv_dir}/${sampleID}.step5_filter.vcf

            $bcftools filter -i "(FORMAT/AF[:0]) >= 0.05" \
            ${snv_dir}/${sampleID}.step5_filter.vcf > \
            ${snv_dir}/${sampleID}.step6_filter.vcf;

            $bcftools filter -i "(FORMAT/AD[:0]+AD[:1]) >= 250" \
            ${snv_dir}/${sampleID}.step6_filter.vcf > \
            ${snv_dir}/${sampleID}.step7_filter.vcf; 

            # merge MNVs
            python3 $merge_mnv ${snv_dir}/${sampleID}.step6_filter.vcf ${ref} \
            -o $snv_dir/${sampleID}.step7_MNV_merged.vcf;
        fi

        # step9 - annotation
        if [[  $do_anno == "on"  ]];
        then
            echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- VEP Annotation";

            $vep --dir ${vep_dir} --cache --offline --cache_version ${cache_version} \
            --assembly GRCh37 --format vcf --fa ${ref} --force_overwrite --vcf \
            --gene_phenotype --use_given_ref --refseq --check_existing \
            --hgvs --hgvsg --transcript_version --max_af \
            --vcf_info_field ANN -i ${snv_dir}/${sampleID}.step7_MNV_merged.vcf \
            -o ${snv_dir}/${sampleID}.step8_anno.vcf;

            python3 $anno_hgvs ${snv_dir}/${sampleID}.step8_anno.vcf \
            $clinic_transcripts $refflat -o $snv_dir/$sampleID.step9_anno.vcf;
        fi
    done
fi

echo "LOGGING: `date --rfc-3339=seconds` -- Analysis finished";
