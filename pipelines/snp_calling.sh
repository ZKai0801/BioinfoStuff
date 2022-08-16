#! /usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
#      Noted: You need to purchase sentieon's license for use of this pipeline.         #
# ------------------------------------------------------------------------------------- #
# Perform:                                                                              #
#   Trimming (fastp) +                                                                  #
#   Alignment (bwa) +                                                                   #
#   statistic analysis +                                                                #
#   Deduplication (optional) +                                                          #
#   Quality Control +                                                                   #
#   Indel-rearrangement +                                                               #
#   BQSR +                                                                              #
#   SNP/INDEL calling (Haplotyper)  +                                                   #
#   Variant Position Normalisation (bcftools) +                                         #
#   Filter based on depth coverage (bcftools) +                                         #
#   Variant Annotation (Ensembl-VEP)                                                    #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash snp_sentieon.sh [input_folder] [output_folder] [BED_file]          #
#                                                                                       #
# input_folder should contain pair-end sequenced Fastq files.                           #
# Each sample is expected to have two Fastq files, with naming convention of:           #
# ${sampleID}_R[1|2].fastq.gz                                                           #
# ------------------------------------------------------------------------------------- #
# 1. Depends on DNA capture methods (i.e. targeted amplicon based or hybrid capture     #
#    based), user can choose to either perform or not perform deduplication step.       #
# ------------------------------------------------------------------------------------- #


# --------------------------- set parameters --------------------------- #
# ---------------------------------------------------------------------- #

sentieon_license="192.168.1.186:8990"
thread=16
# whether perform deduplicate step (true || false)
dedup=true

# path to software 
fastp="/data/ngs/softs/fastp/fastp"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-201808.08/bin/sentieon"
bcftools="/public/software/bcftools-1.9/bcftools"
bgzip="/public/software/htslib-1.9/bin/bgzip"
tabix="/public/software/htslib-1.9/bin/tabix"
vep="/public/software/98vep/ensembl-vep/vep"
bamdst="/public/software/bamdst/bamdst"

# additional files
# reference fasta file
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
# sentieon provided three additional vcf file that can be used in Base Recalibration
k1="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
k2="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
k3="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
# dbSNP database
dbsnp="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
# vep database
vep_dir="/public/software/vep_98/"
cache_version="98"


# switch (on||off)
do_trim="on"
do_align="on"
do_dedup="on"
do_qc="on"
do_realign="on"
do_bqsr="on"
do_snp="on"
do_anno="on"

# ------------------------------ argparser ----------------------------- #
# ---------------------------------------------------------------------- #

if [[  $1 == '-h'  ]]; then
    echo "Usage: ./snp_calling.sh [input_folder] [output_folder] [BED]"
    echo "-------------------------------------------------------------------------"
    echo "[input_folder] should contain fastq files with following naming system:"
    echo "  \${sampleID}_R[1|2].fastq.gz"
    exit 0
fi

input_folder=`realpath $1`
output_folder=`realpath $2`
bed=`realpath $3`

if [[ ! -d $input_folder  ]];
then
    echo "Error: input_folder does not Found!"
    exit 1
fi


if [[ ! -f $bed  ]];
then
    echo "Error: BED file does not Found!"
    exit 1
fi


if [[ ! -d $output_folder ]]; 
then
    mkdir $output_folder
fi


# ----------------------  orgnise output dir  -------------------------- #
# ---------------------------------------------------------------------- #

trim_dir=$output_folder/trim/;
if [[ ! -d $trim_dir ]]; then
    mkdir $trim_dir
fi

align_dir=$output_folder/align/;
if [[ ! -d $align_dir ]]; then
    mkdir $align_dir
fi

qc_dir=$output_folder/qc/;
if [[ ! -d $qc_dir ]]; then
    mkdir $qc_dir
fi

snp_dir=$output_folder/snp/;
if [[ ! -d $snp_dir ]]; then
    mkdir $snp_dir
fi

# ---------------------------  LOGGING  -------------------------------- #
# ---------------------------------------------------------------------- #

echo "LOGGING: `date --rfc-3339=seconds` -- Analysis started"
echo "LOGGING: This is the snp_calling.sh pipeline"
echo "========================================================"
echo "LOGGING: -- settings -- input folder -- ${input_folder}"
echo "LOGGING: -- settings -- output folder -- ${output_folder}"
echo "LOGGING: -- settings -- BED file -- ${bed}"
echo "========================================================"

if [[  $do_qc == "on"  ]]; then
echo "sampleID,fastq_size,raw_reads,raw_bases,clean_reads,clean_bases,\
qc30_rate,mapping_rate(%),on-target_percent(%),\
mean_depth,mean_dedup_depth,dup_rate(%),\
average_insert_size,std_insert_size,\
Uniformity_0.1X(%),Uniformity_0.2X(%),\
Uniformity_0.5X(%),Uniformity_1X(%),\
50x_depth_percent(%),100x_depth_percent(%),\
150x_depth_percent(%),200x_depth_percent(%),\
300x_depth_percent(%),400x_depth_percent(%),\
500x_depth_percent(%)" > $qc_dir/QC_summary.csv
fi

# ---------------------------------------------------------------------- #
# ---------------------------  Pipeline  ------------------------------- #
# ---------------------------------------------------------------------- #

export SENTIEON_LICENSE=${sentieon_license};

for ifile in $input_folder/*_R1.fastq.gz;
do 
    sampleID=`basename ${ifile%%"_R1"*}`

    # step1 - trim reads
    if [[  $do_trim == "on"  ]]; 
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- trimming";

        $fastp --in1 $input_folder/${sampleID}_R1.fastq.gz \
        --in2 $input_folder/${sampleID}_R2.fastq.gz \
        --out1 $trim_dir/${sampleID}_trim_R1.fastq.gz \
        --out2 $trim_dir/${sampleID}_trim_R2.fastq.gz \
        -c -q 25 \
        --length_required 50 --detect_adapter_for_pe -p \
        --thread ${thread} \
        --html $trim_dir/${sampleID}.trim.html \
        --json $trim_dir/${sampleID}.trim.json;
    fi

    # step2 - align & sort
    if [[  $do_align == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- alignment & sorting";

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
    fi

    # step4 (optional) - remove duplicates
    if [[  $do_dedup == "on"  ]] && [[ $dedup == true ]]; then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- deduplication";

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

    # step5 - quality control
    if [[  $do_qc == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- quality control";

        if [[ ! -d $qc_dir/${sampleID} ]]; then
            mkdir $qc_dir/${sampleID}
        fi

        $bamdst -p $bed -o $qc_dir/${sampleID} \
        ${align_dir}/${sampleID}.sorted.dedup.bam;

        $samtools stats -@ ${thread} ${align_dir}/${sampleID}.sorted.dedup.bam > ${qc_dir}/${sampleID}.stats.txt;

        normal_r1=$(du $input_folder/${sampleID}_R1.fastq.gz -shL |awk '{print $1}');
        normal_r2=$(du $input_folder/${sampleID}_R2.fastq.gz -shL |awk '{print $1}');

        normal_raw_reads=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['total_reads'])"`

        normal_clean_reads=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['total_reads'])"`

        normal_raw_bases=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['total_bases'])"`

        normal_clean_bases=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['total_bases'])"`

        normal_qc_rate=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['q30_rate'])"`

        normal_mapping_rate=$(grep "Fraction of Mapped Reads" $qc_dir/${sampleID}/coverage.report | awk -F"\t" '{print $2}');
        
        normal_mean_depth=$(grep "Average depth" $qc_dir/${sampleID}/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
        normal_mean_dedup_depth=$(grep "Average depth(rmdup)" $qc_dir/${sampleID}/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
        normal_dup_rate=$(grep "Fraction of PCR duplicate reads" $qc_dir/${sampleID}/coverage.report |awk -F"\t" '{print $2}');

        normal_on_target=$(grep "Fraction of Target Reads in all reads" $qc_dir/${sampleID}/coverage.report |awk -F"\t" '{print $2}');

        normal_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);
        normal_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);

        normal_50x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 50) count+=1} END {print count/NR*100}');
        normal_100x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 100) count+=1} END {print count/NR*100}');
        normal_150x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 150) count+=1} END {print count/NR*100}');
        normal_200x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 200) count+=1} END {print count/NR*100}');
        normal_300x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 300) count+=1} END {print count/NR*100}');
        normal_400x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 400) count+=1} END {print count/NR*100}');
        normal_500x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 500) count+=1} END {print count/NR*100}');

        normal_01x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.1) count+=1} END {print count/NR*100}');
        normal_02x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.2) count+=1} END {print count/NR*100}');
        normal_05x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.5) count+=1} END {print count/NR*100}');
        normal_1x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${normal_mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth) count+=1} END {print count/NR*100}');


        echo "${sampleID},${normal_r1}/${normal_r2},${normal_raw_reads},${normal_raw_bases},${normal_clean_reads},\
        ${normal_clean_bases},${normal_qc_rate},${normal_mapping_rate},${normal_on_target},\
        ${normal_mean_depth},${normal_mean_dedup_depth},${normal_dup_rate},\
        ${normal_insert_size},${normal_insert_std},${normal_01x},${normal_02x},${normal_05x},${normal_1x},\
        ${normal_50x},${normal_100x},${normal_150x},${normal_200x},${normal_300x},${normal_400x},${normal_500x}" \
        >> $qc_dir/QC_summary.csv
    fi

    # step6 - indel re-alignment
    if [[  $do_realign == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- indel-realignment";
        if [ $dedup == true ]; then
            ${sentieon} driver -t ${thread} -r ${ref} \
            -i ${align_dir}/${sampleID}.sorted.dedup.bam \
            --algo Realigner -k ${k1} -k ${k2} -k ${k3} \
            ${align_dir}/${sampleID}.realigned.bam
        else
            ${sentieon} driver -t ${thread} -r ${ref} \
            -i ${align_dir}/${sampleID}.sorted.bam \
            --algo Realigner -k ${k1} -k ${k2} -k ${k3} \
            ${align_dir}/${sampleID}.realigned.bam
        fi
    fi

    # step7 - BQSR
    if [[  $do_bqsr == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- BQSR";
        ${sentieon} driver -t ${thread} -r ${ref} \
        -i ${align_dir}/${sampleID}.realigned.bam \
        --algo QualCal -k ${k1} -k ${k2} -k ${k3} ${align_dir}/${sampleID}.recal.table;

        ${sentieon} driver -t ${thread} -r ${ref} -i ${align_dir}/${sampleID}.realigned.bam \
        -q ${align_dir}/${sampleID}.recal.table --algo QualCal -k ${k1} -k ${k2} -k ${k3} \
        ${align_dir}/${sampleID}.recal.post.table \
        --algo ReadWriter ${align_dir}/${sampleID}.recal.bam;
    fi

    # step8 - Haplotyper
    if [[  $do_snp == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- snp calling";

        ${sentieon} driver -r ${ref} -t ${thread} \
        -i ${align_dir}/${sampleID}.recal.bam \
        --algo Haplotyper --emit_conf=10 --call_conf=10 \
        -d ${dbsnp} \
        ${snp_dir}/${sampleID}.germline.raw.vcf;

        # step9 - normalise variants positions
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- VCF normalisation & filtration";

        $bcftools norm -m -both -f ${ref} \
        ${snp_dir}/${sampleID}.germline.raw.vcf \
        -o ${snp_dir}/${sampleID}.germline.step1_norm.vcf

        # step10 - filter variants without enough coverage
        $bcftools filter -i "FORMAT/DP>20" \
        ${snp_dir}/${sampleID}.germline.step1_norm.vcf > \
        ${snp_dir}/${sampleID}.germline.step2_filter.vcf;

        $bgzip ${snp_dir}/${sampleID}.germline.step2_filter.vcf;

        $tabix -p vcf ${snp_dir}/${sampleID}.germline.step2_filter.vcf.gz;

        $bcftools view -R $bed \
        ${snp_dir}/${sampleID}.germline.step2_filter.vcf.gz \
        > ${snp_dir}/${sampleID}.germline.step3_on_target.vcf;
    fi

    # step11 - annotation
    if [[  $do_anno == "on"  ]]; then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- variants annotation";

        $vep --dir ${vep_dir} --cache --offline --cache_version ${cache_version} \
        --assembly GRCh37 --format vcf --fa ${ref} --force_overwrite --vcf \
        --gene_phenotype --use_given_ref --refseq --check_existing \
        --hgvs --hgvsg --transcript_version --max_af \
        --vcf_info_field ANN -i ${snp_dir}/${sampleID}.germline.step3_on_target.vcf \
        -o ${snp_dir}/${sampleID}.germline.step4_anno.vcf;
    fi
done

echo "LOGGING: `date --rfc-3339=seconds` -- Analysis completed"
