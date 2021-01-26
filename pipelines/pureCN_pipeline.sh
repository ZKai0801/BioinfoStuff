#! /usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
# Perform:                                                                              #
#   1) prepare normal database                                                          #
#        Trimming (fastp) +                                                             #
#        Alignment (Sentieon-bwa) +                                                     #
#        Deduplication (optional) +                                                     #
#        Quality Control +                                                              #
#        Calculate coverage (PureCN) +                                                  #
#        Build PoN (PureCN)                                                             #
#                                                                                       #
#   2) run PureCN tumor-only pipeline                                                   #
#                                                                                       #
#   Remove low-quality variants +                                                       #
#   VEP annotation +                                                                    #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$ bash pureCN_pipeline.sh [mode] [input_folder] [output_folder] [bed] \  #
#                [normal_ref]                                                           #
#                                                                                       #
# ------------------------------------------------------------------------------------- #


# --------------------------- set parameters --------------------------- #
# ---------------------------------------------------------------------- #

# path to software & scripts
fastp="/data/ngs/softs/fastp/fastp"
PURECN="/public/software/R_3.5.1/lib64/R/library/PureCN/extdata"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-201808.08/bin/sentieon"
bcftools="/public/software/bcftools-1.9/bcftools"
samtools="/public/software/samtools-1.9/samtools"
bamdst="/public/software/bamdst/bamdst"
bgzip="/public/software/htslib-1.9/bin/bgzip"
tabix="/public/software/htslib-1.9/bin/tabix"
vep="/public/software/98vep/ensembl-vep/vep"
Rscript="/public/software/R_3.5.1/bin/Rscript"

merge_mnv="/public/home/kai/BioinfoStuff/tertiary_analysis/merge_mnv.py"
HRDecipher="/public/home/kai/softwares/HRDecipher/HRDecipher.py"

sentieon_license="192.168.1.186:8990"

thread=8

# additional files
# reference fasta file
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
dbsnp="/public/database/GATK_Resource_Bundle/dbsnp_138.hg19.vcf.gz"
k1="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
k2="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
k3="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
refflat="/public/database/GATK_Resource_Bundle/refFlat.txt"
vep_dir="/public/software/vep_98/"
cache_version="98"
clinic_transcripts="/public/home/kai/database/LRG/parsed_LRG.tsv"

db_header="/public/test_data/HRD_project/sources/DB_info.header"

# ------------------------------ argparser ----------------------------- #
# ---------------------------------------------------------------------- #
if [[  $1 == '-h'  ]]; then
    echo "Usage: ./pureCN_pipeline.sh [mode] [input_folder] [output_folder] [bed] [normal_ref]"
    echo "-------------------------------------------------------------------------"
    echo "[mode] should be either 'run' mode or 'prepare' mode"
    echo "  'prepare' mode -- generate PoN from a cohort of normal samples"
    echo "  'run' mode -- run pureCN pipeline"
    echo "[input_folder] should contain fastq files with following naming system:"
    echo "  \${sampleID}_R[1|2].fastq.gz"
    echo "[normal_ref] when 'run' mode is selected, you need to provide path to the previously generated normal reference"
    echo  "However, you do not need to provide [normal_ref] when you use 'prepare' mode"
    exit 0
fi

mode=$1
input_folder=`realpath $2`
output_folder=`realpath $3`
bed=`realpath $4`
normal_ref=$5


if [[  $mode != '-h' ]] && [[  $mode != 'prepare'  ]] && [[  $mode != 'run'  ]]; then
    echo "InputError: Only 'prepare' and 'run' are valid modes!"
    exit 0
fi


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

if [[ ! -d $normal_ref  ]] && [[  $mode == 'run'  ]];
then
    echo "Error: Normal reference cannot be Found!"
    exit 1
fi

# ----------------------  orgnise output dir  -------------------------- #
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

qc_dir=$output_folder/qc/;
if [[ ! -d $qc_dir ]]; then
    mkdir $qc_dir
fi

snv_dir=$output_folder/snv/;
if [[ ! -d $snv_dir ]] && [[  $mode == 'run'  ]]; then
    mkdir $snv_dir
fi

purecn_dir=$output_folder/purecn/;
if [[ ! -d $purecn_dir ]]; then
    mkdir $purecn_dir
fi

hrd_dir=$output_folder/HRD/;
if [[  ! -d $hrd_dir  ]]; then
    mkdir $hrd_dir
fi

# ---------------------------  LOGGING  -------------------------------- #
# ---------------------------------------------------------------------- #
echo "LOGGING: `date --rfc-3339=seconds` -- Analysis started"
echo "LOGGING: This is the pureCN_pipeline.sh"
echo "========================================================"
echo "LOGGING: -- settings -- Mode -- ${mode}"
echo "LOGGING: -- settings -- input folder -- ${input_folder}"
echo "LOGGING: -- settings -- output folder -- ${output_folder}"
echo "LOGGING: -- settings -- BED file -- ${bed}"
echo "========================================================"


echo "sampleID,fastq_size,raw_reads,raw_bases,clean_reads,clean_bases,\
qc30_rate,mapping_rate(%),on-target_percent(%),\
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


if [[  $mode == 'prepare'  ]];
then
    # step1 -- prepare normal BAM files
    echo "LOGGING: `date --rfc-3339=seconds` -- prepare BAM files"

    for ifile in $input_folder/*_R1.fastq.gz;
    do 
        sampleID=`basename ${ifile%%"_normal"*}`

        # step1 - trim reads
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- trimming reads";

        $fastp --in1 $input_folder/${sampleID}_normal_R1.fastq.gz \
        --in2 $input_folder/${sampleID}_normal_R2.fastq.gz \
        --out1 $trim_dir/${sampleID}_trim_R1.fastq.gz \
        --out2 $trim_dir/${sampleID}_trim_R2.fastq.gz \
        -c --length_required 3 --detect_adapter_for_pe -p \
        --thread ${thread} \
        --html $trim_dir/${sampleID}.trim.html \
        --json $trim_dir/${sampleID}.trim.json;

        # step2 - align & sort
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- alignment & sorting";

        ($sentieon bwa mem -M -R "@RG\tID:${sampleID}\tSM:${sampleID}\tPL:illumina" \
        -t ${thread} -K 10000000 ${ref} \
        ${trim_dir}/${sampleID}_trim_R1.fastq.gz \
        ${trim_dir}/${sampleID}_trim_R2.fastq.gz \
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

        # step5 - quality control
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
    done


    # step2 -- generate interval file
    echo "LOGGING: `date --rfc-3339=seconds` -- generate interval file"

    $Rscript $PURECN/IntervalFile.R --infile $bed \
    --fasta $ref --outfile $purecn_dir/intervals.txt \
    --offtarget --genome hg19 \
    --export $purecn_dir/baits_optimized.bed;


    # step3 -- calculate coverage
    echo "LOGGING: `date --rfc-3339=seconds` -- calculate GC-normalised coverages"
    
    ls $align_dir/*dedup.bam > $purecn_dir/normal_bams.list;

    cov_dir=$purecn_dir/coverages/;
    if [[ ! -d $cov_dir ]]; then
        mkdir $cov_dir
    fi

    $Rscript $PURECN/Coverage.R --bam $purecn_dir/normal_bams.list \
    --intervals $purecn_dir/intervals.txt \
    --cpu $thread \
    --outdir $cov_dir;

    # step4 -- build normal database
    echo "LOGGING: `date --rfc-3339=seconds` -- build normal database"

    ls $cov_dir/*loess.txt > $purecn_dir/normal_loess.list;

    ref_dir=$purecn_dir/normal_ref;
    if [[ ! -d $ref_dir ]]; then
        mkdir $ref_dir
    fi

    $Rscript $PURECN/NormalDB.R --outdir $ref_dir \
    --coveragefiles $purecn_dir/normal_loess.list \
    --genome hg19 -f;

elif [[  $mode == 'run'  ]];
then
    for ifile in $input_folder/*_R1.fastq.gz;
    do 
        sampleID=`basename ${ifile%%"_R1"*}`

        # step1 - trim reads
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- trimming reads";

        $fastp --in1 $input_folder/${sampleID}_R1.fastq.gz \
        --in2 $input_folder/${sampleID}_R2.fastq.gz \
        --out1 $trim_dir/${sampleID}_trim_R1.fastq.gz \
        --out2 $trim_dir/${sampleID}_trim_R2.fastq.gz \
        -c --length_required 3 --detect_adapter_for_pe -p \
        --thread ${thread} \
        --html $trim_dir/${sampleID}.trim.html \
        --json $trim_dir/${sampleID}.trim.json;

        # step2 - align & sort
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- alignment & sorting";

        ($sentieon bwa mem -M -R "@RG\tID:${sampleID}\tSM:${sampleID}\tPL:illumina" \
        -t ${thread} -K 10000000 ${ref} \
        ${trim_dir}/${sampleID}_trim_R1.fastq.gz \
        ${trim_dir}/${sampleID}_trim_R2.fastq.gz \
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

        # step4 - remove duplicates
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

        # step5 - quality control
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

        # step6 - BQSR
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- BQSR";

        ${sentieon} driver -t ${thread} -r ${ref} \
        -i ${align_dir}/${sampleID}.sorted.dedup.bam \
         --algo QualCal -k ${k1} -k ${k2} -k ${k3} \
        ${align_dir}/${sampleID}.recal.table;

        # step7 - variant calling
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- variant calling";

        ${sentieon} driver -t ${thread} -r ${ref} \
        -i ${align_dir}/${sampleID}.sorted.dedup.bam \
        -q ${align_dir}/${sampleID}.recal.table \
        --algo TNhaplotyper2 --tumor_sample ${sampleID} \
        --default_af 0.00000001 \
        ${snv_dir}/${sampleID}.tumor.tmp.vcf;

        ${sentieon} tnhapfilter --tumor_sample ${sampleID} \
        -v ${snv_dir}/${sampleID}.tumor.tmp.vcf \
        ${snv_dir}/${sampleID}.tumor.raw.vcf;

        rm ${snv_dir}/${sampleID}.tumor.tmp.vcf*;

        # step8 - normalise
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- normalise VCF + filter low-support";

        # split multiallelic sites + left-alignment
        $bcftools norm -m -both -f ${ref} \
        ${snv_dir}/${sampleID}.tumor.raw.vcf \
        -o ${snv_dir}/${sampleID}.step2_norm.vcf;

        # bgzip and make index file
        $bgzip ${snv_dir}/${sampleID}.step2_norm.vcf;
        $tabix -p vcf ${snv_dir}/${sampleID}.step2_norm.vcf.gz;

        # filter off-target variants
        $bcftools view -R $bed \
        ${snv_dir}/${sampleID}.step2_norm.vcf.gz \
        > ${snv_dir}/${sampleID}.step3_on_target.vcf;

        # filter low-support variants
        grep "#\|PASS" ${snv_dir}/${sampleID}.step3_on_target.vcf > \
        ${snv_dir}/${sampleID}.step4_filter.vcf

        $bcftools filter -i "(FORMAT/AF[0:0]) >= 0.05" \
        ${snv_dir}/${sampleID}.step4_filter.vcf > \
        ${snv_dir}/${sampleID}.step5_filter.vcf;

        $bcftools filter -i "(FORMAT/AD[0:0]+AD[0:1]) >= 250" \
        ${snv_dir}/${sampleID}.step5_filter.vcf > \
        ${snv_dir}/${sampleID}.step6_filter.vcf; 
  
        # merge MNVs
        python3 $merge_mnv ${snv_dir}/${sampleID}.step6_filter.vcf ${ref} \
        -o $snv_dir/${sampleID}.step7_MNV_merged.vcf;

        # annotate variants
        $bgzip $snv_dir/${sampleID}.step7_MNV_merged.vcf;
        $tabix -p vcf $snv_dir/${sampleID}.step7_MNV_merged.vcf.gz;

        $bcftools annotate -c "INFO/DB" -h $db_header -a $dbsnp \
        $snv_dir/${sampleID}.step7_MNV_merged.vcf.gz \
        -o $snv_dir/${sampleID}.step8_final.vcf;

        # step6 - run pureCN
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- calculate coverage";
        
        if [[ ! -d $purecn_dir/${sampleID} ]]; then
            mkdir $purecn_dir/${sampleID};
        fi;

        $Rscript $PURECN/Coverage.R --outdir $purecn_dir/$sampleID \
        --bam ${align_dir}/${sampleID}.sorted.dedup.bam \
        --intervals $normal_ref/intervals.txt;

        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- run PureCN";
        $Rscript $PURECN/PureCN.R --out $purecn_dir/$sampleID \
        --tumor $purecn_dir/${sampleID}/${sampleID}.sorted.dedup_coverage_loess.txt \
        --sampleid $sampleID \
        --vcf $snv_dir/${sampleID}.step8_final.vcf \
        --normaldb $normal_ref/normalDB_hg19.rds \
        --intervals $normal_ref/intervals.txt \
        --intervalweightfile $normal_ref/interval_weights_hg19.txt \
        --genome hg19 --cores $thread \
        --force --postoptimize --seed 100;

        # step7 - parse pureCN results and calculate HRD scores
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- calculate HRD scores";

        ploidy=`awk -F"," 'NR==2 {printf "%.1f", $3}' ${purecn_dir}/${sampleID}/${sampleID}.csv`

        awk -F"," -v ploidy=${ploidy} \
        'BEGIN {OFS="\t"; print "SampleID","Chromosome","Start_position","End_position","total_cn","A_cn","B_cn","ploidy"};
         NR!=1 {printf "%s\t%s\t%s\t%s\t%.0f\t%.0f\t%.0f\t%s\n", $1, $2, $3, $4, $6, $6-$7, $7, ploidy}' \
        ${purecn_dir}/${sampleID}/${sampleID}_loh.csv > $hrd_dir/${sampleID}.pre_hrd.tsv;
        
        python3 $HRDecipher $hrd_dir/${sampleID}.pre_hrd.tsv;
    done

    cat $hrd_dir/*.hrd.tsv | awk 'NR==1 || !/(HRD-sum)/' > $hrd_dir/HRD_results.tsv

fi


echo "LOGGING: `date --rfc-3339=seconds` -- Analysis finished"
