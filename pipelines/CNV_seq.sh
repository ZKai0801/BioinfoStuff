#!/usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
# Perform:                                                                              #
#   Trimming (fastp) +                                                                  #
#   Alignment (Sentieon-bwa) +                                                          #
#   Deduplication +                                                                     #
#   statistic analysis +                                                                #
#   CNV calling (ichorCNA)                                                              #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash CNV_seq.sh [input_folder] [output_folder] [mode]                   #
#                                                                                       #
# input_folder should contain low-depth CNV-seq pair-end sequenced Fastqs               #
# Three modes could be choosed from:                                                    #
# - pair: used paried fastqs to calculate CNV portrait                                  #
# - single: no normal fastqs need to supplied.                                          #
# - pon: You need to provide a PoN generated by createPanelOfNormals.R                  #
#                                                                                       #
# Naming convention:                                                                    #
# tumor:   ${sampleID}_tumor_R[1|2].fastq.gz                                            #
# normal:  ${sampleID}_normal_R[1|2].fastq.gz                                           #
#                                                                                       #
# Note:                                                                                 #
# This script assume you use hg19 reference genome.                                     #
# ------------------------------------------------------------------------------------- #


# --------------------------- set parameters --------------------------- #
# ---------------------------------------------------------------------- #
sentieon_license="192.168.1.186:8990"
thread=8
dedup=true
low_tfx_setting="true"

# software
fastp="/data/ngs/softs/fastp/fastp"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-201808.08/bin/sentieon"
samtools="/public/software/samtools-1.9/samtools"
bamdst="/public/software/bamdst/bamdst"
ichorCNA="/public/software/ichorCNA/"
readCounter="/public/software/hmmcopy_utils/bin/readCounter"
compute_cin="/public/home/kai/BioinfoStuff/compute_cin_score.R"
Rscript="/public/software/R_3.6.3/bin/Rscript"

ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
bed="/public/database/ucsc/hg19_no_gaps.txt"

# PoN
pon="/public/home/kai/CNVseq_project/pon/pon_median.rds"

# switch (on||off)
do_trim="on"
do_align="on"
do_dedup="on"
do_qc="on"
do_readcount="on"
do_icorcna="on"

# low TFx setting
if [[  low_tfx_setting == "true"  ]];
then
    nfx="c(0.95, 0.99, 0.995, 0.999)"
    ploidy="c(2)"
    maxCN="3"
    subclone="c()"
    scpre="FALSE"
    chrs="c(1:22)"
else
    nfx="c(0.5,0.6,0.7,0.8,0.9)"
    ploidy="c(2,3)"
    maxCN="5"
    subclone="c(1,3)"
    scpre="FALSE"
    chrs='c(1:22)'
fi


# ------------------------------ argparser ----------------------------- #
# ---------------------------------------------------------------------- #
if [[  $1 == '-h'  ]]; then
    echo "Usage: ./CNV_seq.sh [input_folder] [output_folder] [mode]"
    echo "-------------------------------------------------------------------------"
    echo "Mode can be choosed from below three options:"
    echo "pair"
    echo "single"
    echo "pon"
    exit 0
fi

input_folder=`realpath $1`
output_folder=`realpath $2`
mode=$3

if [[ ! -d $input_folder  ]];
then
    echo "Error: input_folder does not Found!"
    exit 1
fi

if [[  $mode != "pair"  ]] && [[  $mode != "single"  ]] && [[  $mode != "pon"  ]];
then
    echo "You can only choose mode of 'pair', 'single' or 'pon' "
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

cnv_dir=$output_folder/cnv/;
if [[ ! -d $cnv_dir ]]; then
    mkdir $cnv_dir
fi

qc_dir=$output_folder/qc/;
if [[ ! -d $qc_dir ]]; then
    mkdir $qc_dir
fi


# ---------------------------  LOGGING  -------------------------------- #
# ---------------------------------------------------------------------- #
echo "LOGGING: `date --rfc-3339=seconds` -- Analysis started"
echo "LOGGING: This is the CNV_seq.sh pipeline"
echo "========================================================"
echo "LOGGING: -- settings -- input folder -- ${input_folder}"
echo "LOGGING: -- settings -- output folder -- ${output_folder}"
echo "LOGGING: -- settings -- mode -- ${mode}"
echo "LOGGING: -- settings -- low_tfx_setting -- ${low_tfx_setting}"
echo "========================================================"

if [[  $do_qc == "on"  ]];
then
echo "sampleID,fastq_size,raw_reads,raw_bases,clean_reads,clean_bases,qc30_rate,\
mapping_rate(%),mean_depth,mean_dedup_depth,dup_rate(%),average_insert_size,std_insert_size" \
> $qc_dir/QC_summary.csv
fi

# ---------------------------------------------------------------------- #
# ---------------------------  Pipeline  ------------------------------- #
# ---------------------------------------------------------------------- #

export SENTIEON_LICENSE=${sentieon_license};

for ifile in $input_folder/*_tumor_R1.fastq.gz;
do 
    sampleID=`basename ${ifile%%"_tumor"*}`

    # step1 - trim reads
    if [[  $do_trim == "on"  ]]; 
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- trimming reads";

        $fastp --in1 $input_folder/${sampleID}_tumor_R1.fastq.gz \
        --in2 $input_folder/${sampleID}_tumor_R2.fastq.gz \
        --out1 $trim_dir/${sampleID}_trim_R1.tumor.fastq.gz \
        --out2 $trim_dir/${sampleID}_trim_R2.tumor.fastq.gz \
        -c --length_required 3 --detect_adapter_for_pe -p \
        --thread ${thread} \
        --html $trim_dir/${sampleID}.tumor.trim.html \
        --json $trim_dir/${sampleID}.tumor.trim.json;

        if [[  $mode == "pair"  ]];
        then
            $fastp --in1 $input_folder/${sampleID}_normal_R1.fastq.gz \
            --in2 $input_folder/${sampleID}_normal_R2.fastq.gz \
            --out1 $trim_dir/${sampleID}_trim_R1.normal.fastq.gz \
            --out2 $trim_dir/${sampleID}_trim_R2.normal.fastq.gz \
            -c --length_required 3 --detect_adapter_for_pe -p \
            --thread ${thread} \
            --html $trim_dir/${sampleID}.normal.trim.html \
            --json $trim_dir/${sampleID}.normal.trim.json;
        fi
    fi

    # step2 - align & sort
    if [[  $do_align == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- alignment & sorting";

        ($sentieon bwa mem -M -R "@RG\tID:${sampleID}.tumor\tSM:${sampleID}.tumor\tPL:illumina" \
        -t ${thread} -K 10000000 ${ref} \
        $trim_dir/${sampleID}_trim_R1.tumor.fastq.gz \
        $trim_dir/${sampleID}_trim_R2.tumor.fastq.gz \
        || echo -n 'error' ) \
        | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.tumor.sorted.bam \
        -t ${thread} --sam2bam -i -;

        if [[  $mode == "pair"  ]];
        then
            ($sentieon bwa mem -M -R "@RG\tID:${sampleID}.normal\tSM:${sampleID}.normal\tPL:illumina" \
            -t ${thread} -K 10000000 ${ref} \
            $trim_dir/${sampleID}_trim_R1.normal.fastq.gz \
            $trim_dir/${sampleID}_trim_R2.normal.fastq.gz \
            || echo -n 'error' ) \
            | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.normal.sorted.bam \
            -t ${thread} --sam2bam -i -;
        fi

    fi

    # step3 - mark duplicates
    if [[  $do_dedup == "on"  ]] && [[ $dedup == true ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- deduplication";

        ${sentieon} driver -t ${thread} \
        -i ${align_dir}/${sampleID}.tumor.sorted.bam \
        --algo LocusCollector \
        --fun score_info ${align_dir}/${sampleID}.tumor.score.txt;

        ${sentieon} driver -t ${thread} \
        -i ${align_dir}/${sampleID}.tumor.sorted.bam \
        --algo Dedup --score_info ${align_dir}/${sampleID}.tumor.score.txt \
        --metrics ${align_dir}/${sampleID}.dedup_metrics.txt \
        ${align_dir}/${sampleID}.tumor.sorted.dedup.bam;


        if [[  $mode == "pair"  ]]; then
            ${sentieon} driver -t ${thread} \
            -i ${align_dir}/${sampleID}.normal.sorted.bam \
            --algo LocusCollector \
            --fun score_info ${align_dir}/${sampleID}.normal.score.txt;

            ${sentieon} driver -t ${thread} \
            -i ${align_dir}/${sampleID}.normal.sorted.bam \
            --algo Dedup \
            --score_info ${align_dir}/${sampleID}.normal.score.txt \
            --metrics ${align_dir}/${sampleID}.normal.dedup_metrics.txt \
            ${align_dir}/${sampleID}.normal.sorted.dedup.bam;
        fi
    fi

    # determine bam
    if [[ $dedup == true ]]; then
        normal_bam=${align_dir}/${sampleID}.normal.sorted.dedup.bam
        tumor_bam=${align_dir}/${sampleID}.tumor.sorted.dedup.bam
    else
        normal_bam=${align_dir}/${sampleID}.normal.sorted.bam
        tumor_bam=${align_dir}/${sampleID}.tumor.sorted.bam
    fi

    # step4 - quality control
    if [[  $do_qc == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- qc";

        $samtools stats -@ ${thread} $tumor_bam > ${qc_dir}/${sampleID}.tumor.stats.txt;

        tumor_r1=$(du $input_folder/${sampleID}_tumor_R1.fastq.gz -shL |awk '{print $1}');
        tumor_r2=$(du $input_folder/${sampleID}_tumor_R2.fastq.gz -shL |awk '{print $1}');

        tumor_raw_reads=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['total_reads'])"`

        tumor_clean_reads=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['total_reads'])"`

        tumor_raw_bases=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['total_bases'])"`

        tumor_clean_bases=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['total_bases'])"`

        tumor_qc_rate=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['q30_rate'])"`

        tumor_map_reads=`awk -F"\t" '$2 == "reads mapped:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt`;
        tumor_mapping_rate=`python3 -c "print('{:.2f}%'.format($tumor_map_reads/$tumor_clean_reads*100))"`
        tumor_map_bases=`awk -F"\t" '$2 == "bases mapped (cigar):" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt`;
        tumor_mean_depth=`python3 -c "print('{:.2f}'.format(${tumor_map_bases}/3101788170))"`;
        tumor_dup_bases=`awk -F"\t" '$2 == "bases duplicated:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt`;
        tumor_dedup_depth=`python3 -c "print('{:.2f}'.format((${tumor_map_bases}-${tumor_dup_bases})/3101788170))"`;
        tumor_dup_rate=`python3 -c "print('{:.2f}'.format(100-${tumor_dedup_depth}/${tumor_mean_depth}*100))"`;

        tumor_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt);
        tumor_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt);
        
        echo "${sampleID}.tumor,${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},\
        ${tumor_qc_rate},${tumor_mapping_rate},${tumor_mean_depth},${tumor_dedup_depth},${tumor_dup_rate},${tumor_insert_size},${tumor_insert_std}" \
        >> ${qc_dir}/QC_summary.csv

        if [[  $mode == "pair"  ]];
        then
        $samtools stats -@ ${thread} $normal_bam > ${qc_dir}/${sampleID}.normal.stats.txt;

        normal_r1=$(du $input_folder/${sampleID}_normal_R1.fastq.gz -shL |awk '{print $1}');
        normal_r2=$(du $input_folder/${sampleID}_normal_R2.fastq.gz -shL |awk '{print $1}');

        normal_raw_reads=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['total_reads'])"`

        normal_clean_reads=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['total_reads'])"`

        normal_raw_bases=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['total_bases'])"`

        normal_clean_bases=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
        print(fh['summary']['after_filtering']['total_bases'])"`

        normal_qc_rate=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['q30_rate'])"`

        normal_map_reads=`awk -F"\t" '$2 == "reads mapped:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt`;
        normal_mapping_rate=`python3 -c "print('{:.2f}%'.format($normal_map_reads/$normal_clean_reads*100))"`
        normal_map_bases=`awk -F"\t" '$2 == "bases mapped (cigar):" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt`;
        normal_dup_bases=`awk -F"\t" '$2 == "bases duplicated:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt`;
        normal_mean_depth=`python3 -c "print('{:.2f}'.format(${normal_map_bases}/3101788170))"`;
        normal_dedup_depth=`python3 -c "print('{:.2f}'.format((${normal_map_bases}-${normal_dup_bases})/3101788170))"`;
        normal_dup_rate=`python3 -c "print('{:.2f}'.format(100-${normal_dedup_depth}/${normal_mean_depth}*100))"`;

        normal_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt);
        normal_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt);

        echo "${sampleID}.normal,${normal_r1}/${normal_r2},${normal_raw_reads},${normal_raw_bases},${normal_clean_reads},${normal_clean_bases},\
        ${normal_qc_rate},${normal_mapping_rate},${normal_mean_depth},${normal_dedup_depth},${normal_dup_rate},${normal_insert_size},${normal_insert_std}" \
        >> ${qc_dir}/QC_summary.csv
        fi
    fi  
    # step4 - read count
    if [[  $do_readcount == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- count reads";

        $readCounter --window 1000000 --quality 20 \
        $tumor_bam > ${cnv_dir}/${sampleID}.tumor.wig;
        
        if [[  $mode == "pair"  ]];
        then
            $readCounter --window 1000000 --quality 20 \
            $normal_bam > ${cnv_dir}/${sampleID}.normal.wig;
        fi
    fi
    
    # step5 - run ichorCNA
    if [[  $do_icorcna == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- run ichorCNA";

        if [[ ! -d $cnv_dir/${sampleID}_ichorCNA ]]; then
            mkdir $cnv_dir/${sampleID}_ichorCNA
        fi

        if [[  $mode == "pair"  ]];
        then
            $Rscript ${ichorCNA}/scripts/runIchorCNA.R --id ${sampleID} \
            --WIG ${cnv_dir}/${sampleID}.tumor.wig \
            --NORMWIG ${cnv_dir}/${sampleID}.normal.wig \
            --ploidy $ploidy \
            --normal $nfx \
            --maxCN $maxCN \
            --gcWig ${ichorCNA}/inst/extdata/gc_hg19_1000kb.wig \
            --mapWig ${ichorCNA}/inst/extdata/map_hg19_1000kb.wig \
            --centromere ${ichorCNA}/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
            --scStates $subclone \
            --estimateScPrevalence $scpre \
            --txnStrength 10000 \
            --chrs $chrs \
            --chrTrain $chrs \
            --outDir $cnv_dir/${sampleID}_ichorCNA \
            --libdir $ichorCNA;
        fi

        if [[  $mode == "single"  ]];
        then
            $Rscript ${ichorCNA}/scripts/runIchorCNA.R --id ${sampleID} \
            --WIG ${cnv_dir}/${sampleID}.tumor.wig \
            --ploidy $ploidy \
            --normal $nfx \
            --maxCN $maxCN \
            --gcWig ${ichorCNA}/inst/extdata/gc_hg19_1000kb.wig \
            --mapWig ${ichorCNA}/inst/extdata/map_hg19_1000kb.wig \
            --centromere ${ichorCNA}/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
            --scStates $subclone \
            --estimateScPrevalence $scpre \
            --txnStrength 10000 \
            --chrs $chrs \
            --chrTrain $chrs \
            --outDir $cnv_dir/${sampleID}_ichorCNA \
            --libdir $ichorCNA;
        fi

        if [[  $mode == "pon"  ]];
        then
            $Rscript ${ichorCNA}/scripts/runIchorCNA.R --id ${sampleID} \
            --WIG ${cnv_dir}/${sampleID}.tumor.wig \
            --normalPanel $pon \
            --ploidy $ploidy \
            --normal $nfx \
            --maxCN $maxCN \
            --gcWig ${ichorCNA}/inst/extdata/gc_hg19_1000kb.wig \
            --mapWig ${ichorCNA}/inst/extdata/map_hg19_1000kb.wig \
            --centromere ${ichorCNA}/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
            --scStates $subclone \
            --estimateScPrevalence $scpre \
            --txnStrength 10000 \
            --chrs $chrs \
            --chrTrain $chrs \
            --outDir $cnv_dir/${sampleID}_ichorCNA \
            --libdir $ichorCNA;
        fi
    fi

    # step6 - calculate CIN score
    if [[  $do_cin == "on"  ]];
    then
        $Rscript $compute_cin ${cnv_dir}/${sampleID}_ichorCNA/MRD210430-079.RData
    fi

done


echo "LOGGING: `date --rfc-3339=seconds` -- Analysis finished";
