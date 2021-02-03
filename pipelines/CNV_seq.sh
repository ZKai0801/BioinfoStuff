#!/usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
# Perform:                                                                              #
#   Trimming (fastp) +                                                                  #
#   Alignment (Sentieon-bwa) +                                                          #
#   statistic analysis +                                                                #
#   CNV calling (ichorCNA)                                                              #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash CNV_seq.sh [input_folder] [output_folder]                          #
#                                                                                       #
# input_folder should contain pairs of matched tumor/normal pair-end sequenced Fastqs   #
#                                                                                       #
# Naming convention:                                                                    #
# tumor:   ${sampleID}_tumor_R[1|2].fastq.gz                                            #
# normal:  ${sampleID}_normal_R[1|2].fastq.gz                                           #
#                                                                                       #
# ------------------------------------------------------------------------------------- #


# --------------------------- set parameters --------------------------- #
# ---------------------------------------------------------------------- #
sentieon_license="192.168.1.186:8990"
thread=8

# software
fastp="/data/ngs/softs/fastp/fastp"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-201808.08/bin/sentieon"
samtools="/public/software/samtools-1.9/samtools"
bamdst="/public/software/bamdst/bamdst"
ichorCNA="/public/software/ichorCNA/"
readCounter="/public/software/hmmcopy_utils/bin/readCounter"
Rscript="/public/software/R_3.6.3/bin/Rscript"

ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"

# switch (on||off)
do_trim=""
do_align=""
do_qc=""
do_readcount=""
do_icorcna="on"

# ------------------------------ argparser ----------------------------- #
# ---------------------------------------------------------------------- #
if [[  $1 == '-h'  ]]; then
    echo "Usage: ./CNV_seq.sh [input_folder] [output_folder]"
    echo "-------------------------------------------------------------------------"
    echo "[input_folder] should contain fastq files with following naming system:"
    echo "\${sampleID}_[tumor|normal]_R[1|2].fastq.gz"
    exit 0
fi

input_folder=`realpath $1`
output_folder=`realpath $2`

if [[ ! -d $input_folder  ]];
then
    echo "Error: input_folder does not Found!"
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
echo "========================================================"

if [[  $do_qc == "on"  ]];
then
echo "sampleID,fastq_size,raw_reads,raw_bases,clean_reads,clean_bases,qc30_rate,\
mapping_rate(%),mean_depth,average_insert_size,std_insert_size" \
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

        $fastp --in1 $input_folder/${sampleID}_normal_R1.fastq.gz \
        --in2 $input_folder/${sampleID}_normal_R2.fastq.gz \
        --out1 $trim_dir/${sampleID}_trim_R1.normal.fastq.gz \
        --out2 $trim_dir/${sampleID}_trim_R2.normal.fastq.gz \
        -c --length_required 3 --detect_adapter_for_pe -p \
        --thread ${thread} \
        --html $trim_dir/${sampleID}.normal.trim.html \
        --json $trim_dir/${sampleID}.normal.trim.json;
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

        ($sentieon bwa mem -M -R "@RG\tID:${sampleID}.normal\tSM:${sampleID}.normal\tPL:illumina" \
        -t ${thread} -K 10000000 ${ref} \
        $trim_dir/${sampleID}_trim_R1.normal.fastq.gz \
        $trim_dir/${sampleID}_trim_R2.normal.fastq.gz \
        || echo -n 'error' ) \
        | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.normal.sorted.bam \
        -t ${thread} --sam2bam -i -;

    fi

    # step3 - quality control
    if [[  $do_qc == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- qc";

        $samtools stats -@ ${thread} ${align_dir}/${sampleID}.normal.sorted.bam > ${qc_dir}/${sampleID}.normal.stats.txt;
        $samtools stats -@ ${thread} ${align_dir}/${sampleID}.tumor.sorted.bam > ${qc_dir}/${sampleID}.tumor.stats.txt;

        tumor_r1=$(du $input_folder/${sampleID}_tumor_R1.fastq.gz -shL |awk '{print $1}');
        tumor_r2=$(du $input_folder/${sampleID}_tumor_R2.fastq.gz -shL |awk '{print $1}');
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

        normal_qc_rate=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.normal.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['q30_rate'])"`

        tumor_qc_rate=`python3 -c "import json; \
        fh = json.load(open('$trim_dir/${sampleID}.tumor.trim.json', 'r')); \
        print(fh['summary']['before_filtering']['q30_rate'])"`

        normal_map_reads=`awk -F"\t" '$2 == "reads mapped:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt`;
        normal_mapping_rate=`python3 -c "print('{:.2f}%'.format($normal_map_reads/$normal_clean_reads*100))"`

        tumor_map_reads=`awk -F"\t" '$2 == "reads mapped:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt`;
        tumor_mapping_rate=`python3 -c "print('{:.2f}%'.format($tumor_map_reads/$tumor_clean_reads*100))"`

        normal_map_bases=`awk -F"\t" '$2 == "bases mapped (cigar):" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt`;
        tumor_map_bases=`awk -F"\t" '$2 == "bases mapped (cigar):" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt`;

        normal_mean_depth=`python3 -c "print('{:.2f}'.format(${normal_map_bases}/3101788170))"`;
        tumor_mean_depth=`python3 -c "print('{:.2f}'.format(${tumor_map_bases}/3101788170))"`;

        normal_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt);
        normal_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.normal.stats.txt);

        tumor_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt);
        tumor_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.tumor.stats.txt);
        
        echo "${sampleID}.normal,${normal_r1}/${normal_r2},${normal_raw_reads},${normal_raw_bases},${normal_clean_reads},${normal_clean_bases},\
        ${normal_qc_rate},${normal_mapping_rate},${normal_mean_depth},${normal_insert_size},${normal_insert_std}" \
        >> ${qc_dir}/QC_summary.csv 

        echo "${sampleID}.tumor,${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},\
        ${tumor_qc_rate},${tumor_mapping_rate},${tumor_mean_depth},${tumor_insert_size},${tumor_insert_std}" \
        >> ${qc_dir}/QC_summary.csv
    fi

    # step4 - read count
    if [[  $do_readcount == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- count reads";

        $readCounter --window 1000000 --quality 20 \
        ${align_dir}/${sampleID}.tumor.sorted.bam > \
        ${cnv_dir}/${sampleID}.tumor.wig;
        
        $readCounter --window 1000000 --quality 20 \
        ${align_dir}/${sampleID}.normal.sorted.bam > \
        ${cnv_dir}/${sampleID}.normal.wig;
    fi
    
    # step5 - run ichorCNA
    if [[  $do_icorcna == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- run ichorCNA";

        if [[ ! -d $cnv_dir/${sampleID}_ichorCNA ]]; then
            mkdir $cnv_dir/${sampleID}_ichorCNA
        fi

        $Rscript ${ichorCNA}/scripts/runIchorCNA.R --id ${sampleID} \
        --WIG ${cnv_dir}/${sampleID}.tumor.wig \
        --ploidy "c(2,3)" \
        --normal "c(0.5,0.6,0.7,0.8,0.9)" \
        --maxCN 5 \
        --gcWig ${ichorCNA}/inst/extdata/gc_hg19_1000kb.wig \
        --mapWig ${ichorCNA}/inst/extdata/map_hg19_1000kb.wig \
        --centromere ${ichorCNA}/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
        --scStates "c(1,3)" \
        --txnStrength 10000 \
        --outDir $cnv_dir/${sampleID}_ichorCNA;

    fi

done


echo "LOGGING: `date --rfc-3339=seconds` -- Analysis finished";
