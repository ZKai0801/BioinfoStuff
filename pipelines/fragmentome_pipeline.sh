#!/usr/bin/bash

version="v2.0"

# ----------------------------------- Description ------------------------------------- #
# Perform:                                                                              #
#   Trimming (fastp) +                                                                  #
#   Alignment (bowtie2) +                                                               #
#   Deduplication (picard) +                                                            #
#   Quality Control (fastqc/samtools/multiqc) +                                         #
#   calculate insert size for each pair of reads  +                                     #
#   insert size re-weighting                                                            #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash fragmentome_pipeline.sh [input_folder] [output_folder]             #
#                                                                                       #
# input_folder should contain pair-end sequenced Fastq files.                           #
# Each sample is expected to have two Fastq files, with naming convention of:           #
# ${sampleID}_R[1|2].fastq.gz                                                           #
# ------------------------------------------------------------------------------------- #


# --------------------------- set parameters --------------------------- #
# ---------------------------------------------------------------------- #

sentieon_license="172.16.11.242:8991"
thread=8

fastp="/data/ngs/softs/fastp/fastp"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-202010/bin/sentieon"
samtools="/public/software/samtools-1.9/samtools"
bedtools="/public/software/bedtools2/bin/bedtools"
bowtie2="/public/home/kai/softwares/bowtie2-2.4.4/bowtie2"
picard="/public/software/picard.jar"
fastqc="/public/software/FastQC/fastqc"
python3_multiqc='/public/ngs/softs/MultiQC-master/python3/python3.6.4/bin/python3.6'
Rscript="/public/software/R_3.5.1/bin/Rscript"

ref="/public/database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
bowtie_index="/public/database/GATK_Resource_Bundle/hg19/hg19_bowtie2"
filt_ref="/public/home/kai/projects/fragmentome/sources/hg19.filters.bed"
ref_gc="/public/home/kai/projects/fragmentome/sources/target20.tsv"
templ_bins="/public/home/kai/projects/fragmentome/sources/templ_bins.tsv"


# switch (on||off)
do_trim="on"
do_align="on"
do_dedup="on"
do_qc="on"
do_frag="on"
adjust_gc="on"

# ------------------------------ argparser ----------------------------- #
# ---------------------------------------------------------------------- #

if [[  $1 == '-h'  ]]; then
    echo "Usage: ./fragmentome_pipeline.sh [input_folder] [output_folder]"
    echo "-------------------------------------------------------------------------"
    echo "[input_folder] should contain fastq files with following naming system:"
    echo "  \${sampleID}_R[1|2].fastq.gz"
    exit 0
fi

input_folder=`realpath $1`
output_folder=`realpath $2`

if [[ ! -d $input_folder  ]];
then
    echo "Error: input_folder does not Found!"
    exit 1
fi

if [[ ! -d $output_folder ]]; 
then
    mkdir -p $output_folder
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

frag_dir=$output_folder/fragmentome/;
if [[ ! -d $frag_dir ]]; then
    mkdir -p $frag_dir/tmp/
fi

# ---------------------------  LOGGING  -------------------------------- #
# ---------------------------------------------------------------------- #

echo "LOGGING: `date --rfc-3339=seconds` -- Analysis started"
echo "LOGGING: This is the fragmentome_pipeline.sh pipeline"
echo "========================================================"
echo "LOGGING: -- settings -- input folder -- ${input_folder}"
echo "LOGGING: -- settings -- output folder -- ${output_folder}"
echo "LOGGING: -- settings -- version -- ${version}"
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

for ifile in $input_folder/*_R1.fastq.gz;
do 
    sampleID=`basename ${ifile%%"_R1"*}`

    # step1 - trim reads
    # -Q indicates no quality filtering will be applied
    if [[  $do_trim == "on"  ]]; 
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- trimming";
        
        $fastp -i $input_folder/${sampleID}_R1.fastq.gz \
        -I $input_folder/${sampleID}_R2.fastq.gz \
        -o $trim_dir/${sampleID}_trim_R1.fastq.gz \
        -O $trim_dir/${sampleID}_trim_R2.fastq.gz \
        -Q --detect_adapter_for_pe \
        --html $trim_dir/${sampleID}.trim.html \
        --json $trim_dir/${sampleID}.trim.json \
        --thread $thread;
    fi

    # step2 - align & sort
    if [[  $do_align == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- alignment & sorting";

        $bowtie2 --seed 42 --very-fast --end-to-end --threads $thread \
        -x $bowtie_index -1 $trim_dir/${sampleID}_trim_R1.fastq.gz \
        -2 $trim_dir/${sampleID}_trim_R2.fastq.gz | \
        $samtools sort -@ $thread - | \
        $samtools view -1bS -@ $thread > ${align_dir}/${sampleID}.sorted.bam;

        $samtools index ${align_dir}/${sampleID}.sorted.bam;

        rm $trim_dir/${sampleID}_trim_R*.fastq.gz

        # ($sentieon bwa mem -M -R "@RG\tID:${sampleID}\tSM:${sampleID}\tPL:illumina" \
        # -t ${thread} -K 10000000 ${ref} \
        # $trim_dir/${sampleID}_trim_R1.fastq.gz \
        # $trim_dir/${sampleID}_trim_R2.fastq.gz \
        # || echo -n 'error' ) \
        # | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.sorted.bam \
        # -t ${thread} --sam2bam -i -;
    fi

    # step3 - deduplicates
    if [[  $do_dedup == "on"  ]]; 
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- deduplication";

        java -jar $picard MarkDuplicates \
        I=${align_dir}/${sampleID}.sorted.bam \
        O=${align_dir}/${sampleID}.sorted.dedup.bam \
        M=${align_dir}/${sampleID}.dup_metrics.txt; 

        $samtools index ${align_dir}/${sampleID}.sorted.dedup.bam

        rm ${align_dir}/${sampleID}.sorted.bam*

        # ${sentieon} driver -t ${thread} \
        # -i ${align_dir}/${sampleID}.sorted.bam \
        # --algo LocusCollector \
        # --fun score_info ${align_dir}/${sampleID}.score.txt;

        # ${sentieon} driver -t ${thread} \
        # -i ${align_dir}/${sampleID}.sorted.bam \
        # --algo Dedup --score_info ${align_dir}/${sampleID}.score.txt \
        # --metrics ${align_dir}/${sampleID}.dedup_metrics.txt \
        # ${align_dir}/${sampleID}.sorted.dedup.bam;
    fi

    # step4 - qc
    if [[  $do_qc == "on"  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- qc";

        $fastqc $input_folder/${sampleID}_R1.fastq.gz \
        $input_folder/${sampleID}_R2.fastq.gz \
        $trim_dir/${sampleID}_trim_R1.fastq.gz \
        $trim_dir/${sampleID}_trim_R2.fastq.gz \
        -t ${thread} -o ${qc_dir};

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

        tumor_map_reads=`awk -F"\t" '$2 == "reads mapped:" {print $3}' ${qc_dir}/${sampleID}.stats.txt`;
        tumor_mapping_rate=`python3 -c "print('{:.2f}%'.format($tumor_map_reads/$tumor_clean_reads*100))"`
        tumor_map_bases=`awk -F"\t" '$2 == "bases mapped (cigar):" {print $3}' ${qc_dir}/${sampleID}.stats.txt`;
        tumor_mean_depth=`python3 -c "print('{:.2f}'.format(${tumor_map_bases}/3101788170))"`;
        tumor_dup_bases=`awk -F"\t" '$2 == "bases duplicated:" {print $3}' ${qc_dir}/${sampleID}.stats.txt`;
        tumor_dedup_depth=`python3 -c "print('{:.2f}'.format((${tumor_map_bases}-${tumor_dup_bases})/3101788170))"`;
        tumor_dup_rate=`python3 -c "print('{:.2f}'.format(100-${tumor_dedup_depth}/${tumor_mean_depth}*100))"`;

        tumor_insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);
        tumor_insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);
        
        echo "${sampleID},${tumor_r1}/${tumor_r2},${tumor_raw_reads},${tumor_raw_bases},${tumor_clean_reads},${tumor_clean_bases},\
        ${tumor_qc_rate},${tumor_mapping_rate},${tumor_mean_depth},${tumor_dedup_depth},${tumor_dup_rate},${tumor_insert_size},${tumor_insert_std}" \
        >> ${qc_dir}/QC_summary.csv
    fi

    # step5 - calc fragment lengths
    if [[  $do_frag == "on"  ]];
    then
        echo "LOGGING: $sampleID -- `date --rfc-3339=seconds` -- calc fragment lengths "

        $samtools sort -@ $thread -n ${align_dir}/${sampleID}.sorted.dedup.bam | \
        $samtools fixmate -@ $thread -r - - | \
        $bedtools bamtobed -bedpe -i - | \
        awk '($2 != -1) && ($5 != -1) && ($1 == $4)' | \
        cut -f 1,2,6,8 | \
        awk '($4 >= 30) && ($3 - $2 < 1000)' | \
        sort --parallel 8 -k1,1 -k2,2n -S 32G -T $frag_dir/tmp/ | \
        grep -v "_" | grep -v "chrM" > \
        $frag_dir/$sampleID.bed;
    fi

    # step5 - adjust for gc content
    if [[  $adjust_gc == "on"  ]];
    then
        echo "LOGGING: $sampleID -- `date --rfc-3339=seconds` -- adjust for gc content"
        # get gc content for each fragment
        $bedtools nuc -fi $ref -bed $frag_dir/$sampleID.bed | \
        awk 'NR!=1 {OFS="\t"; print $1, $2, $3, $4, $6}' | \
        awk '(length($1) <= 5) ($3 - $2 >= 100) && ($3 - $2 <= 220)' > $frag_dir/$sampleID.step1_gc.bed;

        # filter blacklisted regions
        $bedtools intersect -v -a $frag_dir/$sampleID.step1_gc.bed -b $filt_ref | \
        awk '{printf "%s\t%s\t%s\t%s\t%.2f\n", $1, $2, $3, $4, $5 }' | \
        > $frag_dir/$sampleID.step2_filt.bed

        Rscript $binning_frag -i $frag_dir/$sampleID.step2_filt.bed -r $ref_gc -t $templ_bins -o $frag_dir/$sampleID.5mb.csv
    fi
done


if [[  $do_qc == "on"  ]];
then
    $python3_multiqc -m multiqc $qc_dir/*zip -o $qc_dir;
fi


echo "LOGGING: $sampleID -- `date --rfc-3339=seconds` -- analysis completed"



