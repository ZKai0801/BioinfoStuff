#!/usr/bin/bash

# ----------------------------------- Description ------------------------------------- #
# Perform:                                                                              #
#   Trimming (Cutadapt) +                                                               #
#   Alignment (Star) +                                                                  #
#   Fusion detection (Star-fusion)                                                      #
# ------------------------------------------------------------------------------------- #
# Usage:                                                                                #
#   [admin@kai]$bash RNA_fusion.sh [input_folder] [output_folder]                       #
#                                                                                       #
# input_folder should contain pairs of pair-end sequenced Fastqs                        #
#                                                                                       #
# Naming convention:                                                                    #
#    ${sampleID}_R[1|2].fastq.gz                                                        #
#                                                                                       #
# ------------------------------------------------------------------------------------- #


# --------------------------- set parameters --------------------------- #
# ---------------------------------------------------------------------- #
thread=8

# software
cutadapt="~/.local/bin/cutadapt"
star=""
starfusion="/public/software/STAR-Fusion-v1.8.0_FULL/STAR-Fusion"
samtools="/public/software/samtools-1.9/samtools"


# databases
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
gtf=""

# switch (on||off)
do_trim="on"
do_align="on"
do_qc="on"


# ------------------------------ argparser ----------------------------- #
# ---------------------------------------------------------------------- #
if [[  $1 == '-h'  ]]; then
    echo "Usage: ./RNA_fusion.sh [input_folder] [output_folder]"
    echo "-------------------------------------------------------------------------"
    echo "[input_folder] should contain fastq files with following naming system:"
    echo "\${sampleID}_R[1|2].fastq.gz"
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

fusion_dir=$output_folder/fusion/;
if [[ ! -d $fusion_dir ]]; then
    mkdir $fusion_dir
fi

qc_dir=$output_folder/qc/;
if [[ ! -d $qc_dir ]]; then
    mkdir $qc_dir
fi

# ---------------------------  LOGGING  -------------------------------- #
# ---------------------------------------------------------------------- #
echo "LOGGING: `date --rfc-3339=seconds` -- Analysis started"
echo "========================================================"
echo "LOGGING: -- settings -- input folder -- ${input_folder}"
echo "LOGGING: -- settings -- output folder -- ${output_folder}"
echo "========================================================"



# ---------------------------------------------------------------------- #
# ---------------------------  Pipeline  ------------------------------- #
# ---------------------------------------------------------------------- #


for ifile in $input_folder/*_R1.fastq.gz;
do 
    sampleID=`basename ${ifile%%"_R1"*}`

    if [[  $do_trim == "on"  ]];
    then
        $cutadapt --cores $thread \
        -a "A{10}" 
    fi
        echo
        echo ============================================
        echo "`date` --- start running STAR + samtools"
        if [[ ! -d `dirname $2`/star_index ]]; then
            echo 'build reference index for STAR'
            mkdir `dirname $2`/star_index
            STAR --runMode genomeGenerate \
            --genomeFastaFiles $2 \
            --sjdbGTFfile $3 \
            --genomeDir `dirname $2` \
            --runThreadN 16
        fi
        STAR --genomeDir `dirname $2`/star_index \
        --runThreadN 16 \
        --readFilesIn ${trimout1::-3} ${trimout2::-3} \
        --outFileNamePrefix $aligndir/`basename ${read1/R1.fastq.gz/}` \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMattributes Standard 

        samtools index $aligndir/`basename ${read1/R1.fastq.gz/Aligned.sortedByCoord.out.bam}`
        echo ============================================
    else
        echo "cannot find the paired fastq file for $read1 ！"
    fi
done
