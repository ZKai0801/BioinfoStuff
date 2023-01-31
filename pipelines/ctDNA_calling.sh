#! /usr/bin/bash

#PBS -l nodes=1:ppn=8

version="v1.3"

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
# umi_templ="4M+T,4M+T"   # capp-seq 原版
# umi_templ="6M1S+T,6M1S+T" # 生产流程默认UMI
umi_templ="6M+T,6M+T"

# path to software & scripts
fastp="/data/ngs/softs/fastp/fastp"
sentieon="/data/ngs/softs/sentieon/sentieon-genomics-202112.06/bin/sentieon"
bcftools="/data/ngs/softs/bcftools/bcftools"
samtools="/data/ngs/softs/samtools-1.16.1/samtools"
bamdst="/data/ngs/softs/bamdst/bamdst"
vep="/data/ngs/softs/98vep/ensembl-vep/vep"
bgzip="/data/ngs/softs/htslib/bgzip"
tabix="/data/ngs/softs/htslib/tabix"
varscan="/data/ngs/softs/varscan/VarScan.v2.4.2.jar"
genefuse="/data/ngs/softs/genefuse/genefuse"

anno_hgvs="/data/ngs/scripts/workflow/mrd/anno_hgvs_mrd.py"
hotspot_filt="/data/ngs/scripts/workflow/script/hotspot_filter.py"

sentieon_license="192.168.1.186:8990"
thread=6

# additional files
ref="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
dbsnp="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
k1="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/dbsnp_138.hg19.vcf.gz"
k2="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz"
k3="/data/ngs/database/soft_database/GATK_Resource_Bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz"
refflat="/data/ngs/database/soft_database/GATK_Resource_Bundle/refFlat.txt"
vep_dir="/data/ngs/database/soft_database/vep_98"
cache_version="98"
clinic_transcripts="/data/ngs/database/publicDataBase/clinic_transcript.tsv"
cbioportal="/data/ngs/database/publicDataBase/cBioportal_hotspot.csv"
fusion_templ="/data/ngs/database/publicDataBase/GRD_fusions.tsv"


# switch (on||off)
do_trim="on"
do_align="on"
do_qc="on"
do_realign="on"
do_tnscope="on"
# do_mpileup="on"
# do_varscan="on"
do_filt="on"
do_anno="on"
do_filt2="on"
do_fusion="on"


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
# output_folder="/data/ngs/files/clean/YDLYX-20230117-L-01-2023-01-191316/$(basename $input_folder)/"

if [[  -f $3  ]];
then
   bed=`realpath $3`
else
   bed="empty"
fi

# bed="/data/ngs/database/bed/GRD_v2_probeCov.bed"

if [[ ! -d $input_folder  ]];
then
    echo "IOError: $input_folder does not Found!"
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

snv_dir=$output_folder/snv/;
tnscope_dir=$snv_dir/tnscope/;
varscan_dir=$snv_dir/varscan/;
if [[  ! -d $snv_dir  ]];
then
    mkdir $snv_dir
    mkdir -p $tnscope_dir
    mkdir -p $varscan_dir
fi

qc_dir=$output_folder/qc/;
if [[ ! -d $qc_dir ]]; then
    mkdir $qc_dir
fi

fusion_dir=$output_folder/fusion/
if [[  ! -d $fusion_dir  ]]; then
    mkdir $fusion_dir
fi

# ---------------------------  LOGGING  -------------------------------- #
# ---------------------------------------------------------------------- #
echo "LOGGING: `date --rfc-3339=seconds` -- Analysis started"
echo "LOGGING: This is the ctDNA_calling.sh pipeline"
echo "========================================================"
echo "LOGGING: -- settings -- input folder -- ${input_folder}"
echo "LOGGING: -- settings -- output folder -- ${output_folder}"
echo "LOGGING: -- settings -- BED file -- ${bed}"
echo "LOGGING: -- settings -- mode -- ${mode}"
echo "========================================================"


echo "sampleID,fastq_size,raw_reads,raw_bases,clean_reads,clean_bases,\
qc30(%),trim_percent(%),mapping_rate(%),on-target_percent(%),\
mean_depth,mean_dedup_depth,dup_rate(%),\
average_insert_size,std_insert_size,\
Uniformity_0.1X(%),Uniformity_0.2X(%),\
Uniformity_0.5X(%),Uniformity_1X(%),\
300x_depth_percent(%),500x_depth_percent(%),\
1000x_depth_percent(%),1500x_depth_percent(%),\
2000x_depth_percent(%)" \
> $qc_dir/QC_summary.csv



# ---------------------------------------------------------------------- #
# -----------------------------  Steps  -------------------------------- #
# ---------------------------------------------------------------------- #

# 1. trimming
function run_trim {
    # remove reads with length < 50bp;
    # remove reads with >30% of bases below Q30;
    # remove low-complexity reads;
    # remove polyG tail; polyX tails;
    # remove read pairs that contain mismatches in overlapped regions
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
}


# 2. alignment 
function run_alignment {
    if [[ ! -n $umi_templ  ]];
    then
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- alignment & sorting";

        ($sentieon bwa mem -M -R "@RG\tID:${sampleID}\tSM:${sampleID}\tPL:illumina" \
        -t ${thread} -K 10000000 ${ref} \
        ${trim_dir}/${sampleID}_trim_R1.fastq.gz \
        ${trim_dir}/${sampleID}_trim_R2.fastq.gz \
        || echo -n 'error' ) \
        | ${sentieon} util sort -r ${ref} -o ${align_dir}/${sampleID}.sorted.bam \
        -t ${thread} --sam2bam -i -;

    else
        echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- umi-consensus & alignment & sorting";

        $sentieon umi extract -d $umi_templ $trim_dir/${sampleID}_trim_R1.fastq.gz $trim_dir/${sampleID}_trim_R2.fastq.gz | \
        ($sentieon bwa mem -t ${thread} -p -C -M -R "@RG\\tID:${sampleID}\\tSM:${sampleID}\\tPL:ILLUMINA\\tLB:lib" -K 10000000 ${ref} - || echo -n 'error') | \
        $sentieon util sort -t ${thread} --sam2bam -i - -o ${align_dir}/${sampleID}.sorted.bam;
        
        $sentieon umi consensus -i ${align_dir}/${sampleID}.sorted.bam --input_format BAM -o ${trim_dir}/${sampleID}.umi_consensus.fq.gz;
        
        ($sentieon bwa mem -t ${thread} -p -C -M -R "@RG\\tID:${sampleID}\\tSM:${sampleID}\\tPL:ILLUMINA\\tLB:lib" -K 10000000 ${ref} ${trim_dir}/${sampleID}.umi_consensus.fq.gz || echo -n 'error' ) | \
        $sentieon util sort --umi_post_process --sam2bam -i - -o ${align_dir}/${sampleID}.sorted.dedup.bam

    fi
}


# 3. qc
function run_qc {
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- qc";
    
    bam=${align_dir}/${sampleID}.sorted.bam

    if [[ ! -f ${qc_dir}/${sampleID}.stats.txt  ]];
    then
        $samtools stats -@ ${thread} $bam > ${qc_dir}/${sampleID}.stats.txt;
    fi

    local r1=$(du $input_folder/${sampleID}_R1.fastq.gz -shL |awk '{print $1}');
    local r2=$(du $input_folder/${sampleID}_R2.fastq.gz -shL |awk '{print $1}');

    local raw_reads=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
    print(fh['summary']['before_filtering']['total_reads'])"`

    local clean_reads=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
    print(fh['summary']['after_filtering']['total_reads'])"`

    local raw_bases=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
    print(fh['summary']['before_filtering']['total_bases'])"`

    local clean_bases=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
    print(fh['summary']['after_filtering']['total_bases'])"`

    local qc_rate=`python3 -c "import json; \
    fh = json.load(open('$trim_dir/${sampleID}.trim.json', 'r')); \
    print(round(fh['summary']['before_filtering']['q30_rate']*100, 2))"`

    local trim_rate=`python3 -c "print(round(100-$clean_bases/$raw_bases*100, 2))"`

    local insert_size=$(awk -F"\t" '$2 == "insert size average:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);
    local insert_std=$(awk -F"\t" '$2 == "insert size standard deviation:" {print $3}' ${qc_dir}/${sampleID}.stats.txt);


    if [[  $bed != "empty"  ]];
    then
        if [[ ! -d $qc_dir/${sampleID} ]]; then
            mkdir $qc_dir/${sampleID};
        fi;

        if [[ ! -f $qc_dir/${sampleID}/coverage.report  ]];
        then
            $bamdst -p $bed -o $qc_dir/${sampleID} $bam;
        fi

        local mapping_rate=$(grep "Fraction of Mapped Reads" $qc_dir/${sampleID}/coverage.report | awk -F"\t" '{print $2}');
        local mean_depth=$(grep "Average depth" $qc_dir/${sampleID}/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
        local mean_dedup_depth=$(grep "Average depth(rmdup)" $qc_dir/${sampleID}/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
        local dup_rate=$(grep "Fraction of PCR duplicate reads" $qc_dir/${sampleID}/coverage.report |awk -F"\t" '{print $2}');
        local on_target=$(grep "Fraction of Target Reads in all reads" $qc_dir/${sampleID}/coverage.report |awk -F"\t" '{print $2}');

        local uniformity_01x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.1) count+=1} END {print count/NR*100}');
        local uniformity_02x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.2) count+=1} END {print count/NR*100}');
        local uniformity_05x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth*0.5) count+=1} END {print count/NR*100}');
        local uniformity_1x=$(less -S $qc_dir/${sampleID}/depth.tsv.gz | awk -v depth=${mean_dedup_depth} 'BEGIN {count=0} {if ($4 > depth) count+=1} END {print count/NR*100}');
    else
        local mapping_rate=".";
        local mean_depth=".";
        local mean_dedup_depth=".";
        local dup_rate=".";
        local on_target=".";
        local cov_1000x=".";
        local cov_2000x=".";
        local cov_3000x=".";
        local cov_4000x=".";
        local cov_5000x=".";
        local uniformity_01x=".";
        local uniformity_02x=".";
        local uniformity_05x=".";
        local uniformity_1x="."
    fi


    if [[ -n $umi_templ  ]];
    then
        if [[ ! -d $qc_dir/${sampleID}_umi ]]; then
            mkdir $qc_dir/${sampleID}_umi;
        fi;

        dedup_bam=${align_dir}/${sampleID}.sorted.dedup.bam

        if [[  $bed != "empty"  ]];
        then

            if [[ ! -f $qc_dir/${sampleID}_umi/coverage.report  ]];
            then
                $bamdst -p $bed -o $qc_dir/${sampleID}_umi $dedup_bam;
            fi
            
            local mean_dedup_depth=$(grep "Average depth" $qc_dir/${sampleID}_umi/coverage.report |head -n 1 |awk -F"\t" '{print $2}');
            local dup_rate=`python3 -c "print(round(100-$mean_dedup_depth/$mean_depth*100, 2))"`

            local cov_300x=$(less -S $qc_dir/${sampleID}_umi/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 300) count+=1} END {print count/NR*100}');
            local cov_500x=$(less -S $qc_dir/${sampleID}_umi/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 500) count+=1} END {print count/NR*100}');
            local cov_1000x=$(less -S $qc_dir/${sampleID}_umi/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 1000) count+=1} END {print count/NR*100}');
            local cov_1500x=$(less -S $qc_dir/${sampleID}_umi/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 1500) count+=1} END {print count/NR*100}');
            local cov_2000x=$(less -S $qc_dir/${sampleID}_umi/depth.tsv.gz | awk 'BEGIN {count=0} {if ($4 > 2000) count+=1} END {print count/NR*100}');
        fi
    fi

    echo "${sampleID},${r1}/${r2},${raw_reads},${raw_bases},${clean_reads},${clean_bases},\
    ${qc_rate},${trim_rate},${mapping_rate},${on_target},${mean_depth},${mean_dedup_depth},${dup_rate},${insert_size},${insert_std},\
    ${uniformity_01x},${uniformity_02x},${uniformity_05x},${uniformity_1x},\
    ${cov_300x},${cov_500x},${cov_1000x},${cov_1500x},${cov_2000x}" \
    >> $qc_dir/QC_summary.csv
}


# 4. Indel realignment
function run_realignment {
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- Indel Realignment";

    $sentieon driver -t ${thread} -r ${ref} -i $bam \
    --algo Realigner -k ${k1} -k ${k2} -k ${k3} \
    ${align_dir}/${sampleID}.realign.bam;
}


# 5. variant calling with TNscope
function run_tnscope {
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- variant calling";
    
    if [[  $mode == "matched"  ]];
    then
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
        
    elif [[  $mode == "single"  ]];
    then
        $sentieon driver -t ${thread} -r ${ref} \
        -i $bam \
        --algo TNscope \
        --tumor_sample ${sampleID} \
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
}


# 6. filter variants
function filter_vars {
    echo "LOGGING: $sampleID -- `date --rfc-3339=seconds` -- filter low quality variants";

    grep -v "SVTYPE=BND" ${tnscope_dir}/${sampleID}.raw.vcf > ${tnscope_dir}/${sampleID}.step0_snv.vcf

    # left normalisation
    $bcftools norm -m -both -f ${ref} \
    ${tnscope_dir}/${sampleID}.step0_snv.vcf \
    -o ${tnscope_dir}/${sampleID}.step1_norm.vcf;

    $bcftools annotate -x "FILTER/triallelic_site" ${tnscope_dir}/${sampleID}.step1_norm.vcf | \
    $bcftools filter -m + -s "low_qual" -e "QUAL < 10" | \
    $bcftools filter -m + -s "short_tandem_repeat" -e "RPA[0]>=5" | \
    $bcftools filter -m + -s "read_pos_bias" -e "FMT/ReadPosRankSumPS[0] < -5" | \
    $bcftools filter -m + -s "base_qual_bias" -e "FMT/BaseQRankSumPS[0] < -5" | \
    $sentieon util vcfconvert - ${tnscope_dir}/${sampleID}.step2_pre_filt.vcf

    # remove variants that non-pass & qual < 100
    awk '$0 ~ /^#/ || ($6 > 100 && $7 == "PASS")' \
    ${tnscope_dir}/${sampleID}.step2_pre_filt.vcf >\
    ${tnscope_dir}/${sampleID}.step3_filt.vcf;

    $bcftools filter -i "(FORMAT/AF[:0]) >= 0.0002 & (FORMAT/AD[:0]+AD[:1]) >= 1000" \
    ${tnscope_dir}/${sampleID}.step3_filt.vcf \
    -o ${tnscope_dir}/${sampleID}.step4_filt.vcf

}


# ?. generate mpileup
function gen_mpileup {
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- generate mpileup";

    if [[  $bed == "empty"  ]];
    then
        $samtools mpileup -B -Q 20 -C 50 -q 30 -d 500000 \
        -f $ref $bam > ${align_dir}/${sampleID}.realign.mpileup;
    else 
        $samtools mpileup -B -Q 20 -C 50 -q 30 -d 500000 -l $bed\
        -f $ref $bam > ${align_dir}/${sampleID}.realign.mpileup;
    fi
}


# ?. variant calling with Varscan
function run_varscan {
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- variant calling (Varscan)";
    java -Xmx3550m -Xms3550m -Xmn2g -jar $varscan \
    mpileup2snp ${align_dir}/${sampleID}.realign.mpileup \
    --output-vcf 1 \
    --min-coverage 100000 \
    --min-reads2 8 \
    --min-avg-qual 30 \
    --min-var-freq 0.0002 \
    --min-strandedness 0.0001 > \
    ${varscan_dir}/${sampleID}.raw.snp.vcf;

    java -Xmx3550m -Xms3550m -Xmn2g -jar $varscan \
    mpileup2indel ${align_dir}/${sampleID}.realign.mpileup \
    --output-vcf 1 \
    --min-coverage 100000 \
    --min-reads2 8 \
    --min-avg-qual 30 \
    --min-var-freq 0.0002 \
    --min-strandedness 0.0001 > \
    ${varscan_dir}/${sampleID}.raw.indel.vcf;
}


# 7. annotation
function anno_vars {
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- VEP Annotation";
    $vep --dir ${vep_dir} --cache --offline --cache_version ${cache_version} \
    --assembly GRCh37 --format vcf --fa ${ref} --force_overwrite --vcf \
    --gene_phenotype --use_given_ref --refseq --check_existing \
    --hgvs --hgvsg --transcript_version --max_af \
    --vcf_info_field ANN -i ${tnscope_dir}/${sampleID}.step4_filt.vcf \
    -o ${tnscope_dir}/${sampleID}.step5_anno.vcf;

    python3 $anno_hgvs ${tnscope_dir}/${sampleID}.step5_anno.vcf \
    $clinic_transcripts $refflat -o $tnscope_dir/$sampleID.step6_anno.vcf;
}


# 8. keep hotspot variants only
function filter_vars2 {
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- keep hotspot variants only";

    # remove variants with max-MAF > 0.1% 
    python3 $hotspot_filt -i $tnscope_dir/$sampleID.step6_anno.vcf -c $cbioportal > $tnscope_dir/$sampleID.step7_filt.vcf

}


# 9. run genefuse
function run_genefuse {
    echo "LOGGING: ${sampleID} -- `date --rfc-3339=seconds` -- running genefuse";

    $genefuse -t ${thread} -r ${ref} \
    -f $fusion_templ \
    -1 $trim_dir/$sampleID.umi_consensus.fq.gz \
    -h $fusion_dir/${sampleID}.fusion.html \
    -j $fusion_dir/${sampleID}.fusion.json  1> $fusion_dir/${sampleID}.fusion.txt
}



# ---------------------------------------------------------------------- #
# ---------------------------  Pipeline  ------------------------------- #
# ---------------------------------------------------------------------- #

export SENTIEON_LICENSE=${sentieon_license};

if [[  $mode == 'matched' ]]; then
    for ifile in $input_folder/*_normal_R1.fastq.gz;
    do 
        sample_basename=`basename ${ifile%%"_normal"*}`

        if [[  $do_trim == "on" ]];
        then
            sampleID=${sample_basename}_tumor;
            run_trim;

            sampleID=${sample_basename}_normal;
            run_trim;
        fi

        # step2 - align & sort
        if [[  $do_align == "on"  ]] 
        then
            sampleID=${sample_basename}_tumor;
            run_alignment;

            sampleID=${sample_basename}_normal;
            run_alignment;
        fi

        if [[  -n $umi_templ  ]];
        then
            tumor_bam=${align_dir}/${sample_basename}_tumor.sorted.dedup.bam
            normal_bam=${align_dir}/${sample_basename}_normal.sorted.dedup.bam
        else
            tumor_bam=${align_dir}/${sample_basename}_tumor.sorted.bam
            normal_bam=${align_dir}/${sample_basename}_normal.sorted.bam
        fi

        # step4 - quality control
        if [[  $do_qc == "on"  ]];
        then
            sampleID=${sample_basename}_tumor;
            run_qc;

            sampleID=${sample_basename}_normal;
            run_qc;
        fi

        # step5 - Indel Realignment
        if [[  $do_realign == "on"  ]];
        then
            sampleID=${sample_basename}_tumor;
            run_realignment;

            sampleID=${sample_basename}_normal;
            run_realignment;
        fi

        tumor_bam=${align_dir}/${sample_basename}_tumor.realign.bam
        normal_bam=${align_dir}/${sample_basename}_normal.realign.bam


        # step6 - variant calling
        if [[  $do_tnscope == "on"  ]];
        then
            run_tnscope;
        fi

        if [[ $do_mpileup == "on"  ]];
        then
            gen_mpileup;
        fi

        if [[ $do_varscan == "on"  ]];
        then
            run_varscan;
        fi

        # step7 - normalisation + remove low quality variants
        if [[  $do_filt == "on"  ]];
        then
            filter_vars;
        fi


        # step8 - annotation
        if [[  $do_anno == "on"  ]];
        then
            anno_vars;
        fi

        if [[  $do_filt2 == "on"  ]];
        then
            filter_vars2;
        fi
    done

elif [[  $mode == 'single'  ]]; then
    for ifile in $input_folder/*_R1.fastq.gz;
    do
        sampleID=`basename ${ifile%%"_R1"*}`

        # step1 - trim reads
        if [[ $do_trim == "on" ]];
        then
            run_trim;
        fi

        # step2 - align & sort
        if [[ $do_align == "on" ]];
        then
            run_alignment;
        fi

        # step5 - quality control
        if [[  $do_qc == "on"  ]];
        then
            run_qc;
        fi

        # step6 - Indel Realignment
        if [[  $do_realign == "on"  ]];
        then
            run_realignment;
        fi

        bam=${align_dir}/${sampleID}.realign.bam


        # step7 - variant calling
        if [[ $do_tnscope == "on"  ]];
        then
            run_tnscope;
        fi

        if [[ $do_mpileup == "on"  ]];
        then
            gen_mpileup;
        fi

        if [[ $do_varscan == "on"  ]];
        then
            run_varscan;
        fi

        if [[ $do_filt == "on"  ]];
        then
            filter_vars;
        fi

        if [[  $do_anno == "on"  ]];
        then
            anno_vars;
        fi

        if [[  $do_filt2 == "on"  ]];
        then
            filter_vars2;
        fi

        if [[  $do_fusion == "on"  ]];
        then
            run_genefuse;
        fi
    done
fi

echo "LOGGING: `date --rfc-3339=seconds` -- Analysis finished";
