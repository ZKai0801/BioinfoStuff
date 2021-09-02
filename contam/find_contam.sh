#!/usr/bin/bash

# ----------------------------------------------------
# Usage:
# [kai@admin] bash find_contam.sh *raw.vcf *bam
# 
# make sure config lines are set correctly;
# and samtools bcftools bgzip are added to your path
# ----------------------------------------------------

ref="/public/database/GATK_Resource_Bundle/hg19/ucsc.hg19.fasta"
contam_py="/public/home/kai/projects/contam/find_contam.py"


args=$@

vcfs=()
bams=()

for fname in $args;
do
    if [[ $fname == *"vcf" ]];
    then
        vcfs+=($fname)
    elif [[ $fname == *"bam" ]];
    then
        bams+=($fname)
    else
        echo "Warning: pls check you input args";
    fi
done


gzs=()
for vcf in ${vcfs[@]};
do
    bgzip $vcf;
    tabix -p vcf $vcf.gz;
    gzs+=($vcf.gz)
done


bcftools merge --force-samples -o merged.vcf ${gzs[@]}

# keep only valid snp sites
# indel, low-quality sites, str are all removed from analysis
grep -v "#\|STR" merged.vcf |\
 grep "PASS\|t_lod\|clustered_events" |\
 awk -F"\t" '!($7 ~ /;/) && (length($4) == 1) && (length($5) == 1) {OFS="\t"; print $1, $2-1, $2, $4, $5}' \
 > avail_sites.bed;


# count the secondary altered allele in each samples
samtools mpileup -Q 20 -C 50 -q 20 -d 30000 \
-l avail_sites.bed \ 
-f $ref ${bams[@]} | \
python3 $contam_py \
avail_sites.bed -i ${bams[@]} -o contam


# remove intermediate files
rm merged.vcf avail_sites.bed;
for gz in ${gzs[@]};
do
    bgzip -d $gz;
    rm $gz.tbi;
done

