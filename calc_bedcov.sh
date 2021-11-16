#!/usr/bin/bash

bed=$1
bam=${@:2}

ls $bam | awk 'function basename(file, a, n) {
        n = split(file, a, "/")
        return a[n]
    } BEGIN {printf "chrom\tstart\tend\ttarget"}; {printf "\t%s", basename($1)}; END {printf "\n"}' 

samtools bedcov $bed $bam |awk -F"\t" '{printf "%s\t%s\t%s\t%s\t", $1, $2, $3, $4; for(i=5;i<=NF;i++){printf "%s\t", $i/($3-$2)}; printf "\n"'} 