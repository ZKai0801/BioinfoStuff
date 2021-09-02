# README

A simple script detects cross-samples contamination. 

In general, this script merges all variants from all samples in a single run into a VCF. For each variants in the VCF,  the script calculates the fractions of the secondary most common alleles for each samples. For each samples, we extracts all its homozygous sites and determined if a single site is contaminated via evaluating its calculated allele fractions. Contaminated samples will be expected to have high fraction of contaminated homozygous sites. 

