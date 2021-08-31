import pysam 
import argparse
import os

def read_vcf(vcf):
    """
    keep variants in autosomal chroms,
    and with qual > 50
    return: (chrom, pos, alt)
    """
    chroms = ["chr" + str(i) for i in range(1,23)]
    variants = []
    with open(vcf, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.strip().split()
            if (line[0] in chroms) and (float(line[5]) >= 50):
                variants.append((line[0], line[1], line[4]))

    return variants


def find_reads(bam, variants, ofname):
    """ mutant reads are save to a file  """
    with open(ofname, "w") as fw:
        for var in variants:
            with pysam.AlignmentFile(bam, "rb") as samin:
                for read in samin.fetch(var[0], int(var[1]), int(var[1])+1):
                    if read.query_sequence[int(var[1]) - read.pos -1: int(var[1]) - read.pos] == var[2]:
                        fw.write(read.query_name+"\n")


def filt_reads(bam, ofname, outbam):
    """ filter bam by read ID & index the bam """
    picard = "/public/software/picard.jar"
    os.system(f"java -jar {picard} FilterSamReads I={bam} O={outbam} READ_LIST_FILE={ofname} FILTER=excludeReadList")
    os.system(f"samtools index -@ 8 {outbam}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bam", help="path to bam")
    parser.add_argument("vcf", help="path to vcf")
    parser.add_argument("ofname", help="prefix to output file")
    args = parser.parse_args()

    vcf = args.vcf
    bam = args.bam
    oprefix = args.ofname

    variants = read_vcf(vcf)
    find_reads(bam, variants, oprefix+"_mutreads.txt")
    filt_reads(bam, oprefix+"_mutreads.txt", oprefix+".filt.bam")