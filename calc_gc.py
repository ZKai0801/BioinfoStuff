"""
calculate GC contents for each genomic position
GC content = count(G + C) / count(A + T + G + C) * 100%
By default, GC contents for a given snp, is calculated from its surrounding 50bp
"""

import sys
from pyfaidx import Fasta
import argparse
from collections import defaultdict



def read_pos(infile):
    """
    input file should be a tab-delimited file containing two columns, 
    with chrom on the first column and position on the second column
    """
    sites = []
    with open(infile, "r") as fh:
        for line in fh:
            chrom, pos = line.strip().split("\t")[:2]
            sites.append((chrom, int(pos)))
    return sites


def calc_gc(site, fasta, surrounding):
    forth = round(surrounding / 2)
    back = surrounding - forth

    seq = fasta[site[0]][site[1]-forth: site[1]+back].seq.upper()
    gc = round(seq.count("G") + seq.count("C") / surrounding * 100, 4)
    return gc


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-s", "--surrounding", default=50)
    parser.add_argument("-o", "--output", type=argparse.FileType("w"), default=sys.stdout)
    args = parser.parse_args()

    fasta = Fasta(args.fasta)

    print("chrom", "pos", "GC_contents(%)", sep="\t", file=args.output)
    for site in read_pos(args.input):
        gc = calc_gc(site, fasta, args.surrounding)
        print(site[0], site[1], gc, sep="\t", file=args.output)
   