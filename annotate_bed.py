__author__ = 'Kai'
__date__ = '21/03/2019'
__email__ = 'zhentian.kai@outlook.com'

__doc__ = '''
    Annotate bed file based on refFlat.txt file from UCSC
    
    Usage:
    [kai@admin]$ python3 annotate_bed.py -i [unannotated_bed] -o [output_filename] -f [refFlat.txt] -l [gene_list]
    
    Required: 1) refFlat.txt; 2) BED file that you want to annotate
    optional: 3) gene_list of that panel
'''

import argparse
import os


def read_bed(bed):
    '''
        read BED file into an array: [(chr, start, end), ...]
    '''
    bedinfo = []
    with open(bed, 'r') as f:
        for line in f:
            line = line.strip().split()
            bedinfo.append((line[0], int(line[1]), int(line[2])))
    return bedinfo


def read_gene_list(gene_list):
    '''
        gene_list should be a text file that each line represents a single gene
        @return: genes -- [gene1, gene2, ...]
    '''
    genes = []
    try:
        with open(gene_list, 'r') as f:
            for line in f:
                genes.append(line.strip())
        return genes
    except:
        return None



def read_refFlat(refFlat, genes):
    '''
        read refFlat.txt into a dictionary: {(chr, start, end):gene, ... }
    '''
    annotation = {}
    with open(refFlat, 'r') as f:
        for line in f:
            line = line.strip().split()
            if genes:
                if line[0] in genes:
                    annotation[(line[2], int(line[4]), int(line[5]))] = line[0]
            else:
                annotation[(line[2], int(line[4]), int(line[5]))] = line[0]
    return annotation


def rewrite_bed(outname, bedinfo, annotation):
    '''
        Annotate BED intervals if that interval overlap with the interval from refFlat.txt
        There might be more than one annotation for each interval
    '''
    fw = open(outname, 'w')
    for interval in bedinfo:
        # interval - (chr, start, end)
        annotated_genes = set()
        # set flag to check if annotation can be found eventually
        flag = False
        for anno in annotation:
            # check if on the same chromosome
            if interval[0] == anno[0]:
                # check if overlapped
                if anno[1] <= interval[1] <= anno[2] or anno[1] <= interval[2] <= anno[2] or interval[1] <= anno[1] <= anno[2] <= interval[2]:
                    annotated_genes.add(annotation[anno])
                    flag = True
        if flag:
            fw.write('{}\t{}\t{}\t{}\n'.format(interval[0], interval[1], interval[2], ','.join(annotated_genes)))
        else:
            # write a dot if no annotation can be found
            fw.write('{}\t{}\t{}\t.\n'.format(interval[0], interval[1], interval[2]))
    fw.close()



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help = 'Path to the unannotated BED file', required = True)
    parser.add_argument('-o', '--output', help = 'The output file name' )
    parser.add_argument('-f', '--refFlat', help = 'Path to the refFlat.txt', required = True)
    parser.add_argument('-l', '--list', help = 'Path to a text file containing genes for the BED')
    args = parser.parse_args()
    
    bedinfo = read_bed(args.input)
    annotation = read_refFlat(args.refFlat, read_gene_list(args.list))
    if args.output:
        rewrite_bed(args.output, bedinfo, annotation)
    else:
        outname = os.path.basename(args.input).split('.')[0] + '.annotated.bed'
        rewrite_bed(outname, bedinfo, annotation)





