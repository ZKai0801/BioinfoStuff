import sys
import os
import pandas as pd
import argparse


def count_frac(stdin, samples, contSites, ofname):
    header = ["chrom", "pos", "ref", "alt"]
    header += samples

    rows = []
    index = 0

    # this pos does not contain reads, therefore loop over stdin
    # until reads can be found
    for line in stdin:
        nline = line.strip().split("\t")

        while index <= contSites.shape[0] - 1:
            if ((nline[0] == contSites.loc[index, "chrom"]) and 
                (int(nline[1]) == contSites.loc[index, "pos"]) and 
                (nline[2].upper() == contSites.loc[index, "ref"])):
                
                alt = contSites.loc[index, "alt"]
                each_row = [nline[0], nline[1], nline[2], alt]

                samples_index = list(range(4, len(nline)+1, 3))
                for i in samples_index:
                    depth = len(nline[i])
                    alt_depth = nline[i].upper().count(alt)
                    frac = alt_depth/depth
                    each_row.append(frac)
                
                rows.append(each_row)
                index += 1
                break
            else:
                index += 1


    df = pd.DataFrame(data = rows, columns = header)
    df.to_csv(ofname, sep="\t", index=False)
    return df


def determine_contam(samples, df, ofname, homo_frac=0.15, contam_frac=0.05):
    with open(ofname, "w") as fw:
        header = ["sampleID", "contaminated_sites", "all_homo_sites", "Frac"]
        fw.write("\t".join(header)+"\n")
        for sample in samples:
            tmp = df[(df[sample] != 0) & ((df[sample] < homo_frac) | (df[sample] > 1-homo_frac))]
            homosites = tmp.shape[0]
            tmp = tmp[(tmp[sample] >= contam_frac) & (tmp[sample] <= 1-contam_frac)]
            contamsite = tmp.shape[0]
            af = 0 if homosites == 0 else contamsite/homosites
            res = [sample, str(contamsite), str(homosites), str(af)]
            fw.write("\t".join(res)+"\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="avail_sites.bed")
    parser.add_argument("-i", "--ids", nargs="*", help="sampleIDs", required=True)
    parser.add_argument("-o", "--ofname", required=True)
    args = parser.parse_args()

    samples = [os.path.basename(i).split(".")[0] for i in args.ids]
    contSites = pd.read_csv(args.input, sep = "\t", usecols = [0, 2, 3, 4], header = None)
    contSites.columns = ["chrom", "pos", "ref", "alt"]

    ofname = args.ofname + "_adjAF.tsv"
    df = count_frac(sys.stdin, samples, contSites, ofname)

    ofname = args.ofname + "_adj_results.tsv"
    determine_contam(samples, df, ofname)
