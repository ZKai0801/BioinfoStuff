#' Run QDNAseq for 30x WGS without paired samples

library(QDNAseq)
library(argparse)

parser <- argparse::ArgumentParser()
parser$add_argument("-i", "--input", required = T, help = "input bam file")
parser$add_argument("-t", "--threads", default = 8)
parser$add_argument("-o", "--ofname")
args <- parser$parse_args()


if (!is.null(args$ofname)){
    ofname <- args$ofname
} else {
    ofname <- file.path(dirname(args$input), 
                        gsub("\\..+", "", basename(args$input)))
}


future::plan("multisession", workers = args$threads)

bins <-  getBinAnnotations(binSize=10)

# get read counts
rc <- binReadCounts(bins, bamfiles = args$input)

# rm blacklisted regions & apply gc/mappability correction
filt_rc <- applyFilters(rc, residual = T, blacklist = T)
filt_rc <- estimateCorrection(filt_rc)
correct_rc <- correctBins(filt_rc)
normalised_rc <- normalizeBins(correct_rc)
smooth <- smoothOutlierBins(normalised_rc)

exportBins(smooth, file = paste0(ofname, ".txt"))

# CBS segmentation
seg <- segmentBins(smooth)
calls <- callBins(seg, ncpus = args$threads)
exportBins(calls, format = "vcf")


file.rename(from = list.files(pattern = "*vcf"),
            to = file.path(dirname(ofname), 
                           list.files(pattern = "*vcf")))

png(filename = paste0(ofname, ".png"), width = 2000)
plot(calls)
dev.off()

