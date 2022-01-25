suppressMessages(library(data.table))
suppressMessages(library(argparse))

parser <- argparse::ArgumentParser()
parser$add_argument("-i", "--input", 
                    help = "the filtered fragments")
parser$add_argument("-r", "--ref", 
                    default = "", 
                    help = "the reference GC content")
parser$add_argument("-t", "--templ", 
                    default = "",
                    help = "template bins")
parser$add_argument("-o", "--ofname",
                    help = "output filename")
args <- parser$parse_args()


ref <- args$ref
templ <- args$templ
frag <- args$input

if (!is.null(args$ofname)){
    ofname <- args$ofname
} else {
    ofname <- file.path(dirname(frag),
                        paste0(gsub("\\..+", "", basename(frag)), ".5mb.csv"))
}

# for test
# setwd("\\\\172.16.11.242/kai/projects/fragmentome/")
# ref <- "sources/target20.tsv"
# templ <- "sources/templ_bins.tsv"
# frag <- "results/MRD211129-505/fragmentome/MRD211129-505_tumor.step2_filt.bed"


ref <- fread(ref)
templ <- fread(templ)
frag <- fread(frag, 
              col.names = c("chr", "start", "end", "mapq", "fraggc"))


# ------------- adjust for GC content -------------- #
#' reference GC baseline is calculated from a cohort
#' of healthy people.

gc_dt <- frag[,.(n=.N), by = .(fraggc, chr)]
gc_dt <- gc_dt[fraggc >= 0.2 & fraggc <= 0.8]
gc_dt <- gc_dt[order(fraggc, chr)]

setkey(gc_dt, fraggc, chr)
setkey(ref, gc, chr)

gc_dt <- gc_dt[ref]
gc_dt[, w:=gcmed/n]

frag[gc_dt, on = .(chr, fraggc), weight := i.w]
frag <- frag[!is.na(weight)]


# get overlapped fragments for each bins
templ <- templ[, chr := factor(chr, paste0("chr", c(1:22, "X")))]
frag <- frag[, chr := factor(chr, paste0("chr", c(1:22, "X")))]

setkey(frag, chr, start, end)
setkey(templ, chr, start, end)

chroms <- paste0("chr",c(1:22, "X"))

fragbins <- foverlaps(frag[chr %in% chroms],
                      templ,
                      type = "within",
                      nomatch = NULL)

fragbins <- fragbins[, width:= i.end - i.start]


# for each bin, calculate the number of short/long fragments
bins <- fragbins[, .(arm = unique(arm), gc = gc[1], map = map[1],
                     short = sum(width >= 100 & width <= 150),
                     long = sum(width > 150 & width <= 250),
                     short.cor = sum(weight[width > 100 & width <= 150]),
                     long.cor = sum(weight[width > 150 & width <= 250]),
                     ultrashort = sum(width < 100),
                     ultrashort.cor = sum(weight[width < 100]), 
                     multinucs = sum(width > 250),
                     multinucs.cor = sum(weight[width > 250]),
                     mediansize = as.integer(median(width)),
                     frag.gc = mean(fraggc)),
                 by = .(chr, start, end)]

# 
setkey(bins, chr, start, end)
setkey(templ, chr, start, end)

bins <- bins[templ]
bins <- bins[is.na(i.gc), which(grepl("short|long|multi", colnames(bins))):=0]
bins[, `:=` (gc = i.gc, map = i.map, arm = i.arm)]
bins[,which(grepl("^i.", colnames(bins))):=NULL]
bins[, bin := 1:.N]
setcolorder(bins, c("chr", "start", "end", "bin"))

fwrite(bins, file = ofname)
