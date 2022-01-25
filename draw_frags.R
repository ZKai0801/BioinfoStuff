library(ggplot2)
library(dplyr)
library(data.table)
library(argparse)

parser <- argparse::ArgumentParser()
parser$add_argument("input", help = "input file should be the 5mb binned csv file")
args <- parser$parse_args()

bins5mb <- args$input

# for test
# bins5mb <- "\\\\172.16.11.242/kai/projects/fragmentome/results/tumor_5mb.csv"

sampleID <- gsub("\\..+", "", basename(bins5mb))
ofname  <- file.path(dirname(bins5mb), paste0(sampleID, ".frag.png"))

bins5mb <- fread(bins5mb)

arms <- c(rbind(paste0(c(seq(1, 22), "X"), "p"), 
                paste0(c(seq(1, 22), "X"), "q")))


bins5mb <- bins5mb %>% 
    mutate(arm = factor(arm, levels = arms)) %>% 
    arrange(arm)  %>%
    group_by(id) %>%  
    mutate(bin = row_number())

breaks <- bins5mb %>% group_by(arm) %>% filter(row_number() == ceiling(n()/2))  %>% pull(bin)
labels <- bins5mb %>% group_by(arm) %>% filter(row_number() == ceiling(n()/2))  %>% pull(arm)

bins5mb  %>% 
    rowwise() %>% 
    mutate(ratio.cor = short.cor/long.cor)  %>% 
    group_by(id) %>% 
    mutate(ratio.centered = scale(ratio.cor, scale = F)[, 1]) %>% 
    ggplot(aes(x = bin, y = ratio.centered, group = id)) + 
    geom_line(color = "grey50", size = .3, alpha = .7) + 
    theme_bw() + 
    labs(x = "", y = "ratio.center", title = sampleID) + 
    scale_x_continuous(breaks = breaks, labels = labels) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ggsave(file = ofname,
       width = 14, height = 5)  
