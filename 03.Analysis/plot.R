#!/bin/env Rscript
# parse parameter ---------------------------------------------------------
library(argparser, quietly=TRUE)
# Create a parser
p <- arg_parser("a program for background analysis")

## Add command line arguments
p <- add_argument(parser = p, arg = "--createParameter", short = "-C", help = "create a Parameter csv file", flag = TRUE)
p <- add_argument(parser = p, arg = "--parameter", help = "Parameter csv file containing parameters, if the same parameter is listed in CMD line, this CMD line parameter will be omitted", type = "character")

p <- add_argument(parser = p, arg = "--donor", short = "-d", help = "donor parent name", type = "character")
p <- add_argument(parser = p, arg = "--recurrent", short = "-r", help = "recurrent parent name", type = "character")
p <- add_argument(parser = p, arg = "--sample", short = "-s", help = "sample parent name", type = "character")

p <- add_argument(parser = p, arg = "--chromosome", short = "-c", help = "tab-separated file contain chromosome ID and label shown in figure. if not assigned, all chromosome in genome will be shown", type = "character")
p <- add_argument(parser = p, arg = "--length", short = "-l", help = "tab-separated file contain chromosome ID and length. if not assigned, the max position of SNP for each chromosome in GATK Table file will be regarded as chromosome length", type = "character")

p <- add_argument(p, "--width", help = "Plot width", type = "numeric", default = 8)
p <- add_argument(p, "--height", help = "Plot height", type = "numeric", default = 7)

# Parse the command line arguments
argv <- parse_args(p)

library(tidyverse)
library(cowplot)

input <- argv$input
chromosome <- argv$chromosome
length <- argv$length

# 
if (is.na(argv$parameter)) {
  cat("Get parameters from CMD line.\n")
  parameter <- data.frame(sample = argv$sample, donor = argv$donor, recurrent = argv$recurrent, 
                          minDdp = argv$minDdp, maxDdp = argv$maxDdp, 
                          minRdp = argv$minRdp, maxRdp = argv$maxRdp, 
                          width = argv$width, height = argv$height)
} else {
  cat("Get parameters from parameter csv file.\n")
  parameter <- read.csv(argv$parameter, header = T)
}

for (i in seq_along(parameter[,1])) {
  sample <- parameter$sample[i]
  donor <- parameter$donor[i]
  recurrent <- parameter$recurrent[i]
  width <- parameter$width[i]
  height <- parameter$height[i]
  
  snpbinner <- read_csv(file = paste("./", sample, "/", sample, ".genotype.csv", sep = "")) %>% rename(CHROM = chr)
  bins <- read_csv(file = paste("./", sample, "/", sample, ".bg.csv", sep = "")) %>% rename(CHROM = chr)
  
  if (is.na(chromosome)) {
    chr <- snpbinner %>% distinct(CHROM) %>% mutate(LABEL = CHROM, y = rev(seq_along(rownames(.))*3))
  } else {
    chr <- read_tsv(chromosome, col_names = c("CHROM", "LABEL"), show_col_types = FALSE) %>% mutate(y = rev(seq_along(rownames(.))*3))
  }
  if (is.na(length)) {
    len <- snpbinner %>% group_by(CHROM) %>% summarise(Len = max(`position(bp)`)) %>% right_join(chr, by = "CHROM")
  } else {
    len <- read_tsv(file = length, col_names = c("CHROM", "Len"), show_col_types = FALSE) %>% right_join(chr, by = "CHROM")
  }
  
  genoforplot <- snpbinner %>% 
    select(CHROM, x = `position(bp)`, genotype = all_of(sample)) %>% 
    right_join(chr, by = "CHROM")
  binsforplot <- bins %>% 
    select(CHROM, START = start, END = end, genotype) %>% 
    mutate(genotype = if_else(genotype == "Donor", "a", 
                              if_else(genotype == "Heterozygous", "h", "-"))) %>% 
    filter(genotype == "a" | genotype == "h") %>% 
    left_join(chr, by = "CHROM")
  
  p <- ggplot(len) + 
    geom_rect(mapping = aes(xmin = 0, xmax = Len, ymin = y-1, ymax = y+1), fill = "gray60", color = "black") + 
    #geom_segment(data = genoforplot, mapping = aes(x = x, xend = x, y = y-1, yend = y, color = genotype), linewidth = 0.1) + 
    geom_rect(data = binsforplot, mapping = aes(xmin = START, xmax = END, ymin = y, ymax = y+1, fill = genotype)) + 
    scale_x_continuous(breaks = seq(0, max(len$Len) + 1000000, 10000000), 
                       labels = paste(seq(0, (max(len$Len) %/% 1000000) + 1, 10), "Mb", sep = " "),
                       limits = c(0, (max(len$Len+1000000) %/% 1000000) * 1000000), 
                       expand = c(0.025, 0)) + 
    scale_y_continuous(breaks = chr$y, labels = chr$LABEL, expand = c(0.025, 0)) + 
    scale_fill_manual(breaks = c("b", "h", "a"), values = c("blue", "green", "red")) + 
    scale_color_manual(breaks = c("b", "h", "a"), values = c("blue", "green", "red"), labels = c(recurrent, "Het", donor)) + 
    guides(x = guide_axis(cap = TRUE), y = guide_axis(cap = TRUE), fill = "none", 
           color = guide_legend(override.aes = list(linewidth = 3))) + 
    labs(x = NULL, y = NULL, color = "Genotype", fill = NULL) + 
    theme_half_open()
  if (require(ggrastr, quietly = TRUE)) {
    p <- p + rasterise(geom_segment(data = genoforplot, mapping = aes(x = x, xend = x, y = y-1, yend = y, color = genotype), linewidth = 0.1), dpi = 500)
  } else {
    p <- p + geom_segment(data = genoforplot, mapping = aes(x = x, xend = x, y = y-1, yend = y, color = genotype), linewidth = 0.1)
  }
  ggsave(p, filename = paste("./", sample, "/", sample, ".bg.corrected.png", sep = ""), width = 8, height = 7, dpi = 500)
  ggsave(p, filename = paste("./", sample, "/", sample, ".bg.corrected.pdf", sep = ""), width = 8, height = 7)
}





