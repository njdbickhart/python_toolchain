#!/usr/bin/env Rscript

library(ggridges)
library(ggplot2)
library(gridExtra)
library(dplyr)

outfile <- "summary_taxonomic_plot.pdf"
# Should be alternating "names" and files
# ie. "ASM1" "asm1.blobtools.table.tab"
args = commandArgs(trailingOnly=TRUE)


# Input file tabs
# LEN GC KING
alltabs <- list()

for(l : seq(1, length(args))){
  if(l %% 2 == 0){
    next
  }
  df <- read.delim(args[l + 1], header = TRUE)
  df <- df %>% mutate(ASM = c(args[l]))
  alltabs[args[l]] <- df
}

fulldata <- do.call(rbind, alltabs)

fulldata$ASM <- as.factor(fulldata$ASM)
fulldata$KING <- factor(fulldata$KING, levels = c("no-hit", "Viruses", "Eukaryota", "Archaea", "Bacteria"))

pdf(file=outfile, useDingBats=FALSE)
p1 <- ggplot(fulldata, aes(y=KING, x=LEN, fill=ASM)) + geom_density_ridges(aes(alpha = 0.4)) + theme_bw() + theme(axis.title=element_text(size=13, face="bold"), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_brewer(palette="Dark2") + scale_x_log10(breaks=c(1000, 10000, 100000, 1000000, 6000000), limits=c(100, 6000000), labels=c("1000", "10,000", "100,000", "1,000,000", "6,000,000")) + xlab(label = "Log10 Contig Lengths (bp)") + ylab(label= "Contig Superkingdom Taxonomic Assignment")
p2 <- ggplot(fulldata, aes(x=ASM, y=GC, fill=ASM)) + geom_violin(trim=FALSE) + geom_boxplot(width=0.1, outlier.shape = NA, fill = "white") + theme_bw() + theme(axis.title=element_text(size=13, face="bold"), axis.text.x = element_text(size=11), axis.text.y = element_text(size=11)) + scale_fill_brewer(palette="Dark2") + labs(x="Assembly", y = "Avg GC ratio per contig")
grid.arrange(p1, p2, nrow=2)
dev.off()
