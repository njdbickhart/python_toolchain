#!/usr/bin/env Rscript

library(ggridges)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Should be alternating "names" and files
# ie. "ASM1" "asm1.blobtools.table.tab"
args = commandArgs(trailingOnly=TRUE)

first <- "figures/contig_facet_cov_by_comp_scatter.pdf"
second <- "figures/contig_top80_cov_by_comp.pdf"



# Input file tabs
# Contig, Coverage, Completeness, Assembly
#alltabs <- list()
fulldata <- data.frame(Contig=NA, Coverage=NA, Completeness=NA, Contamination=NA, Assembly=NA)[-1,]

for(l in seq(1, length(args))){
  if(l %% 2 == 0){
    next
  }
  print(args[l])
  df <- read.delim(args[l + 1], header = TRUE)
  df <- df %>% mutate(ASM = c(args[l]))
  print(head(df))
  fulldata <- bind_rows(fulldata, df)
}

fulldata$Assembly <- as.factor(fulldata$Assembly)

pdf(file=first, useDingbats = FALSE)
ggplot(data = fulldata, aes(x=Completeness, y=Coverage, colour=Assembly)) + geom_point() + scale_colour_brewer(palette='Set2') + theme_bw() + scale_y_log10() + facet_wrap(~ Assembly)
dev.off()

pdf(file=second, useDingbats = FALSE)
top80 <- fulldata[fulldata$Completeness >= 80,]
ggplot(data = top80, aes(x=Completeness, y=Coverage, colour=Assembly)) + geom_point() + scale_colour_brewer(palette='Set2') + theme_bw() + scale_y_log10() + facet_wrap(~ Assembly)
dev.off()
