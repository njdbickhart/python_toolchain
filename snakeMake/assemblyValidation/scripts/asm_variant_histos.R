#!/usr/bin/env Rscript
# Sections cribbed from Mike Schatz and Maria Nattestad's AssemblyTics

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RColorBrewer))

option_list <- list(
  make_option(c("-i","--input"), type="character", default=NULL,
              help="Bed file from between alignments program",
              dest="input_filename"),
  make_option(c("-o","--output"), type="character", default="out",
              help="output filename prefix [default %default]",
              dest="output_prefix"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print out all parameter settings [default]"),
  make_option(c("-q", "--min-length"), type="numeric", default=0,
              help="The minimum length to plot",
              dest="min_query_aln"),
  make_option(c("-m", "--max-length"), type="numeric", default=1000000,
              help="The maximum length to plot",
              dest="max_align"),
  make_option(c("-p", "--to_png"), action="store_true", default=FALSE,
              help="Plot as a PNG instead of a PDF [default: PDF]",
              dest="to_png")
)

options(error=traceback)

parser <- OptionParser(usage = "%prog -i alignments.coords -o out [options]",option_list=option_list)
opt = parse_args(parser)

theme_set(theme_bw(base_size = 12))

color_palette_name <- "Dark2"
big_palette<-brewer.pal(6,"Dark2")

# Nature-style formatting for publication using commas (e.g.: 7,654,321)
comma_format<-function(num) {
    formatC(abs(num),format="f",big.mark=",",drop0trailing = TRUE)
}

bed <- read.csv(opt$input_filename, sep="\t", quote='', header=TRUE)

names(bed)[1:11] <- c("chrom","start","stop","name","size","strand","type","ref.dist","query.dist","contig_position","method.found")

bed$type <- recode(bed$type, "Repeat_expansion"="Repeat expansion", "Repeat_contraction"="Repeat contraction", "Tandem_expansion"="Tandem expansion", "Tandem_contraction"="Tandem contraction")

types.allowed <- c("Insertion","Deletion","Repeat expansion","Repeat contraction","Tandem expansion","Tandem contraction")
bed$type <- factor(bed$type, levels = types.allowed)

var_size_cutoffs <- c(opt$min_query_aln, 50, 1000, opt$max_align)
var_size_cutoffs <- var_size_cutoffs[var_size_cutoffs>=opt$min_query_aln & var_size_cutoffs<=opt$max_align]


for (i in seq(1,length(var_size_cutoffs)-1)) {
    min_var <- var_size_cutoffs[i]
    max_var <- var_size_cutoffs[i+1]
    if (opt$verbose){
        print(paste0("Working on variant sizes: ", min_var, " and ", max_var))
    }
    if (min_var < opt$max_align && max_var > opt$min_query_aln)
    {
        types_to_plot = types.allowed
        filtered_bed <- bed[bed$size>=min_var &
                        bed$size<=max_var &
                        bed$type %in% types_to_plot,]
        filtered_bed$type <- factor(filtered_bed$type,levels=types_to_plot)
        binwidth <- max_var/100
        if (binwidth < 1) {
            binwidth <- 1
        }
        if (opt$verbose){
            summary(filtered_bed)
        }

        if (nrow(filtered_bed)>0) {
            if (opt$to_png) {
                if (opt$verbose){ print("Printing png...")}
                png(paste(opt$output_prefix,".", format(min_var, scientific=FALSE), "-",format(max_var, scientific=FALSE), ".png", sep=""),1000,1000,res=200)
            } else {
                if (opt$verbose){ print("Printing pdf...")}
                pdf(paste(opt$output_prefix,".", format(min_var, scientific=FALSE), "-",format(max_var, scientific=FALSE), ".pdf", sep=""))
            }

            print(ggplot(filtered_bed,aes(x=size, fill=type)) +
                geom_histogram(binwidth=binwidth) +
                scale_fill_manual(values=big_palette,drop=FALSE) +
                facet_grid(type ~ .,drop=FALSE) +
                labs(fill="Variant type",x="Variant size",y="Count",title=paste("Variants",comma_format(min_var),"to", comma_format(max_var),"bp")) +
                        scale_x_continuous(labels=comma_format,expand=c(0,0), limits=c(min_var-1,max_var)) +
                        scale_y_continuous(labels=comma_format,expand=c(0,0)) +
                theme(
                    strip.text=element_blank(),strip.background=element_blank(),
                    plot.title = element_text(vjust=3),
                    axis.text=element_text(size=8),
                    panel.grid.minor = element_line(colour = NA),
                    panel.grid.major = element_line(colour = NA)
                )
            )

            dev.off()
        } else {
            print("No variants in plot:")
            print(paste("min_var=",min_var))
            print(paste("max_var=",max_var))
        }

    }

}

# Prep data for log-scaled plot
alt <- bed
# Next part is required to subset rows for the Type column
alt <- alt[!is.na(alt$type),]

alt$Type <- "None"
if (nrow(alt[alt$type %in% c("Insertion","Deletion"),]) > 0) {
    alt[alt$type=="Deletion",]$size <- -1*alt[alt$type=="Deletion",]$size
    alt[alt$type %in% c("Insertion","Deletion"),]$Type <- "Indel"
}
if (nrow(alt[alt$type %in% c("Tandem expansion","Tandem contraction"),]) > 0) {
    alt[alt$type=="Tandem contraction",]$size <- -1*alt[alt$type=="Tandem contraction",]$size
    alt[alt$type %in% c("Tandem expansion","Tandem contraction"),]$Type <- "Tandem"
}
if (nrow(alt[alt$type %in% c("Repeat expansion","Repeat contraction"),]) > 0) {
    alt[alt$type=="Repeat contraction",]$size <- -1*alt[alt$type=="Repeat contraction",]$size
    alt[alt$type %in% c("Repeat expansion","Repeat contraction"),]$Type <- "Repeat"
}

if(opt$verbose){
  print("Now printing log scale plot")
}
# Save log-scaled plot to PNG
if (opt$to_png) {
  if (opt$verbose){ print("Printing large png...")}
    png(paste(opt$output_prefix,".log_all_sizes.png", sep=""),width=2000,height=1000,res=200)
} else {
  if (opt$verbose){ print("Printing large pdf...")}
    pdf(paste(opt$output_prefix,".log_all_sizes.pdf", sep=""))
}

print(ggplot(alt,aes(x=size, fill=type,y=..count..+1)) +
    geom_histogram(binwidth=opt$max_align/100, position="identity",alpha=0.7) +
    scale_fill_manual(values=big_palette,drop=FALSE) +
    facet_grid(Type ~ .,drop=FALSE) +
    labs(fill="Variant type",x="Variant size",y="Log(count + 1)",title=paste("Variants",comma_format(opt$min_query_aln),"to", comma_format(opt$max_align),"bp")) +
    scale_x_continuous(labels=comma_format,expand=c(0,0),limits=c(-1*opt$max_align,opt$max_align)) +
    scale_y_log10(labels=comma_format,expand=c(0,0)) +
    annotation_logticks(sides="l") +
    theme(
        strip.text=element_blank(),strip.background=element_blank(),
        plot.title = element_text(vjust=3),
        axis.text=element_text(size=8),
        panel.grid.minor = element_line(colour = NA),
        panel.grid.major = element_line(colour = NA)
    )
)
dev.off()
