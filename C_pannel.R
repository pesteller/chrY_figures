# Modified version that only considers chrY

setwd("~/Desktop/Link2/Analyses/mapQ60/")

plot_by_seqClass_CpGannot <- function(seqClass, gene_annot) {

  require(preprocessCore) # Normalisation
  require(ggplot2) # Plotting
  require(ggridges)
  require(tidyr)
  
  load("RData/02_Methylation_data_raw_alt2C_removed.RData") # Loads cov and meth data frames
  
  # Import the colors and seqClasses
  cols <- read.csv("../allmapQ/Samples_cols2.csv", header = T, sep = ",")
  seqClasses <- read.table("../../seqClass.hg38.bed", col.names = c("Chr","Start","End","Type","Length"))

  if (gene_annot == "Intragenic") {
    Geneannots <- read.table("../../GeneAnnotation_chrY.ranked.multiinter.bothUTRs.hg38.bed", header = F, col.names = c("chr", "start", "end", "length", "type"))  
    Geneannots$type[which(Geneannots$type %in% c("Intron", "Exon"))] <- "Intragenic"
    title <- "Intragenic"
  } else if (gene_annot == "TSS_200up") {
    Geneannots <- read.table("../../01.5_TSS200.hg38.bed", header = F, col.names = c("chr", "start", "end", "strand", "type"))
    Geneannots$type <- "TSS_200up"
    title <- "TSS (200 bp upstream)"
  } else if (gene_annot == "5UTR") {
    Geneannots <- read.table("../../GeneAnnotation_chrY.ranked.multiinter.bothUTRs.hg38.bed", header = F, col.names = c("chr", "start", "end", "length", "type"))  
    title <- "5'UTR"
  } else if (gene_annot == "3UTR") {
    Geneannots <- read.table("../../GeneAnnotation_chrY.ranked.multiinter.bothUTRs.hg38.bed", header = F, col.names = c("chr", "start", "end", "length", "type"))  
    title <- "3'UTR"
  }
  
  # Normalise frequency data
  meth4_qnorm_long <- meth4[,1:9]
  meth4_qnorm_long[,3:ncol(meth4_qnorm_long)] <- normalize.quantiles(as.matrix(meth4[,3:ncol(meth4)]))
  colnames(meth4_qnorm_long)[3:ncol(meth4_qnorm_long)] <- cols$Haplogroup[sapply(3:ncol(meth4_qnorm_long), function(x) which(cols$Cell_line == names(meth4_qnorm_long)[x]))]
  
  start <- seqClasses$Start[which(seqClasses$Type == seqClass)]
  end <- seqClasses$End[which(seqClasses$Type == seqClass)]
  
  coords_seqClass <- c()
  for (i in 1:length(start)) {
    coords_seqClass <- c(coords_seqClass, which(meth4_qnorm_long$chromosome == "chrY" & meth4_qnorm_long$start >= start[i] & meth4_qnorm_long$start <= end[i])) # Consider chrY only
  }
  
  meth4_qnorm_long_seqClass <- meth4_qnorm_long[coords_seqClass,]
  
  GeneAnnots_GeneAnnot <- Geneannots[Geneannots$type == gene_annot & Geneannots$chr == "chrY",]
  
  coords_GeneAnnot <- c()
  for (i in 1:nrow(GeneAnnots_GeneAnnot)) {
    coords_GeneAnnot <- c(coords_GeneAnnot, which(meth4_qnorm_long_seqClass$start >= GeneAnnots_GeneAnnot$start[i] & meth4_qnorm_long_seqClass$start <= GeneAnnots_GeneAnnot$end[i]))
  }
  
  meth4_qnorm_long_seqClass_GeneAnnot <- meth4_qnorm_long_seqClass[coords_GeneAnnot,]
  
  meth4_qnorm_long_seqClass_GeneAnnot_g <- pivot_longer(meth4_qnorm_long_seqClass_GeneAnnot, cols = 3:9, names_to = "Sample", values_to = "Freq_meth")
  
  plt_dens_ridge <- ggplot(meth4_qnorm_long_seqClass_GeneAnnot_g, aes(y = Sample, x = Freq_meth, group = Sample, fill = Sample)) +
    geom_density_ridges2(show.legend = FALSE) +
    stat_density_ridges(quantile_lines = TRUE, quantiles = 2, show.legend = FALSE) +
    scale_fill_manual(breaks = cols$Haplogroup, values = cols$Colour) +
    scale_y_discrete(expand = expansion(0.05)) +
    labs(x = "5mC frequency", y = NULL)+
    ggtitle(title) +
    lims(x = c(0,1))+
    theme_bw(base_size = 8) +
    theme(panel.border = element_rect(size = 0.5, fill = NA),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank())
  
  plt_dens <- ggplot(meth4_qnorm_long_seqClass_GeneAnnot_g, aes(x = Freq_meth, group = Sample, col = Sample)) +
    geom_density(show.legend = FALSE, size = 0.75) +
    scale_color_manual(breaks = cols$Haplogroup, values = cols$Colour) +
    labs(x = "5mC frequency", y = "Density")+
    ggtitle(title) +    
    scale_x_continuous(limits = c(0-0.1,1.1), expand = expansion(0)) +
    theme_bw(base_size = 8) +
    theme(panel.border = element_rect(size = 0.5, fill = NA),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank())

  if (length(which(is.na(meth4_qnorm_long_seqClass_GeneAnnot_g$Freq_meth)) > 0)) {
    dat <- meth4_qnorm_long_seqClass_GeneAnnot_g[-which(is.na(meth4_qnorm_long_seqClass_GeneAnnot_g$Freq_meth)),]
  } else {
    dat <- meth4_qnorm_long_seqClass_GeneAnnot_g
  }
  
  plt_count <- dat |>
    na.omit() |>
    ggplot(aes(y = Sample, fill = Sample)) +
    geom_bar(col = "black", show.legend = FALSE) +
    stat_count(geom = "text", colour = "black", size = 4,
               aes(label = ..count..),
               position=position_stack(vjust=0.5)) +
    scale_fill_manual(breaks = cols$Haplogroup, values = cols$Colour) +
    labs(x = "Count", y = NULL)+
    ggtitle(title) +
    theme_bw(base_size = 8) +
    theme(panel.border = element_rect(size = 0.5, fill = NA),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank())
  
  return(list(ridge = plt_dens_ridge, 
              dens = plt_dens, 
              count = plt_count, 
              counts = table(dat$Sample),
              dat = meth4_qnorm_long_seqClass_GeneAnnot_g))
  
}

Xdeg_TSS <- plot_by_seqClass_CpGannot(seqClass = "X-degenerate", gene_annot = "TSS_200up")
Xdeg_5UTR <- plot_by_seqClass_CpGannot(seqClass = "X-degenerate", gene_annot = "5UTR")
Xdeg_intragenic <- plot_by_seqClass_CpGannot(seqClass = "X-degenerate", gene_annot = "Intragenic")
Xdeg_3UTR <- plot_by_seqClass_CpGannot(seqClass = "X-degenerate", gene_annot = "3UTR")

Amp_TSS <- plot_by_seqClass_CpGannot(seqClass = "Ampliconic", gene_annot = "TSS_200up")
Amp_5UTR <- plot_by_seqClass_CpGannot(seqClass = "Ampliconic", gene_annot = "5UTR")
Amp_intragenic <- plot_by_seqClass_CpGannot(seqClass = "Ampliconic", gene_annot = "Intragenic")
Amp_3UTR <- plot_by_seqClass_CpGannot(seqClass = "Ampliconic", gene_annot = "3UTR")

library(ggpubr)
plt_Xdeg <- ggarrange(Xdeg_TSS$ridge, Xdeg_5UTR$ridge, Xdeg_intragenic$ridge, Xdeg_3UTR$ridge, ncol = 4)
plt_Xdeg <- annotate_figure(plt_Xdeg, top = text_grob("X-degenerate", size = 12))

plt_Amp <- ggarrange(Amp_TSS$ridge, Amp_5UTR$ridge, Amp_intragenic$ridge, Amp_3UTR$ridge, ncol = 4)
plt_Amp <- annotate_figure(plt_Amp, top = text_grob("Ampliconic", size = 12))

plt_counts_Xdeg <- ggarrange(Xdeg_TSS$count, Xdeg_5UTR$count, Xdeg_intragenic$count, Xdeg_3UTR$count, ncol = 4)
plt_counts_Xdeg <- annotate_figure(plt_counts_Xdeg, top = text_grob("X-degenerate", size = 12))

plt_counts_Amp <- ggarrange(Amp_TSS$count, Amp_5UTR$count, Amp_intragenic$count, Amp_3UTR$count, ncol = 4)
plt_counts_Amp <- annotate_figure(plt_counts_Amp, top = text_grob("Ampliconic", size = 12))

ggarrange(plt_Xdeg, plt_Amp, nrow = 1)
