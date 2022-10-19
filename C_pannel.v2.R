# Modified version that only considers chrY
setwd("~/Desktop/Link2/Analyses/mapQ60/")

plot_data_per_seqClass <- function(seqClass) {
  
  library(dplyr)
  require(ggplot2) # Plotting
  require(ggridges)

  # Load function to retrieve data per seqClass and gene annotation
  get_data_seqClass_gene_annot <- function(seqClass, gene_annot) {
    
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
      name <- "Intragenic"
    } else if (gene_annot == "TSS_200up") {
      Geneannots <- read.table("../../01.5_TSS200.hg38.bed", header = F, col.names = c("chr", "start", "end", "strand", "type"))
      Geneannots$type <- "TSS_200up"
      name <- "TSS"
    } else if (gene_annot == "5UTR") {
      Geneannots <- read.table("../../GeneAnnotation_chrY.ranked.multiinter.bothUTRs.hg38.bed", header = F, col.names = c("chr", "start", "end", "length", "type"))  
      name <- "5'UTR"
    } else if (gene_annot == "3UTR") {
      Geneannots <- read.table("../../GeneAnnotation_chrY.ranked.multiinter.bothUTRs.hg38.bed", header = F, col.names = c("chr", "start", "end", "length", "type"))  
      name <- "3'UTR"
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
    
    meth4_qnorm_long_seqClass_GeneAnnot_g$gene_annot <- name
    
    return(meth4_qnorm_long_seqClass_GeneAnnot_g)
  }
  
  # Apply it to our data
  dat_TSS <- get_data_seqClass_gene_annot(seqClass = seqClass, gene_annot = "TSS_200up")
  dat_5UTR <- get_data_seqClass_gene_annot(seqClass = seqClass, gene_annot = "5UTR")
  dat_intragenic <- get_data_seqClass_gene_annot(seqClass = seqClass, gene_annot = "Intragenic")
  dat_3UTR <- get_data_seqClass_gene_annot(seqClass = seqClass, gene_annot = "3UTR")
  
  cols <- read.csv("../allmapQ/Samples_cols2.csv", header = T, sep = ",")
  
  plt <- rbind(dat_TSS, dat_5UTR, dat_intragenic, dat_3UTR) |>
    mutate(gene_annot = factor(gene_annot, levels = c("TSS", "5'UTR", "Intragenic", "3'UTR"))) |>
    ggplot(aes(y = Sample, x = Freq_meth, group = Sample, fill = Sample)) +
    geom_density_ridges2(show.legend = FALSE) +
    facet_grid(~gene_annot) +
    stat_density_ridges(quantile_lines = TRUE, quantiles = 2, show.legend = FALSE, size = 0.35) +
    scale_fill_manual(breaks = cols$Haplogroup, values = cols$Colour) +
    scale_y_discrete(expand = expansion(0.05)) +
    scale_x_continuous(limits = c(-0.05,1.05), labels = c("0","0.25","0.5","0.75","1"), breaks = c(0,0.25,0.5,0.75,1)) +
    labs(x = "5mC frequency", y = NULL)+
    ggtitle(seqClass) +
#    lims(x = c(-0.05,1.05))+
    theme_minimal(base_size = 8) +
    theme(strip.text.x = element_text(size = 8),
          panel.border = element_rect(size = 0.75, fill = NA),
          axis.ticks = element_line(colour = "black"),
          axis.ticks.length = unit(.1, "cm"),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank())
  
  return(plt)
}

pannelC1 <- plot_data_per_seqClass(seqClass = "X-degenerate")
pannelC2 <- plot_data_per_seqClass(seqClass = "Ampliconic")
