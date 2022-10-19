setwd("~/Desktop/Link2/Analyses/mapQ60/")

# Load data where tracks are stored
load("RData/20_CpGannot_retrieve_mean_median_freq_counts.RData")

plt_tiles_points <- function(lower_lim, upper_lim, track) {
  
  # lower_lim: lower limit to use for plotting (in bp)
  # upper_lim: upper limit to use for plotting (in bp)
  # track: track generated in 20_CpGannot_retrive_mean_median_freq_counts.R and stored in .RData
  
  require(ggplot2)
  require(tidyr)
  require(scales)
  require(ggrepel)
  require(dplyr)
  require(preprocessCore) # Normalisation
  
  # Convert limits from bp to Mb
  lower_lim <- lower_lim/1000000
  upper_lim <- upper_lim/1000000
  
  # Load seqClasses and cols
  seqClasses <- read.table("../../seqClass.hg38.bed", col.names = c("Chr","Start","End","Type","Length"))
  cols <- read.csv("../allmapQ/Samples_cols2.csv", header = T, sep = ",")
  
  # Function to calculate limits
  setLims <- function(lower_lim, upper_lim) {
    
    seqClasses_out <- seqClasses
    seqClasses_out$Start <- seqClasses_out$Start/1000000
    seqClasses_out$End <- seqClasses_out$End/1000000
    
    seqClasses_out$Start[which(seqClasses_out$End >= lower_lim & seqClasses_out$Start <= lower_lim)] <- lower_lim
    seqClasses_out$End[which(seqClasses_out$End >= upper_lim & seqClasses_out$Start <= upper_lim)] <- upper_lim
    
    if (lower_lim != 0) {
      seqClasses_out <- seqClasses_out[which(seqClasses_out$Start == lower_lim):which(seqClasses_out$End == upper_lim),] # Remove all other rows after upper_lim
    } else {
      seqClasses_out <- seqClasses_out[1:which(seqClasses_out$End == upper_lim),] # Remove all other rows after upper_lim
    }
    
    return(seqClasses_out)
    
  }
  
  # 1. SeqClasses 
  # 1.1. Reformat seqClasses for plotting
  seqClasses_filt <- setLims(lower_lim = lower_lim, upper_lim = upper_lim) 
  
  cols_seq <- c("#F8E4A0", "#EDBFE1", "#F68E79", "#A1E3BB", "#C8CFE4", "grey85")
  names(cols_seq) <- c("Ampliconic", "Heterochromatic", "Others", "X-degenerate", "X-transposed", "Pseudo-autosomal")
  
  # 1.2. Generate the plot
  plt_seq <- ggplot(seqClasses_filt, aes(xmin = Start, xmax = End, ymin = 0, ymax = 1, fill = Type, label = Type)) +
    geom_rect(stat = "identity", position = "identity") +
    geom_text(aes(x = Start + (End - Start)/2, y = 0.5), size = 3, show.legend = F) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_fill_manual(values = cols_seq[which(names(cols_seq) %in% unique(seqClasses_filt$Type))]) +
    # lims(x = c(lower_lim,upper_lim)) +
    #  coord_cartesian(expand = FALSE) +
    theme_void() +
    theme(legend.position = "none",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0)) 
  
  # 2. Plot for genes
  genes <- read.table("~/Desktop/chrY/hg38.ensembl.chrY.gene.coords.bed", col.names = c("Chr", "Start", "End", "Strand", "Gene_id", "Gene_name", "Gene_type"))
  
  genes$type <- NA
  
  for (i in 1:nrow(genes)) {
    if (genes$Gene_type[i]  == "unprocessed_pseudogene" | genes$Gene_type[i]  == "processed_pseudogene" | 
        genes$Gene_type[i]  == "transcribed_unprocessed_pseudogene" | genes$Gene_type[i]  == "transcribed_processed_pseudogene") {
      genes$type[i] <- "Pseudogenes"
    } else if (genes$Gene_type[i] == "lincRNA" | genes$Gene_type[i] == "snRNA" | genes$Gene_type[i] == "rRNA" | genes$Gene_type[i] == "misc_RNA" |
               genes$Gene_type[i] == "antisense_RNA" | genes$Gene_type[i] == "snoRNA" | genes$Gene_type[i] == "processed_transcript" | genes$Gene_type[i] == "sense_intronic") {
      genes$type[i] <- "Non_coding"
    } else if (genes$Gene_type[i] == "protein_coding") {
      genes$type[i] <- "Protein_coding"
    } else {
      stop(paste0("Gene type not classified for ", i, " line."))
    }
  }
  
  # Plot of genes coloured by type of gene
  plt_genes <- ggplot() + 
    geom_rect(data = genes[genes$Strand == "+",], aes(xmin = Start/1000000, xmax = End/1000000, ymin = 1, ymax = 2, fill = type), show.legend = F) +
    geom_rect(data = genes[genes$Strand == "-",], aes(xmin = Start/1000000, xmax = End/1000000, ymin = 0, ymax = 1, fill = type), show.legend = F) +
    scale_fill_manual(values = c("seagreen4", "plum4", "sienna3"), breaks = c("Protein_coding", "Pseudogenes", "Non_coding")) +
    scale_y_continuous(breaks = c(0.5, 1.5), labels = c("-", "+")) +
    lims(x = c(lower_lim,upper_lim)) +
    labs(y = "Genes", fill = "Gene type") +
    theme_void()+
    theme(legend.position = "none",
          axis.text.y = element_text(colour = "black", size = 8),
          axis.title.y = element_text(colour = "black", size = 8, hjust = 0))
  
  plt_track <- ggplot(track$freq, aes(xmin = start/1000000, xmax = end/1000000, ymin = 0, ymax = 1)) + 
    geom_rect(stat = "identity", position = "identity", fill = "black") +
    geom_point(aes(x = (start + (end-start)/2)/1000000, y = 0), pch = 24, fill = "black", size = 1) +
    geom_point(aes(x = (start + (end-start)/2)/1000000, y = 1), pch = 25, fill = "black", size = 1) +
    lims(x = c(lower_lim,upper_lim), y = c(0,1)) +
    labs(y = strsplit(track$freq$type[1],"_")[[1]][1]) +
    theme_void() +
    theme(axis.title.y = element_text(colour = "black", size = 8, hjust = 0))

  plt_meth <- track$freq |>    
    pivot_longer(cols = 5:11, names_to = "Haplogroup", values_to = "Freq_meth") |>
    mutate(mean = start + (end-start)/2) |>
    filter(type %in% c("CGI_18", "CGI_19", "CGI_20")) |>
    ggplot(aes(x = mean/1000000, y = Freq_meth, fill = Haplogroup, shape = type)) +
    geom_jitter(size = 3.5, width = 0.02, height = 0.02, stroke = .35, show.legend = F) +
    #geom_jitter(pch = 21, size = 4, width = 0.02, height = 0.02, stroke = .35) +
    scale_fill_manual(values = alpha(cols$Colour, 0.8), breaks = cols$Haplogroup) +
    scale_shape_manual(values = c(21,22,21), breaks = c("CGI_18", "CGI_19", "CGI_20")) +
    guides(fill = guide_legend(nrow = 1)) +
    lims(x = c(lower_lim,upper_lim), y = c(0,1)) +
    labs(x = "Position (Mb)", y = "5mC freq (median)") +
    theme_bw(base_size = 10) +
    theme(panel.border = element_rect(size = 0.5, fill = NA),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank(),
          legend.position = "none")
    
  dat_mean <- track$freq |>    
    pivot_longer(cols = 5:11, names_to = "Haplogroup", values_to = "Freq_meth") |>
    mutate(mean = start + (end-start)/2) |>
    filter(type %in% c("CGI_18", "CGI_19", "CGI_20")) |>
    mutate(type = recode_factor(type, "CGI_18" = "CGI_1", "CGI_19" = "CGI_2", "CGI_20" = "CGI_3")) |>
    group_by(type) |>
    summarise(chr = unique(chr), 
              mean = unique(mean), 
              Freq_meth = mean(Freq_meth))
  
  plt_meth2 <- track$freq |>    
    pivot_longer(cols = 5:11, names_to = "Haplogroup", values_to = "Freq_meth") |>
    mutate(mean = start + (end-start)/2) |>
    filter(type %in% c("CGI_18", "CGI_19", "CGI_20")) |>
    ggplot() +
    geom_jitter(aes(x = mean/1000000, y = Freq_meth, col = Haplogroup), pch = 20, size = 1.5, width = 0.015, height = 0.015, stroke = .35, show.legend = F) +
    geom_point(data = dat_mean, aes(x = mean/1000000, y = Freq_meth), pch = 21, size = 3, stroke = .35, show.legend = F) +
    geom_text_repel(data = dat_mean, aes(x = mean/1000000, y = Freq_meth, label = type), size = 3) +
    scale_colour_manual(values = alpha(cols$Colour, 0.8), breaks = cols$Haplogroup) +
    scale_y_continuous(limits = c(0,1), labels = c("0","0.25","0.5","0.75","1"), breaks = c(0,0.25,0.5,0.75,1)) +
    lims(x = c(lower_lim,upper_lim)) +
    labs(x = "Position (Mb)", y = "5mC frequency") +
    theme_bw(base_size = 10) +
    theme(panel.border = element_rect(size = 0.75, fill = NA),
          axis.text = element_text(colour = "black"),
          panel.grid = element_blank())
  
  library(ggpubr)
  
  plt <- ggarrange(plt_seq, plt_genes, plt_track, plt_meth2,
                   ncol = 1, align = "v", heights = c(0.25,0.5,0.25,2.5))

  return(plt)
}

pannelD <- plt_tiles_points(lower_lim = 14.5*1e6, upper_lim = 14.9*1e6, track = Xdeg_CGI)
