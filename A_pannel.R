setwd("~/Desktop/Link2/Analyses/mapQ60/")
library(cowplot)
library(ggpubr)

plt_tiles_points <- function(window_size, min_number, operation, plt, lower_lim, upper_lim) {
  
  # window_size: window size (in bp)
  # min_number: minimum number of positions with methylation information to keep the window 
  # plt (options: 1,2,both): to print the first(tiles), second(points) or both plots
  # lower_lim: lower limit to use for plotting (and also for doing the windows, use 1 for 1Mbp)
  # upper_lim: upper limit to use for plotting (and also for doing the windows, use 1 for 1Mbp)
  
  require(ggplot2)
  require(tidyr)
  require(scales)
  require(preprocessCore) # Normalisation
  require(RColorBrewer)
  
  seqClasses <- read.table("../../seqClass.hg38.bed", col.names = c("Chr","Start","End","Type","Length"))
  
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
  
  if (plt != 1 & plt != 2 & plt != "both") {
    stop("plt only accepts 1, 2 or \"both\".")
  }
  if (operation != "mean" & operation != "median") {
    stop("operation can only be \"mean\" or \"median\".")
  }
  
  # 1. SeqClasses 
  
  # 1.1. Reformat seqClasses for plotting
  seqClasses_filt <- setLims(lower_lim = lower_lim, upper_lim = upper_lim) 
  
  cols_seq <- c("#F8E4A0", "#EDBFE1", "#F68E79", "#A1E3BB", "#C8CFE4", "grey85")
  names(cols_seq) <- c("Ampliconic", "Heterochromatic", "Others", "X-degenerate", "X-transposed", "Pseudo-autosomal")
  
  # 1.2. Generate the plot
  plt_seq <- ggplot(seqClasses_filt, aes(xmin = Start, xmax = End, ymin = 0, ymax = 1, fill = Type)) +
    geom_rect(stat = "identity",position = "identity") +
    guides(fill = guide_legend(nrow = 1), ) +
    scale_fill_manual(values = cols_seq[which(names(cols_seq) %in% unique(seqClasses_filt$Type))]) +
    scale_x_continuous(limits = c(lower_lim-0.5,upper_lim), expand = expansion(0)) +
    labs(y = "Sequence classes", fill = "Sequence class") +
    #  coord_cartesian(expand = FALSE) +
    theme_void() +
    theme(legend.position = "top",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          legend.justification = c(0,0),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          axis.title.y = element_text(colour = "black", size = 8, hjust = 0)) 
  
  # 2. Plot for genes
  genes <- read.table("../../hg38.ensembl.chrY.gene.coords.bed", col.names = c("Chr", "Start", "End", "Strand", "Gene_id", "Gene_name", "Gene_type"))
  
  genes$type <- NA
  
  for (i in 1:nrow(genes)) {
    if (genes$Gene_type[i]  == "unprocessed_pseudogene" | genes$Gene_type[i]  == "processed_pseudogene" | 
        genes$Gene_type[i]  == "transcribed_unprocessed_pseudogene" | genes$Gene_type[i]  == "transcribed_processed_pseudogene") {
      genes$type[i] <- "Pseudogenes"
    } else if (genes$Gene_type[i] == "lincRNA" | genes$Gene_type[i] == "snRNA" | genes$Gene_type[i] == "rRNA" | genes$Gene_type[i] == "misc_RNA" |
               genes$Gene_type[i] == "antisense_RNA" | genes$Gene_type[i] == "snoRNA" | genes$Gene_type[i] == "processed_transcript" | genes$Gene_type[i] == "sense_intronic") {
      genes$type[i] <- "Non-coding"
    } else if (genes$Gene_type[i] == "protein_coding") {
      genes$type[i] <- "Protein coding"
    } else {
      stop(paste0("Gene type not classified for ", i, " line."))
    }
  }
  
  # Plot of genes coloured by type of gene
  plt_genes <- ggplot() + 
    geom_rect(data = genes[genes$Strand == "+",], aes(xmin = Start/1000000, xmax = End/1000000, ymin = 1, ymax = 2, fill = type)) +
    geom_rect(data = genes[genes$Strand == "-",], aes(xmin = Start/1000000, xmax = End/1000000, ymin = 0, ymax = 1, fill = type)) +
    scale_fill_manual(values = c("seagreen4", "plum4", "sienna3"), breaks = c("Protein coding", "Pseudogenes", "Non-coding")) +
    scale_y_continuous(breaks = c(0.5, 1.5), labels = c("-", "+")) +
    scale_x_continuous(limits = c(lower_lim-0.5,upper_lim), expand = expansion(0)) +
    labs(y = "Genes", fill = "Gene type") +
    theme_void()+
    theme(legend.position = "bottom",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          legend.justification = c(0,0),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 7),
          axis.text.y = element_text(colour = "black", size = 8),
          axis.title.y = element_text(colour = "black", size = 8, hjust = 0))
          
  # 2. Calculate windows
  
  load("RData/02_Methylation_data_raw_alt2C_removed.RData") # Loads cov and meth data frames
  cols <- read.csv("../allmapQ/Samples_cols2.csv", header = T, sep = ",")
  
  #max <- max(cov4$start) # From previous version (calculates the max number of methylation)
  
  max <- upper_lim*1000000 # Go from bp to Mpb
  
  window = window_size 
  sliding = window/2 
  
  # Regardless of the boundary you want to plot, it always starts computing the windows from the start of the chromosome
  #start <- seq(0,max, by = sliding)
  #end <- seq(window,max, by = sliding)
  
  # Start computing the windows from the lower limit given
  start <- seq((lower_lim)*1000000,max, by = sliding)
  end <- seq((lower_lim)*1000000+window,max, by = sliding)
  
  length(start)
  length(end)
  
  dat <- data.frame(matrix(nrow = length(start), ncol = 11))
  colnames(dat) <- c("Chrom", "Start", "End",colnames(cov4)[3:9],"mad")
  dat$Chrom <- "chrY"
  dat$Start <- start
  dat$End <- c(end, max, max)
  dat$mad <- NA
  dat$sd <- NA
  
  meth4 <- meth4[1:9] # chr, start, 7 cell lines
  meth4[,3:ncol(meth4)] <- normalize.quantiles(as.matrix(meth4[,3:ncol(meth4)])) # Normalise all the samples but H_AS and H_CS
  
  dat_counts <- dat
  for (i in 1:nrow(dat)) {
    coords <- which(meth4$chromosome == "chrY" & meth4$start >= dat$Start[i] & meth4$start < dat$End[i]) # Consider chrY only
    for (j in 4:(4+7)) {
      if (sum(!is.na(meth4[coords,j-1])) >= min_number) {
        if (operation == "median") {
          dat[i,j] <- median(meth4[coords,j-1], na.rm = T)
        } else if (operation == "mean") {
          dat[i,j] <- mean(meth4[coords,j-1], na.rm = T)
        }
      }
      #    dat_counts[i,j] <- sum(!is.na(meth4[coords,j-1]))
    }
    if(sum(!is.na(dat[i,4:10])) >= 3) {
      dat$mad[i] <- mad(dat[i,4:10], na.rm = T)
      dat$sd[i] <- sd(dat[i,4:10], na.rm = T)
      
    }
  }
  
  dat_g <- gather(dat, key = "Sample", value = "Freq_meth", 4:10)
  
  dat_gm <- merge(x = dat_g, y = cols, by.x = "Sample", by.y = "Cell_line")
  
  # Geom mad: this is a tile plot that show the dispersion per window
  plt_mad <- ggplot(dat, aes(x = Start/1000000, y = 0, fill = mad)) +
    geom_tile(lwd = 0.25, linetype = 1, show.legend = F) +
    scale_fill_gradientn(limits=c(0,max(dat$mad, na.rm = T)),
                         colors = brewer.pal(9,"Blues"), na.value="grey75") +
    lims(x = c(lower_lim,upper_lim)) +
    theme_void()
  
  # Geom sd: this is a tile plot that show the dispersion per window
  plt_sd <- ggplot(dat, aes(x = Start/1000000, y = 0, fill = sd)) +
    geom_tile(lwd = 0.2, linetype = 1, show.legend = F) +
    scale_fill_gradientn(limits=c(0,max(dat$sd, na.rm = T)),
                         colors = brewer.pal(9,"Greys"), na.value="white") +
                  #       colors = brewer.pal(9,"Blues"), na.value="grey75") +
    scale_x_continuous(limits = c(lower_lim-0.5,upper_lim), expand = expansion(0)) +
    labs(y = "5mC freq. sd") +
    theme_void() +
    theme(axis.title.y = element_text(colour = "black", size = 8, hjust = 0, vjust = 0.5))

  # Geom points
  plt_points <- ggplot(dat_gm, aes(x = Start/1000000, y = Freq_meth, col = Haplogroup)) +
    geom_hline(yintercept = c(0.5,1), lty = 2, col = "grey80") +
    geom_point(size = 0.75) +
    geom_line(size = 0.25) +
    scale_y_continuous(limits = c(0,1), labels = c("0","0.25","0.5","0.75","1"), breaks = c(0,0.25,0.5,0.75,1)) +
    scale_x_continuous(breaks = c(5,10,15,20,25), labels = c(5,10,15,20,25), limits = c(lower_lim-0.5,upper_lim), expand = expansion(0)) +
    scale_color_manual(breaks = cols$Haplogroup, values = cols$Colour) +
    labs(y = "5mC frequency", x = "Position (Mb)") +#,
    guides(colour = guide_legend(nrow = 1)) +
    theme_bw(base_size = 12) +
    theme(panel.border = element_rect(size = 0.75, fill = NA),
          panel.grid = element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          axis.text = element_text(colour = "black"),
          legend.position = "bottom",
          legend.text = element_text(colour = "black", size = 9),
          legend.justification = c(0,0),
          legend.title = element_text(size = 8))
  
  # Heat map with geom_tile
  plt_tiles <- ggplot(dat_gm, aes(x = Start/1000000, y = Haplogroup, fill = Freq_meth)) +
    geom_tile(lwd = 0.25, linetype = 1) +
    labs(x = "Position (Mb)")+ 
    #      title = paste0("Windows size: ", window_size, " bp - Min counts: ", min_number),
    #       subtitle = paste0("Operation: ", operation)) +
    scale_fill_gradientn(limits=c(0,1),
                         colors = brewer.pal(9,"RdBu"),
                         values = rescale(c(0,0.5,1)), na.value="grey75") +
    #coord_fixed() +
    #  coord_cartesian(expand = FALSE) +
    lims(x = c(lower_lim,upper_lim)) +
    guides(fill = guide_colourbar(barwidth = 8,
                                  barheight = 0.5,
                                  title = "5mC frequency",
                                  title.vjust = 1 )) +
    theme_classic(base_size = 15, base_line_size = 0.25) +
    theme(
      title = element_text(size = 10),
      axis.text = element_text(colour = "black"),
      legend.direction="horizontal",
      legend.text = element_text(size = 8),
      legend.position = "bottom",
      legend.justification=c(0.4,0),
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,0,0,0)) 
  
  
  
  if (plt == 1) {
    return(list(tiles = plt_tiles, seq = plt_seq, mad = plt_mad, sd = plt_sd, genes = plt_genes))
  } else if (plt == 2) {
    return(list(points = plt_points, seq = plt_seq, mad = plt_mad, sd = plt_sd, genes = plt_genes))
  } else if (plt == "both") {
    return(list(tiles = plt_tiles, points = plt_points, seq = plt_seq, mad = plt_mad, sd = plt_sd, genes = plt_genes))
  }
  
}

# Do plot of the resolved MSY of the hg38: 2781333 bp to 27078487 bp
plt <- plt_tiles_points(window_size = 250000, min_number = 10, operation = "median", plt = 2, lower_lim = 2781333/1000000, upper_lim = 27078487/1000000)

plt_seq_empty <- plt$seq + theme(legend.position = "none")
legend_seq <- as_ggplot(get_legend(plt$seq))

plt_genes_empty <- plt$genes + theme(legend.position = "none")
legend_genes <- as_ggplot(get_legend(plt$genes))

plt_points_empty <- plt$points + theme(legend.position = "none")
legend_points <- as_ggplot(get_legend(plt$points))

ggarrange(legend_seq, legend_genes, legend_points, plt_seq_empty, plt_genes_empty, plt$sd, plt_points_empty, align = "v", ncol = 1, heights = c(0.3,0.3,0.3,0.4,0.8,0.4,5))

#pannelA <- ggarrange(legend_seq, legend_genes, legend_points, plt_seq_empty, plt_genes_empty, plt$sd, plt_points_empty, align = "v", ncol = 1, heights = c(0.4,0.4,0.4,0.5,1,0.5,5))

# Title not needed
#title <- ggdraw() + draw_label("Resolved MSY of  hg38")

