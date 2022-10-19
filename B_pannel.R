# Modified version that only considers chrY
setwd("~/Desktop/Link2/Analyses/mapQ60/")

load("RData/02_Methylation_data_raw_alt2C_removed.RData") # Loads cov and meth data frames

library(preprocessCore) # Normalisation
library(tidyr)
library(ggplot2)
library(fmsb) # For doing the radar plot

cols <- read.csv("../allmapQ/Samples_cols2.csv", header = T, sep = ",")
seqClasses <- read.table("../../seqClass.hg38.bed", col.names = c("Chr","Start","End","Type","Length"))
CpGannots <- read.table("../../CpG_annotatr.hg38.bed", header = T)  
CpGannots <- CpGannots[CpGannots$chr == "chrY",]

split_by_seqClass <- function(seqClass) {
  
  # Normalise frequency data
  meth4_qnorm_long <- meth4[,1:9]
  meth4_qnorm_long[,3:ncol(meth4_qnorm_long)] <- normalize.quantiles(as.matrix(meth4[,3:ncol(meth4)]))
  colnames(meth4_qnorm_long)[3:ncol(meth4_qnorm_long)] <- cols$Haplogroup[sapply(3:ncol(meth4_qnorm_long), function(x) which(cols$Cell_line == names(meth4_qnorm_long)[x]))]
  
  cov4_long <- cov4[,1:9]
  
  colnames(cov4_long) <- colnames(meth4_qnorm_long)
  
  coords <- which(seqClasses$Type == seqClass)
  positions <- c()
  for (i in 1:length(coords)) {
    chr <- seqClasses$Chr[coords[i]]
    start <- seqClasses$Start[coords[i]]
    end <- seqClasses$End[coords[i]]
    positions <- c(positions, which(meth4_qnorm_long$chromosome == chr & meth4_qnorm_long$start >= start & meth4_qnorm_long$start <= end))
  }
  cov4_SeqClass <- cov4_long[positions,]
  cov4_SeqClass_g <- gather(cov4_SeqClass, haplogroup, coverage, 3:9)
  cov4_SeqClass_g_m <- merge(x = cov4_SeqClass_g, y = cols, by.x = "haplogroup", by.y = "Haplogroup")
  cov4_SeqClass_g_m$Type <- seqClass
  meth4_SeqClass <- meth4_qnorm_long[positions,]
  meth4_SeqClass_g <- gather(meth4_SeqClass, haplogroup, frequency, 3:9)
  meth4_SeqClass_g_m <- merge(x = meth4_SeqClass_g, y = cols, by.x = "haplogroup", by.y = "Haplogroup")
  meth4_SeqClass_g_m$Type <- seqClass
  
  cov4_SeqClass_g_m <- cov4_SeqClass_g_m |>
    drop_na()
  
  meth4_SeqClass_g_m <- meth4_SeqClass_g_m |>
    drop_na()
  
  return(list(cov = cov4_SeqClass_g_m, meth=meth4_SeqClass_g_m))
}

pseudo <- split_by_seqClass(seqClass = "Pseudo-autosomal")
xdeg <- split_by_seqClass(seqClass = "X-degenerate")
xtrans <- split_by_seqClass(seqClass = "X-transposed")
ampliconic <- split_by_seqClass(seqClass = "Ampliconic")
others <- split_by_seqClass(seqClass = "Others")
heterochom <- split_by_seqClass(seqClass = "Heterochromatic")

merge_seqClasses_perCpGannot <- function(CpGannot, operation) {
  
  # Load function that calculates the values for the radarplot
  
  data2CpGannot2radarplot <- function(seqClass, name_data, operation) {
    
    CpGannots_CpGannot <- CpGannots[CpGannots$type == CpGannot,]
    
    coords_CpGannot <- c()
    for (i in 1:nrow(CpGannots_CpGannot)) {
      coords_CpGannot <- c(coords_CpGannot, which(name_data$meth$start >= CpGannots_CpGannot$start[i] & name_data$meth$start <= CpGannots_CpGannot$end[i])) # To consider chrY only
    }
    
    name_data_filt <- name_data
    name_data_filt$meth <- name_data$meth[coords_CpGannot,]
    name_data_filt$cov <- name_data$cov[coords_CpGannot,]
    
    if (operation == "mean") {
      a <- aggregate(.~haplogroup, name_data_filt[[2]][,c(1,4)], mean)
    } else if (operation == "median") {
      a <- aggregate(.~haplogroup, name_data_filt[[2]][,c(1,4)], median)
    }
    
    b <- data.frame(seqClass=rep(NA, 7), row.names = unique(name_data$cov$haplogroup))
    colnames(b) <- seqClass
    b[which(rownames(b) %in% a$haplogroup),] <- data.frame(a[,2])
    return(b)
  }
  
  # Now apply it to each seqClass
  data2radarplotFinal <- cbind(data2CpGannot2radarplot(seqClass = "Heterochromatic", name_data = heterochom, operation = operation),
                               data2CpGannot2radarplot(seqClass = "Pseudo-autosomal", name_data = pseudo, operation = operation),
                               data2CpGannot2radarplot(seqClass = "X-degenerate", name_data = xdeg, operation = operation),
                               data2CpGannot2radarplot(seqClass = "X-transposed", name_data = xtrans, operation = operation),
                               data2CpGannot2radarplot(seqClass = "Ampliconic", name_data = ampliconic, operation = operation),
                               data2CpGannot2radarplot(seqClass = "Others", name_data = others, operation = operation))
  
  min_max <- matrix(nrow = 2, ncol = 6, byrow = TRUE, data = c(rep(1,6), rep(0,6)))
  colnames(min_max) <- colnames(data2radarplotFinal)
  rownames(min_max) <- c("Max", "Min")
  
  data2radarplotFinal <- rbind(min_max, data2radarplotFinal)
  
  return(data2radarplotFinal)
}

CGI_median <- merge_seqClasses_perCpGannot(CpGannot = "CGI", operation = "median")
Shore_median <- merge_seqClasses_perCpGannot(CpGannot = "CpGshores", operation = "median")
Shelf_median <- merge_seqClasses_perCpGannot(CpGannot = "CpGshelf", operation = "median")
Sea_median <- merge_seqClasses_perCpGannot(CpGannot = "openSea", operation = "median")

legend(x = 1.35, y = 0.65,
       horiz = F,
       legend = cols$Haplogroup,
       bty = "n", pch = 20, col = cols$Colour,
       text.col = "grey25", pt.cex = 1.5, cex = 1.1)


#pdf("~/Desktop/Link2/Analyses/mapQ60/Figure_main/B_pannel.pdf", width = 3.5, height = 3)
library(gridGraphics)
library(grid)
par(mar = c(0, 0, 0.75, 0))

radarchart(CGI_median[,c(3:5,2,1,6)], # data
           vlcex = 1, # Size for labels
           axistype=1, # axistype=1 means center axis label only.
           seg = 4, # Number of segments for each axis
           # Custom the points of the plot (p*)
           pty = 19, pcol = cols$Colour, plwd = 1.5, plty = 1,
           # Custom the grid (cgl*)
           cglcol="grey50", cglty=2, cglwd=0.8,
           na.itp = F,
           # Axis labels 
           caxislabels=c("0 ", "0.25", "0.5 ", "0.75", "1 "), axislabcol="grey25", calcex = 1,
           title = "CGI")


radarchart(Shore_median[,c(3:5,2,1,6)], # data
           vlcex = 1, # Size for labels
           axistype=1, # axistype=1 means center axis label only.
           seg = 4, # Number of segments for each axis
           # Custom the points of the plot (p*)
           pty = 19, pcol = cols$Colour, plwd = 1.5, plty = 1,
           # Custom the grid (cgl*)
           cglcol="grey50", cglty=2, cglwd=0.8,
           na.itp = F,
           # Axis labels 
           caxislabels=c("0 ", "0.25", "0.5 ", "0.75", "1 "), axislabcol="grey25", calcex = 1,
           title = "CpG shore")

radarchart(Shelf_median[,c(3:5,2,1,6)], # data
           vlcex = 1, # Size for labels
           axistype=1, # axistype=1 means center axis label only.
           seg = 4, # Number of segments for each axis
           # Custom the points of the plot (p*)
           pty = 19, pcol = cols$Colour, plwd = 1.5, plty = 1,
           # Custom the grid (cgl*)
           cglcol="grey50", cglty=2, cglwd=0.8,
           na.itp = F,
           # Axis labels 
           caxislabels=c("0 ", "0.25", "0.5 ", "0.75", "1 "), axislabcol="grey25", calcex = 1,
           title = "CpG shelf")

radarchart(Sea_median[,c(3:5,2,1,6)], # data
           vlcex = 1, # Size for labels
           axistype=1, # axistype=1 means center axis label only.
           seg = 4, # Number of segments for each axis
           # Custom the points of the plot (p*)
           pty = 19, pcol = cols$Colour, plwd = 1.5, plty = 1,
           # Custom the grid (cgl*)
           cglcol="grey50", cglty=2, cglwd=0.8,
           na.itp = F,
           # Axis labels 
           caxislabels=c("0 ", "0.25", "0.5 ", "0.75", "1 "), axislabcol="grey25", calcex = 1,
           title = "Open sea")
#dev.off()
