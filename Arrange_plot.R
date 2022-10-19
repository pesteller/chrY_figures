# Arrange all plots in one

library(ggpubr)

pannelC <- ggarrange(pannelC1, pannelC2, nrow = 1)

mid <- ggarrange(pannelD, pannelD, nrow = 1, widths = c(2,1))

ggarrange(mid, pannelC, nrow = 2, align = "v")

ggarrange(pannelA, mid, pannelC, nrow = 3, align = "v", heights = c(3,2,2))

library(gridExtra)
library(grid)
library(ggplot2)

grid.arrange(pannelA, pannelD, pannelD, pannelC, nrow=3)

lay <- rbind(c(1,1,1),
             c(1,1,1),
             c(1,1,1),
             c(2,2,3),
             c(2,2,3),
             c(4,4,4))

grid.arrange(grobs = list(pannelA, pannelD, pannelD, pannelC), 
                  layout_matrix = lay)

grid.arrange(grobs =  pl[1:5], layout_matrix = lay)

gs <- lapply(1:9, function(ii) 
  grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))

pl <- lapply(1:11, function(.x) qplot(1:10, rnorm(10), main=paste("plot", .x)))

