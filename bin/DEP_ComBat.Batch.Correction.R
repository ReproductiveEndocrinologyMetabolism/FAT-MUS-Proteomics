#library(DEqMS)
library(matrixStats)
library(ggplot2)
library(cowplot)
library(readxl)
#library(ComBat)
#library(gmm)
library(DEP)
library(dplyr)
library(ggfortify)
library(reshape2)
library(sva)
library(data.table)
library(enrichR)

# Set default theme to cowplot
theme_set(theme_cowplot())

Batch.Correction <- function(P.matrix, batch.order, batch.colour, batch.shape, output.dir, name) {
  
  P.corr = as.data.frame(sva::ComBat(dat = P.matrix, batch = batch.order))
  
  t.P.corr = as.data.frame(t(P.corr))
  t.P.matrix = as.data.frame(t(P.matrix))
  
  t.P.corr$Batch = batch.colour
  t.P.matrix$Batch = batch.colour
  
  t.P.corr$Group = batch.shape
  t.P.matrix$Group = batch.shape
  
  pca.raw = prcomp(t.P.matrix[,1:(ncol(t.P.matrix)-2)], scale. = TRUE)
  pca.corr = prcomp(t.P.corr[,1:(ncol(t.P.matrix)-2)], scale. = TRUE)
  
  autoplot(pca.raw, data = t.P.matrix, colour =  "Batch", shape = "Group", label = TRUE, size = 3) + theme_cowplot()
  ggsave(paste0(output.dir,name,"_PCA_raw_Labelled.pdf"))
  
  autoplot(pca.corr, data = t.P.corr, colour = "Batch", shape = "Group", label = TRUE, size = 3) + theme_cowplot()
  ggsave(paste0(output.dir,name,"_PCA_BatchCorr_Labelled.pdf"))
  
  autoplot(pca.raw, data = t.P.matrix, colour =  "Batch", shape = "Group", size = 3) + theme_cowplot()
  ggsave(paste0(output.dir,name,"_PCA_raw_Grouped.pdf"))
  
  autoplot(pca.corr, data = t.P.corr, colour = "Batch", shape = "Group", size = 3) + theme_cowplot()
  ggsave(paste0(output.dir,name,"_PCA_BatchCorr_Grouped.pdf"))
  
  return(P.corr)
  
}
