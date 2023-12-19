library(matrixStats)
library(ggplot2)
library(cowplot)
library(readxl)
library(DEP)
library(dplyr)
library(ggfortify)
library(reshape2)
library(sva)
library(data.table)
library(enrichR)

# Set default theme to cowplot
theme_set(theme_cowplot())

DEP_analysis <- function(P.corr, ids, output.dir, name, exp.design, p.cut = 0.01, q.cut = 0.05, log2fc = log2(1.5)) {
  
  P.output = list()
  
  P.corr$name = rownames(P.corr)
  P.corr$ID = rownames(ids)
  values.col = (ncol(P.corr)-2)
  P.table = P.corr[,c(ncol(P.corr),(ncol(P.corr)-1), 1:values.col)]
  P.output[["P.table"]] = P.table
  
  P.corr[,1:values.col] = 2^P.corr[,1:values.col]
  P.unique = make_unique(proteins = P.corr, names = "name", ids = "ID")
  P.design <- readxl::read_xlsx(exp.design)
  P.design$label = colnames(P.corr[,1:values.col])
  P.se = make_se(proteins_unique = P.unique, columns = c(1:values.col), expdesign = P.design)
  P.output[["P.se.corr"]] = P.se
  
  # PCOS w5 vs. PCOS w0, adjusted p-value
  P.diff <- test_diff(P.se, type = "manual", test = "PCOS_w5_vs_PCOS_w0")
  P.diff@elementMetadata$PCOS_w5_vs_PCOS_w0_p.adj = p.adjust(p = P.diff@elementMetadata$PCOS_w5_vs_PCOS_w0_p.val,
                                                             method = "BH")
  dep <- add_rejections(P.diff, alpha = q.cut, lfc = log2fc)
  P.res = get_results(dep)
  P.output[["PCOS_w5_w0_p.adj"]] = P.res
  plot_volcano(dep, contrast = "PCOS_w5_vs_PCOS_w0", label_size = 2, add_names = TRUE)
  ggsave2(paste0(output.dir,name,"_Volcano_PCOS_w5-w0_p.adj.pdf"))
  plot = plot_pca(dep, x = 1, y = 2, n = 100, point_size = 4, indicate = c("replicate", "condition"))
  ggsave2(paste0(output.dir,name,"_PCA_Labelled_PCOS_w5-w0.pdf"))
  plot_pca(dep, x = 1, y = 2, n = 100, point_size = 4, indicate = "condition")
  ggsave2(paste0(output.dir,name,"_PCA_PCOS_w5-w0.pdf"))
  plot_cor(dep, significant = FALSE, lower = 0, upper = 1, pal = "Reds")
  ggsave2(paste0(output.dir,name,"_PCOS_w5-w0_Pearson_correlation_matrix.pdf"))
  
  # PCOS w5 vs. PCOS w0, normal p-value
  P.diff <- test_diff(P.se, type = "manual", test = "PCOS_w5_vs_PCOS_w0")
  P.diff@elementMetadata$PCOS_w5_vs_PCOS_w0_p.adj = P.diff@elementMetadata$PCOS_w5_vs_PCOS_w0_p.val
  
  dep <- add_rejections(P.diff, alpha = p.cut, lfc = log2fc)
  P.res = get_results(dep)
  P.output[["PCOS_w5_w0_p.value"]] = P.res
  plot_volcano(dep, contrast = "PCOS_w5_vs_PCOS_w0", label_size = 2, add_names = TRUE)
  ggsave2(paste0(output.dir,name,"_Volcano_PCOS_w5-w0_p.value.pdf"))
  
  # PCOS w0 vs. Control, adjusted p-value
  P.diff <- test_diff(P.se, type = "manual", test = "PCOS_w0_vs_Ctrl")
  P.diff@elementMetadata$PCOS_w0_vs_Ctrl_p.adj = p.adjust(p = P.diff@elementMetadata$PCOS_w0_vs_Ctrl_p.val,
                                                             method = "BH")
  dep <- add_rejections(P.diff, alpha = q.cut, lfc = log2fc)
  P.res = get_results(dep)
  P.output[["PCOS_w0-Ctrl_p.adj"]] = P.res
  
  plot_volcano(dep, contrast = "PCOS_w0_vs_Ctrl", label_size = 2, add_names = TRUE)
  ggsave2(paste0(output.dir,name,"_Volcano_PCOS_w0-Ctrl_p.adj.pdf"))
  plot_pca(dep, x = 1, y = 2, n = 100, point_size = 4, indicate = c("replicate", "condition"))
  ggsave2(paste0(output.dir,name,"_PCA_Labelled_PCOS_w0-Ctrl.pdf"))
  plot_pca(dep, x = 1, y = 2, n = 100, point_size = 4, indicate = "condition")
  ggsave2(paste0(output.dir,name,"_PCA_PCOS_w0-Ctrl.pdf"))
  plot_cor(dep, significant = FALSE, lower = 0, upper = 1, pal = "Reds")
  ggsave2(paste0(output.dir,name,"_PCOS_w0-Ctrl_Pearson_correlation_matrix.pdf"))
  
  # PCOS w0 vs. Control, normal p-value
  P.diff <- test_diff(P.se, type = "manual", test = "PCOS_w0_vs_Ctrl")
  P.diff@elementMetadata$PCOS_w0_vs_Ctrl_p.adj = P.diff@elementMetadata$PCOS_w0_vs_Ctrl_p.val
  
  dep <- add_rejections(P.diff, alpha = p.cut, lfc = log2fc)
  P.res = get_results(dep)
  P.output[["PCOS_w0-Ctrl_p.value"]] = P.res
  plot_volcano(dep, contrast = "PCOS_w0_vs_Ctrl", label_size = 2, add_names = TRUE)
  ggsave2(paste0(output.dir,name,"_Volcano_PCOS_w0-Ctrl_p.value.pdf"))
  
  return(P.output)
}

Generate_DEP_table <- function(x, name, output.dir, alpha = 0.05, log2fc = log2fc_set) {
  
  x.out = x[,c(1:4, 7)]
  colnames(x.out) = c("name", "ID", "p.value", "p.adj", "Log2FC")
  
  x.out$significant = (x.out$p.value < alpha & x.out$Log2FC < -log2fc) | 
    (x.out$p.value < alpha & x.out$Log2FC > log2fc)
  
  writexl::write_xlsx(x = x.out, path = paste0(output.dir, name, "_DEP_table.xlsx"))
  
  return(x.out)
  
}
