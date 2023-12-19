# Load the packages
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

DEP_Preprocess <- function(P.matrix, exp.design, output.dir, name, vsn = FALSE, 
                           keep.duplicates = TRUE) {
  
  P.output = list()
  
  # Are there any duplicated gene names?
  Flag.duplicates = P.matrix$`Gene Symbol` %>% duplicated() %>% any()
  
  if (Flag.duplicates == TRUE) {
    print("Keeping renamed duplicates")
    # Make a table of duplicated gene names
    Duplicates = P.matrix %>% group_by(`Gene Symbol`) %>% summarize(frequency = n()) %>% 
      arrange(desc(frequency)) %>% filter(frequency > 1)
    
    
    # Make unique names using the annotation in the "Gene Symbol" column as 
    #primary names and the annotation in "Accession" as name for those 
    #that do not have an gene name.
    P.unique <- make_unique(P.matrix, "Gene Symbol", "Accession", delim = ";")
    
    # Save unique dataframe
    P.output[["unique.df"]] = P.unique
  }
  
  # Values below 1 lead to negative log2 value. Therefore, make <1 TMT to NA
  P.unique[P.unique < 1] = NA
  
  # Generate a SummarizedExperiment object using an experimental design
  P.TMT <- grep("TMT_", colnames(P.unique)) # get TMT column numbers
  P.design <- readxl::read_xlsx(exp.design)
  P.se <- make_se(P.unique, P.TMT, P.design)
  
  # Plot a barplot of the protein identification overlap between samples
  plot_frequency(P.se)
  ggsave(paste0(output.dir,name,"_protein-sample_overlap.pdf"))
  
  # Less stringent filtering:
  # Filter for proteins that are identified in 2 out of 3 replicates of at least one condition
  P.filt <- filter_missval(P.se, thr = 1)

  # Turn SE to dataframe and keep duplicates with the highest value
  if (keep.duplicates == FALSE) {
    
    t.P.corr = as.data.frame(t(P.corr))
    
    P.corr$name = rownames(P.corr)
    P.corr$ID = rownames(ids)
    P.table = P.corr[,c(32,31, 1:30)]
    P.output[["P.table"]] = P.table
    
    P.corr[,1:30] = 2^P.corr[,1:30]
    P.unique = make_unique(proteins = P.corr, names = "name", ids = "ID")
    P.design <- readxl::read_xlsx(exp.design)
    P.design$label = colnames(P.corr[,1:30])
    P.se = make_se(proteins_unique = P.unique, columns = c(1:30), expdesign = P.design)
    P.output[["P.se.corr"]] = P.se
    
    print("Keeping duplicate with highest value")
    # Order protein matrix to remove duplicate with lowest absolute value
    P.matrix = P.matrix[order(P.matrix$`Gene Symbol`, -rowSums(P.matrix[,17:46], na.rm = TRUE)),]
    P.matrix = P.matrix[!duplicated(P.matrix$`Gene Symbol`), ]
    P.unique <- make_unique(P.matrix, "Gene Symbol", "Accession", delim = ";")
    
  }
  
  # Plot a barplot of the protein identification overlap between samples
  plot_frequency(P.filt)
  ggsave(paste0(output.dir,name,"_filtered_protein-sample_overlap.pdf"))
  
  # Visualize log2 transformation by boxplots for all samples per replicate
  plot_normalization(P.se, P.filt)
  ggsave(paste0(output.dir,name,"_boxplot_transformation.pdf"))
  
  # Plot intensity distributions and cumulative fraction of proteins with and without missing values
  plot_detect(P.filt)
  ggsave(paste0(output.dir,name,"_Plot_Missing-Values.pdf"))
  
  # If TRUE, do variance stabilizing normalisation (vsn) on the data
  if (vsn == TRUE) {
    P.no_vsn = P.filt
    P.filt <- normalize_vsn(P.filt)
    plot_normalization(P.no_vsn, P.filt)
    ggsave(paste0(output.dir,name,"_boxplot_after-vsn.pdf"))
    
    P.output[["SD.plot_No-VSN"]] = P.no_vsn
    P.output[["SD.plot_VSN"]] = P.filt
    
  }
  
  # The missing values seems to be concentrated to proteins with low intensities i.e.,
  # the data is missing not at random (MNAR). Thus, MNAR mode of imputting is 
  # chosen to impute data. The mode is ti impute data from gaussian distribution 
  # centered around a minimal value
  P.imp <- impute(P.filt, fun = "MinProb", q = 0.01)
  
  # Plot data before and after imputation
  plot_imputation(P.filt, P.imp)
  ggsave(paste0(output.dir,name,"_Imputation_plot.pdf"))
  
  plot_normalization(P.filt, P.imp)
  ggsave(paste0(output.dir,name,"_Imputation_boxplot_vsn-normalisation.pdf"))
  
  # Cbind log2 transformed and filtered data with accession number
  P.log.name = as.data.frame(P.imp@assays@data@listData)
  P.log.ID = P.log.name
  P.log.ID$name = rownames(P.log.name)
  P.log.ID = merge(P.unique[,(ncol(P.unique)-1):ncol(P.unique)],P.log.ID, by = "name")
  rownames(P.log.ID) = P.log.ID$name ### CHANGED TO NAME BECAUSE OF ERROR IN PHOSPHO DATA
  P.log.ID = P.log.ID[,-c(1,2)]
  P.output[["ID"]] = P.log.ID
  P.output[["name"]] = P.log.name
  
  # Return list with matrix of log2 formated values with either gene name or ID
  return(P.output)
}


