# Author: Gustaw Eriksson
# Date: 2021-09-16
# Description: Aalysis of proteomics data of fat and muscle received from Anna Benrick

# Load packages
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
library(org.Hs.eg.db)

# Set default theme to cowplot
theme_set(theme_cowplot())

# Set work directory to directory where script is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load function scripts in /bin
source(file = "bin/DEP-Preprocessing.R")
source(file = "bin/DEP-ComBat.Batch.Correction.R")
source(file = "bin/DEP_Differential.Expression.Analysis.R")
source(file = "bin/DEP.Enrichment_Analysis.R")

# Setting parameters

# Is VSN data normalisation to be done (Currently not used)
VSN_norm = FALSE

# Set alpha threshold (p-value significans cut-off)
p_set = 0.05
q_set = 0.1

#Set log2 fold change threshold 
log2fc_set = 0.5 # 50% difference

# Load either fat or muscle data and reformatting the data
Load_muscle = FALSE
Load_muscle_phospho = FALSE
Load_fat = FALSE
Load_fat_phospho = TRUE
Load_muscle_protein_outlier_removed = FALSE
Load_muscle_phospho_outlier_removed = FALSE

# Set output directory in project fold
data_output = ""

if (dir.exists(data_output) == FALSE) {
  print(paste("Creating directory", data_output))
  dir.create(path = data_output)
} else if (dir.exists(data_output) == TRUE) {
  print(paste(data_output, "directory exists"))
}

if (Load_muscle == TRUE) {
  print("Loading muscle protein")
  data_raw = readxl::read_xlsx("Data/MUSCLE_TMT_Intensity_edited.xlsx")
  data_raw = subset(data_raw, select = -c(3,18:47, 63, 79))
  data_raw = cbind(data_raw[,1:16], data_raw[,17:46][order(sub('.*m ', '', colnames(data_raw[,17:46])))])
  colnames(data_raw)[17:46] = c(paste0("TMT_m", seq(30)))
  
  # Setting muscle parameters
  exp.design.matrix = "Data/Muscle_experimental_design.xlsx"
  tissue.type = "Muscle"
  
} else if (Load_muscle_protein_outlier_removed == TRUE) {
  print("Loading muscle protein")
  data_raw = readxl::read_xlsx("Data/MUSCLE_TMT_Intensity_edited.xlsx")
  data_raw = subset(data_raw, select = -c(3,18:47, 63, 79))
  data_raw = cbind(data_raw[,1:16], data_raw[,17:46][order(sub('.*m ', '', colnames(data_raw[,17:46])))])
  colnames(data_raw)[17:46] = c(paste0("TMT_m", seq(30)))
  
  # Removing outlier
  data_raw = subset(data_raw, select = -c(34))
  
  # Setting muscle parameters
  exp.design.matrix = "Data/Muscle_experimental_design_NO_OUTLIER.xlsx"
  tissue.type = "Muscle"
  
} else if (Load_muscle_phospho == TRUE) {
  print("Loading muscle phosphorilation")
  data_raw = readxl::read_xlsx("Data/MUSCLE_TMT_Phospho_edited.xlsx")
  
  
  data_raw = subset(data_raw, select = -c(11:42, 58, 74:75))
  data_raw = cbind(data_raw[,1:10], data_raw[,11:40][order(sub('.*m ', '', colnames(data_raw[,11:40])))])
  colnames(data_raw)[11:40] = c(paste0("TMT_m", seq(30)))
  
  # Remove rows of which one phosposites maps to several protein
  data_raw = data_raw[!grepl(";", data_raw$`Master Protein Accessions`),]
  
  # Comvert Uniprot accession numbers to gene symbols
  gene_symbols = select(org.Hs.eg.db, data_raw$`Master Protein Accessions`, "SYMBOL", "UNIPROT")
  gene_symbols = gene_symbols[!duplicated(gene_symbols$UNIPROT),]

  #data_raw$`Gene Symbol` = data_raw$`Master Protein Accessions`
  data_raw = merge(data_raw, gene_symbols, 
                      by.x = "Master Protein Accessions", 
                      by.y = "UNIPROT", 
                      all.x = TRUE)
  
  # Put the accession numbers first
  colnames(data_raw)[1] = "Accession"
  colnames(data_raw)[41] = "Gene Symbol"
  data_raw = data_raw[,c(1, 41, 2:40)]
  data_raw$`Gene Symbol`[is.na(data_raw$`Gene Symbol`)] = data_raw$Accession[is.na(data_raw$`Gene Symbol`)]
  
  # Setting muscle phospho parameters
  exp.design.matrix = "Data/Muscle_experimental_design.xlsx"
  tissue.type = "Muscle_Phospo"
  
} else if (Load_muscle_phospho_outlier_removed == TRUE) {
  print("Loading muscle phosphorilation and removing outlier")
  data_raw = readxl::read_xlsx("Data/MUSCLE_TMT_Phospho_edited.xlsx")
  
  
  data_raw = subset(data_raw, select = -c(11:42, 58, 74:75))
  data_raw = cbind(data_raw[,1:10], data_raw[,11:40][order(sub('.*m ', '', colnames(data_raw[,11:40])))])
  colnames(data_raw)[11:40] = c(paste0("TMT_m", seq(30)))
  
  # Removing outlier
  data_raw = subset(data_raw, select = -c(28))
  
  # Remove rows of which one phosposites maps to several protein
  data_raw = data_raw[!grepl(";", data_raw$`Master Protein Accessions`),]
  
  # Comvert Uniprot accession numbers to gene symbols
  gene_symbols = select(org.Hs.eg.db, data_raw$`Master Protein Accessions`, "SYMBOL", "UNIPROT")
  gene_symbols = gene_symbols[!duplicated(gene_symbols$UNIPROT),]
  
  #data_raw$`Gene Symbol` = data_raw$`Master Protein Accessions`
  data_raw = merge(data_raw, gene_symbols, 
                   by.x = "Master Protein Accessions", 
                   by.y = "UNIPROT", 
                   all.x = TRUE)
  
  # Put the accession numbers first
  colnames(data_raw)[1] = "Accession"
  colnames(data_raw)[40] = "Gene Symbol"
  data_raw = data_raw[,c(1, 40, 2:39)]
  data_raw$`Gene Symbol`[is.na(data_raw$`Gene Symbol`)] = data_raw$Accession[is.na(data_raw$`Gene Symbol`)]
  
  # Setting muscle phospho parameters
  exp.design.matrix = "Data/Muscle_experimental_design_NO_OUTLIER.xlsx"
  tissue.type = "Muscle_Phospo"
  
}else if (Load_fat == TRUE) {
  print("Loading fat protein")
  data_raw = readxl::read_xlsx("Data/FAT_TMT_Intensity_edited.xlsx")
  
  
  data_raw = subset(data_raw, select = -c(3,18:47, 63, 79))
  data_raw = cbind(data_raw[,1:16], data_raw[,17:46][order(sub('.*f ', '', colnames(data_raw[,17:46])))])
  colnames(data_raw)[17:46] = c(paste0("TMT_f", seq(30)))
  
  # Setting fat parameters
  exp.design.matrix = "Data/Fat_experimental_design.xlsx"
  tissue.type = "Fat"

} else if (Load_fat_phospho == TRUE) {
  print("Loading fat phosphorilation")
  data_raw = readxl::read_xlsx("Data/FAT_TMT_Phospho_edited.xlsx")
  
  data_raw = subset(data_raw, select = -c(1, 12:43, 59, 75:76))
  data_raw = cbind(data_raw[,1:10], data_raw[,11:40][order(sub('.*f ', '', colnames(data_raw[,11:40])))])
  colnames(data_raw)[11:40] = c(paste0("TMT_f", seq(30)))
  
  # Remove rows of which one phosposites maps to several protein
  ## sum(grepl(";", data_raw$`Master Protein Accessions`))
  data_raw = data_raw[!grepl(";", data_raw$`Master Protein Accessions`),]
  
  # Comvert Uniprot accession numbers to gene symbols
  gene_symbols = select(org.Hs.eg.db, data_raw$`Master Protein Accessions`, "SYMBOL", "UNIPROT")
  gene_symbols = gene_symbols[!duplicated(gene_symbols$UNIPROT),]
  
  #data_raw$`Gene Symbol` = data_raw$`Master Protein Accessions`
  data_raw = merge(data_raw, gene_symbols, 
                   by.x = "Master Protein Accessions", 
                   by.y = "UNIPROT", 
                   all.x = TRUE)
  
  # Put the accession numbers first
  colnames(data_raw)[1] = "Accession"
  colnames(data_raw)[41] = "Gene Symbol"
  data_raw = data_raw[,c(1, 41, 2:40)]
  data_raw$`Gene Symbol`[is.na(data_raw$`Gene Symbol`)] = data_raw$Accession[is.na(data_raw$`Gene Symbol`)]
  

  # Setting fat phospho parameters
  exp.design.matrix = "Data/Fat_experimental_design.xlsx"
  tissue.type = "Fat"
}

# Preprocessing fo the data to rename duplicates, filter and impute missing values, and plot the results
data.preprocess = DEP_Preprocess(P.matrix = data_raw, 
                                   exp.design = exp.design.matrix, 
                                   output.dir = data_output, 
                                   name = tissue.type,
                                   vsn = VSN_norm, 
                                   keep.duplicates = TRUE)

# Batch correct the data using ComBat batch correction
if (Load_muscle_phospho_outlier_removed == TRUE | Load_muscle_protein_outlier_removed == TRUE) {
  print("Processing after removing outlier")
  data.corr = Batch.Correction(P.matrix = data.preprocess$name, 
                               batch.order = c(rep(1,5), rep(2,5), rep(1,5), rep(2,4), rep(1,5), rep(2,5)), 
                               batch.colour = c(rep("Batch 1", 5), rep("Batch 2", 5), rep("Batch 1", 5),
                                                rep("Batch 2", 4), rep("Batch 1", 5), rep("Batch 2", 5)), 
                               batch.shape = c(rep("Control", 10), rep("PCOS w0", 9), rep("PCOS w5", 10)), 
                               name = tissue.type,
                               output.dir = data_output)
  
} else if (Load_muscle_phospho_outlier_removed == FALSE) {
  data.corr = Batch.Correction(P.matrix = data.preprocess$name, 
                               batch.order = c(rep(1,5), rep(2,5), rep(1,5), rep(2,5), rep(1,5), rep(2,5)), 
                               batch.colour = c(rep("Batch 1", 5), rep("Batch 2", 5), rep("Batch 1", 5),
                                                rep("Batch 2", 5), rep("Batchha 1", 5), rep("Batch 2", 5)), 
                               batch.shape = c(rep("Control", 10), rep("PCOS w0", 10), rep("PCOS w5", 10)), 
                               name = tissue.type,
                               output.dir = data_output)
}

# Do differential expression analysis of protein data
data.DEP = DEP_analysis(P.corr = data.corr,
                       ids = data.preprocess$ID,
                       name = tissue.type,
                       output.dir = data_output,
                       exp.design = exp.design.matrix,
                       p.cut = p_set,
                       q.cut = q_set,
                       log2fc = log2fc_set)

data.w5_w0 = Generate_DEP_table(x = data.DEP$PCOS_w5_w0_p.adj, name = paste0(tissue.type, "_PCOS_w0_w5"),
                               output.dir = data_output, alpha = p_set)
data.Ctrl_w0 = Generate_DEP_table(x = data.DEP$`PCOS_w0-Ctrl_p.adj`, name = paste0(tissue.type, "_Ctrl_PCOS-w0"),
                                 output.dir = data_output, alpha = p_set)

# Enrichment analysis on significant samples based on p-value
Enrichment.Analysis.w5_w0 = Enrichment.Analysis(DEP.table = data.w5_w0, 
                                   name = paste0(tissue.type, "_PCOS_w0_w5"),
                                   output.dir = data_output, 
                                   selected.db = c("GO_Biological_Process_2021", "GO_Cellular_Component_2021",
                                                   "GO_Molecular_Function_2021"))

Enrichment.Analysis.Ctrl_w0 = Enrichment.Analysis(DEP.table = data.Ctrl_w0, 
                                     name = paste0(tissue.type, "_Ctrl_PCOS-w0"),
                                     output.dir = data_output, 
                                     selected.db = c("GO_Biological_Process_2021", "GO_Cellular_Component_2021",
                                                     "GO_Molecular_Function_2021"))

# Enrichment analysis on significant samples based on q-value
data.w5_w0.Qsign = data.w5_w0
data.w5_w0.Qsign$significant = "FALSE"
data.w5_w0.Qsign$significant[which(data.w5_w0.Qsign$p.adj < q_set)] = "TRUE"

data.Ctrl_w0.Qsign = data.Ctrl_w0
data.Ctrl_w0.Qsign$significant = "FALSE"
data.Ctrl_w0.Qsign$significant[which(data.Ctrl_w0.Qsign$p.adj < q_set)] = "TRUE"

Enrichment.Analysis.w5_w0 = Enrichment.Analysis(DEP.table = data.w5_w0.Qsign, 
                                                name = paste0(tissue.type, "_PCOS_w0_w5_padj"),
                                                output.dir = data_output, 
                                                selected.db = c("GO_Biological_Process_2021", "GO_Cellular_Component_2021",
                                                                "GO_Molecular_Function_2021"))

Enrichment.Analysis.Ctrl_w0 = Enrichment.Analysis(DEP.table = data.Ctrl_w0.Qsign, 
                                                  name = paste0(tissue.type, "_Ctrl_PCOS-w0_padj"),
                                                  output.dir = data_output, 
                                                  selected.db = c("GO_Biological_Process_2021", "GO_Cellular_Component_2021",
                                                                  "GO_Molecular_Function_2021"))

# Output and save batch corrected table
data.corr_df = setDT(data.corr, keep.rownames = TRUE)
colnames(data.corr_df)[1] = "Protein" 
writexl::write_xlsx(x = data.corr_df, 
                    path = paste0(data_output, tissue.type,"_ComBat_Protein_table.xlsx"))

if (Load_muscle_phospho | Load_fat_phospho == TRUE) {
  
  # Extract phospho data from initial raw data frame
  DEP_phospho.Ctrl_w0 = merge(data.Ctrl_w0, data.preprocess$unique.df, by = "name", all.x = TRUE)
  DEP_phospho.Ctrl_w0 = merge(DEP_phospho.Ctrl_w0, data.corr_df, by.x = "name", by.y = "Protein", all.x = TRUE)
  DEP_phospho.Ctrl_w0 = subset(DEP_phospho.Ctrl_w0, select = c(8, 1, 7, 3:6, 9:17, 49:78))
  colnames(DEP_phospho.Ctrl_w0)[2] = "Gene Symbol unique"
  writexl::write_xlsx(x = DEP_phospho.Ctrl_w0, 
                      path = paste0(data_output, tissue.type,"_phospho_ComBat_merged_Ctrl_W0.xlsx"))

  DEP_phospho.w5_w0 = merge(data.w5_w0, data.preprocess$unique.df, by = "name", all.x = TRUE)
  DEP_phospho.w5_w0 = merge(DEP_phospho.w5_w0, data.corr_df, by.x = "name", by.y = "Protein", all.x = TRUE)
  DEP_phospho.w5_w0 = subset(DEP_phospho.w5_w0, select = c(8, 1, 7, 3:6, 9:17, 49:78))
  colnames(DEP_phospho.w5_w0)[2] = "Gene Symbol unique"
  writexl::write_xlsx(x = DEP_phospho.w5_w0, 
                      path = paste0(data_output, tissue.type,"_phospho_ComBat_merged_w5_W0.xlsx"))

} else if (Load_muscle | Load_fat == TRUE) {
  # Extract protein data from initial raw data frame
  DEP_protein.Ctrl_w0 = merge(data.Ctrl_w0, data.preprocess$unique.df, by = "name", all.x = TRUE)
  DEP_protein.Ctrl_w0 = merge(DEP_protein.Ctrl_w0, data.corr_df, by.x = "name", by.y = "Protein", all.x = TRUE)
  DEP_protein.Ctrl_w0 = subset(DEP_protein.Ctrl_w0, select = c(1, 7, 3:6, 8, 54:83))
  colnames(DEP_protein.Ctrl_w0)[1] = "Gene Symbol"
  writexl::write_xlsx(x = DEP_protein.Ctrl_w0, 
                      path = paste0(data_output, tissue.type,"_protein_ComBat_merged_Ctrl_W0.xlsx"))

  DEP_protein.w5_w0 = merge(data.w5_w0, data.preprocess$unique.df, by = "name", all.x = TRUE)
  DEP_protein.w5_w0 = merge(DEP_protein.w5_w0, data.corr_df, by.x = "name", by.y = "Protein", all.x = TRUE)
  DEP_protein.w5_w0 = subset(DEP_protein.w5_w0, select = c(1, 7, 3:6, 8, 54:83))
  colnames(DEP_protein.w5_w0)[2] = "Gene Symbol"
  writexl::write_xlsx(x = DEP_protein.w5_w0, 
                      path = paste0(data_output, tissue.type,"_protein_ComBat_merged_w5_W0.xlsx"))
} 
