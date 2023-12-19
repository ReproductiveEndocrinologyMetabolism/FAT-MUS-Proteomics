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
library(enrichplot)

# Set default theme to cowplot
theme_set(theme_cowplot())

Enrichment.Analysis <- function(DEP.table, selected.db, output.dir, name) {
  
  setEnrichrSite("Enrichr")
  DEP.signif = DEP.table$name[DEP.table$significant == TRUE]
  DEP.signif.pos = DEP.table$name[DEP.table$significant == TRUE & DEP.table$Log2FC >= 0.5]
  DEP.signif.neg = DEP.table$name[DEP.table$significant == TRUE & DEP.table$Log2FC <= -0.5]
  enriched.all = enrichr(DEP.signif, selected.db)
  enriched.pos = enrichr(DEP.signif.pos, selected.db)
  enriched.neg = enrichr(DEP.signif.neg, selected.db)
  
  reg = "all"
  for (i in names(enriched.all)) {
    print(i)
    Enrich.plot = plotEnrich(enriched.all[[i]], showTerms = 20, numChar = 200, y = "Count", orderBy = "P.value",
                             title = i)
    ggsave(paste0(output.dir, name, "_",i, "_", reg, "_enrichR.pdf"), plot = Enrich.plot, 
           width = 10)
    writexl::write_xlsx(x = enriched.all[[i]], path = paste0(output.dir, name, "_",i, "_", reg, "_enrichR_table.xlsx"))
    
  }
  
  reg = "pos"
  for (i in names(enriched.pos)) {
    print(i)
    Enrich.plot = plotEnrich(enriched.pos[[i]], showTerms = 20, numChar = 200, y = "Count", orderBy = "P.value",
                             title = i)
    ggsave(paste0(output.dir, name, "_",i, "_", reg, "_enrichR.pdf"), plot = Enrich.plot, 
           width = 10)
    writexl::write_xlsx(x = enriched.pos[[i]], path = paste0(output.dir, name, "_",i, "_", reg, "_enrichR_table.xlsx"))
    
  }
  
  reg = "neg"
  for (i in names(enriched.neg)) {
    print(i)
    Enrich.plot = plotEnrich(enriched.neg[[i]], showTerms = 20, numChar = 200, y = "Count", orderBy = "P.value",
                             title = i)
    ggsave(paste0(output.dir, name, "_",i, "_", reg, "_enrichR.pdf"), plot = Enrich.plot, 
           width = 10)
    writexl::write_xlsx(x = enriched.neg[[i]], path = paste0(output.dir, name, "_",i, "_", reg, "_enrichR_table.xlsx"))
    
  }
  
}