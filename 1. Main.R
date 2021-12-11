# === GWAS tutorial ============================================================
# authors: Umayeer Milky, Severin Santana

# This code is the same as the one demostrated by Reed et al. (2015). The 
# purpose of running this code is to familiarize ourselves with the different 
# steps involved in a GWAS analysis in R. Minor modifications to the original
# code has been made as required. 

# === Installing Packages ======================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("snpStats")
BiocManager::install("SNPRelate")
BiocManager::install("rtracklayer")
BiocManager::install("biomaRt")
install.packages(c("plyr", "GenABEL", "LDheatmap", "doParallel", "ggplot2"))

library("snpStats")
library("SNPRelate")
library("biomaRt")
library("plyr", "GenABEL", "LDheatmap", "doParallel", "ggplot2")
library("rtracklayer")

# === customizing data directory ===============================================

dir.path <- getwd()

# === folder management ========================================================

folder.names <- c("Raw Data", "Filtered Data", "Generated Data",
                  "GWAS", "Graphics")

for(i in 1:length(folder.names)){ 
  if(file.exists(folder.names[i]) == FALSE){
    dir.create(folder.names[i])
  } 
}

# ******************************************************************************

#path names to variables
p.data.raw <- paste(dir.path, "/", folder.names[1], "/", sep = "")
p.data.filtered <- paste(dir.path, "/", folder.names[2], "/", sep = "")
p.data.generated <- paste(dir.path, "/", folder.names[3], "/", sep = "")
p.gwas <- paste(dir.path, "/", folder.names[4], "/", sep = "")
p.graphics <- paste(dir.path, "/", folder.names[5], "/", sep = "")


# === Input Files ==============================================================
#returns character objects that are a formatted combination of input values
 
gwas.data <- lapply(c(bed="bed", bim="bim", fam="fam", gds="gds"), 
                      function(n) sprintf("%s/GWAStutorial.%s",
                                          p.data.raw, n)) 

clinical.data <- sprintf("%s/GWAStutorial_clinical.csv", p.data.raw)

onethou.fn <- lapply(c(info='info',ped='ped'), 
                     function(n) sprintf("%s/chr16_1000g_CEU.%s",
                                         p.data.raw, n))

protein.coding.coords.fname <- sprintf("%s/ProCodgene_coords.csv", p.data.raw)

# on to script 2.Reading.raw.data.R

