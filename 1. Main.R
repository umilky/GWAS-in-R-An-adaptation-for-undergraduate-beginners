# === GWAS tutorial ============================================================
# authors: Umayeer Milky, Severin Santana

# This code is an adaptation of the codes demonstrated by  Reed et al. (2015), 
# Foulkes (2016), Lima (2017) & Marees et al.(2018). The 
# purpose of this code is to familiarize ourselves with the different 
# steps involved in a GWAS analysis in R. 

# === Installing Packages ======================================================
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("snpStats")
BiocManager::install("SNPRelate")
BiocManager::install("rtracklayer")
BiocManager::install("biomaRt")
install.packages(c("plyr", "LDheatmap", "doParallel", "ggplot2"))
install.packages(c("GenABEL.data_1.0.0.tar.gz", "GenABEL_1.8-0.tar.gz"), 
                 repos = NULL, type="source")

library("snpStats")
library("SNPRelate")
library("biomaRt")
library("plyr")
library("GenABEL")
library("LDheatmap")
library("doParallel")
library("ggplot2")
library("survival")

library("rtracklayer")

# === customizing data directory ===============================================

dir.path <- getwd()

# === folder management ========================================================

folder.names <- c("a.raw.data", "b.GWAS", "c.graphics", "d.simulated.data")

for(i in 1:length(folder.names)){ 
  if(file.exists(folder.names[i]) == FALSE){
    dir.create(folder.names[i])
  } 
}

# ******************************************************************************

#path names to variables
p.data.raw <- paste(dir.path, "/", folder.names[1], "/", sep = "")
p.gwas <- paste(dir.path, "/", folder.names[2], "/", sep = "")
p.graphics <- paste(dir.path, "/", folder.names[3], "/", sep = "")
p.simulated.data <- paste(dir.path, "/", folder.names[4], "/", sep = "")


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

simulated.bed <- sprintf("%s/simulated.bed", p.data.raw)


# === Output Files =============================================================

gwaa.out <- sprintf("%s/GWASout.txt", p.gwas)
impute.out.fname <- sprintf("%s/Imputation.csv", p.gwas)

# on to script 2.Reading.raw.data.R

