# === GWAS =====================================================================
# authors: Umayeer Milky, Severin Santana

# This code is an adaptation of the codes demonstrated by  Reed et al. (2015), 
# Foulkes (2016), Lima (2017) & Marees et al.(2018). 
# The code runs a genome-wide association study on genomic and phenotype data.
# The GWA results are visualized in a Manhattan plot sorted to look for 
# correlations that cross a Bonferroni corrected threshold or a less stringent 
# suggestive association threshold.
# The code is designed to run one GWA analysis at a time and 
# should not be used to analyze data from multiple populations parallel.

# GWAS main version 
# This version of the GWAS analysis runs on data downloaded from 
# Reed et al. (2015). The downloaded data is collected from humans and 
# measures genomic data along with phenotype data for coronary artery 
# disease along with lipid level of triglycerides, HDL and LDL. 

# R version
R.version.string
# "R version 4.1.2 (2021-11-01)" 

# NOTE: run the 1.Main.R before starting your session.
# === Notes ====================================================================

# • The code must be run sequentially, following the numerical order,
#   the .R script files without numerical names do not need to be run. 
#   However the code saves objects between scripts such that the environment can 
#   be closed or restarted if needed or if a crash happens. 
# • Ensure all the packages listed below are installed and loaded before 
#   running the code.
# • After running the main code please continue to the 
#   "Test" & "Simulated Data" folders to run the code with a independently 
#   downloaded data set and a simulated data set, respectively.
# • After running 3.2.Data.simulation.R, the simulated data files simulated.bed, 
#   simulated.bim, simulated.fam need to be manually moved from the 
#   "Simulated Data" folder to "a.simulated.data"

# === Script index =============================================================

# 1.Main.R        
# 2.Reading.raw.data.R
# 3.Data.filtering.R
# 4.Dataa.generation.R
# 5.Data.analysis.R
# 6.Visualization.R 

# For "Test" folder
# 2.1.Main.R        
# 2.2.Reading.raw.data.R
# 2.3.Data.filtering.R
# 2.4.Dataa.generation.R
# 2.5.Data.analysis.R
# 2.6.Visualization.R 

# For "Simulated Data" folder
# 3.1.Main.R     
# 3.2.Data.simulation.R
# 3.3.Reading.raw.data.R
# 3.4.Data.filtering.R
# 3.5.Dataa.generation.R
# 3.6.Data.analysis.R
# 3.7.Visualization.R 

# Other scripts included in this repository are required to run some functions 
# used in the aforementioned scripts. 

# === Installing Packages ======================================================
# The following packages are required to run the code 

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

# === Customizing data directory ===============================================

main.dir <- getwd()

# === Folder management ========================================================

folder.names <- c("a.raw.data", "b.working.data", "c.GWAS", "d.graphics")

for(i in 1:length(folder.names)){ 
  if(file.exists(folder.names[i]) == FALSE){
    dir.create(folder.names[i])
  } 
}

# ******************************************************************************

# Creation of paths to each folder in repository 
p.data.raw <- paste(main.dir, "/", folder.names[1], "/", sep = "")
p.work.data <- paste(main.dir, "/", folder.names[2], "/", sep = "")
p.gwas <- paste(main.dir, "/", folder.names[3], "/", sep = "")
p.graphics <- paste(main.dir, "/", folder.names[4], "/", sep = "")



# === Input Files ==============================================================
# Create character objects that are a formatted combination of input values
 
gwas.data <- lapply(c(bed="bed", bim="bim", fam="fam", gds="gds"), 
                      function(n) sprintf("%s/GWAStutorial.%s",
                                          p.data.raw, n)) 

clinical.data <- sprintf("%s/GWAStutorial_clinical.csv", p.data.raw)

onethou.fn <- lapply(c(info='info',ped='ped'), 
                     function(n) sprintf("%s/chr16_1000g_CEU.%s",
                                         p.data.raw, n))


# === Output Files =============================================================

gwaa.out <- sprintf("%s/GWASout.txt", p.gwas)
impute.out.fname <- sprintf("%s/Imputation.csv", p.gwas)


# Working data saved between each script so each script can run independently.
# Use save(data, file=working.data.fname(num)) to save data between scripts.
working.data.fname <- function(num) { sprintf("%s/working.%s.Rdata", 
                                              p.work.data, num) }
# on to script 2.Reading.raw.data.R

#___ end _______________________________________________________________________

