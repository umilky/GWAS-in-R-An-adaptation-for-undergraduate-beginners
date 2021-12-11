# === Data analysis ============================================================

## installing packages 
# They are in our GitHub Repository 

install.packages(c("GenABEL.data_1.0.0.tar.gz", "GenABEL_1.8-0.tar.gz"), 
                 repos = NULL, type="source")

# association analysis of typed SNPs 
# this is specific to the data and here we use the HDL-cholestrol as the trait

# merge clinical data and principal components to create phenotype table
pheno.sub <- merge(phenotype,pcs) # data frame => [FamID CAD sex age hdl...pc10]

# rank-based inverse normal transformation of hdl
pheno.sub$phenotype <- rntransform(pheno.sub$hdl, family = "gaussian")

# remove columns that are unnecessary from the table 
pheno.sub$hdl <- NULL
pheno.sub$ldl <- NULL
pheno.sub$tg <- NULL 
pheno.sub$CAD <- NULL 

# Rename columns to match names necessary for GWAS() function 
pheno.sub <- rename(pheno.sub, replace = c(FamID = "id"))

# Include only subjects with hdl data 
pheno.sub <- pheno.sub[!is.na(pheno.sub$phenotype),]
# 1309 subjects identified to have phenotype data

print (head(pheno.sub))

# Run GWAS analysis using parallel processing 

# loading required packages 
library(plyr)

# Loading GWAS function. This function was created by Reed et al. (2015).
# Gives output as a .txt file 

source("GWAA.R")

# The GWAA function apparently takes an hour or two to run depending on data 
# Don't freak out if you don't see anything for a bit.

Start <- Sys.time()

GWAA(genodata = genotype.s, phenodata = pheno.sub, filename = gwaa.out)

End <- Sys.time()
print(End-Start)

# Reading GWAS output 

GWASoutput <- read.table("GWAS/GWASout.txt", 
                         header=TRUE, 
                         colClasses=c("character", rep("numeric",4)))

head(GWASoutput)

# on to script 6.Visualization 


