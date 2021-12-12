# === Data analysis ============================================================

## installing packages 
# They are in our GitHub Repository 

### we could create a folder called "packages" and add them in there? 

# association analysis of typed SNPs 
# this is specific to the data and here we use the HDL-cholestrol as the trait

# merge clinical data and principal components to create phenotype table

pheno.sub <- merge.data.frame(phenotype,pcs) # data frame => [FamID CAD sex age hdl...pc10]

# rank-based inverse normal transformation of hdl
pheno.sub$phenotype <- rntransform(pheno.sub$hdl, family = "gaussian")

# show that the assumptions of normality met after transformation
par(mfrow=c(1,2))
hist(pheno.sub$hdl, main = "Histogram of HDL", xlab = "HDL")
hist(pheno.sub$phenotype, main = "Histogram of Transformed HDL",
     xlab = "Transformed HDL")

# remove columns that are unnecessary from the table 
pheno.sub$hdl <- NULL
pheno.sub$ldl <- NULL
pheno.sub$tg <- NULL 
pheno.sub$CAD <- NULL 


# Rename columns to match names necessary for GWAS() function 
names(pheno.sub)[names(pheno.sub ) == "FamID"] <- "id"

# Include only subjects with hdl data 
pheno.sub <- pheno.sub[!is.na(pheno.sub$phenotype),]
# 1309 subjects identified to have phenotype data

print (head(pheno.sub))


# Run GWAS analysis using parallel processing 

# loading required packages 
#### we have this in 1.main 
#library(plyr)

# Loading GWAS function. This function was created by Reed et al. (2015).
# Gives output as a .txt file 

source("GWAA.R")

# The GWAA function takes an hour or two to run depending on data and your computer
# Don't get discouraged if you don't see anything for a bit.

Start <- Sys.time()

GWAA(genodata = genotype.s, phenodata = pheno.sub, filename = gwaa.out)

End <- Sys.time()
print(End-Start)

# Reading GWAS output 
GWASoutput <- read.table(gwaa.out, sep = "",header=TRUE, colClasses=c("character"))
head(GWASoutput)

# Calculating the -log10 of the p-values 
GWASoutput$p.value <- as.numeric(GWASoutput$p.value)
GWASoutput$neg.logp <- log10(GWASoutput$p.value)
head(GWASoutput)

# Merge output with geno.bim by SNP name to add position and chromosome number
# this does not work and kills the whole thing :(
GWASoutput <- merge(GWASoutput, geno.bim[,c("SNP", "chr", "position")], by = "SNP")
head(GWASoutput)

# Order SNPs by significance 
GWASoutput <- arrange(GWASoutput, neg.logp)
print(head(GWASoutput))

# on to script 6.Visualization 


