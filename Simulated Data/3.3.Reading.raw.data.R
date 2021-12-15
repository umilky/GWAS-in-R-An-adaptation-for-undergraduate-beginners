# === Reading Data in R ========================================================

# *Note 3.1.Main.R is required to run first to access folder paths and to input 
# raw data required to run this script.
# All scripts must be run in numerical order 

# This scrips creates and formats objects out of raw data for 
# analysis in future scripts. 
# Objects involving genomic data and phenotype data are created in this script.

# 1) Genomic data creation
# 2) Phenotype data creation
# 3) Data manipulation 


# === 1) Genomic data creation =================================================

# Reading of PLINK files to create a list of genomic data.
# PLINK is a genome association analysis toolkit.
# .bed files contain individual genotypes.
# .bin files contain the locations of all SNPs in the data in the genome.
# .fam files contain information about the individual including family pedigrees.

# Creating genomic data matrix from PLINK
geno.me <- read.plink(gwas.data$bed, gwas.data$bim, gwas.data$fam)

# Creating the genotype matrix from the geno.me matrix
genotype <- geno.me$genotype 
head(genotype)

# Creating a data frame for the SNP information from geno.me matrix
geno.bim <- geno.me$map
# renaming of column names for later data frame merge
colnames(geno.bim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
print(head(geno.bim))


# Creating SNP dataframe from the geno.me list 
SNP.genome <- geno.me$map
# Renaming of column names for later data frame merge
colnames(SNP.genome) <- c("chromosome", "SNP", "gen.dist", "position", "A1", 
                          "A2")
print(head(SNP.genome))

# === 2) Phenotype data creation================================================

# Reading phenotype data from the clinical data csv file. 
# Using colClasses saves data importing time
# colClasses creates a vector of column classes used for tabular reading. 
# This is required for subsetting genotype in the next step. 

# Creating data frame from phenotype csv file
phenotype <- read.csv(clinical.data,
                     colClasses = c("character", rep("numeric", 4)))
rownames(phenotype) <- phenotype$id
head(phenotype)


# === 3) Data manipulation =====================================================

# Subsetting genotype object to only include individuals with known phenotype 
genotype.s <- genotype[phenotype$id, ]


# ******************************************************************************

# Write genotype.s, geno.bim & phenotype for future use.
# *Note this function can take a few minutes
save(genotype.s, geno.bim, phenotype, file = working.data.fname(3.3))

# on to script 3.4.Data.filtering.R 

#___ end _______________________________________________________________________
