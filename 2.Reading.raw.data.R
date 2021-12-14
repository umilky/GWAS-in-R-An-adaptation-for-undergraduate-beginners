# === Reading Data in R ========================================================

# Script 1 is required to run first to access folder paths and to input raw 
# data required to run this script.


# Read PLINK files to create a list
# PLINK is a genome association analysis toolkit.
# bed files contain individual genotypes.
# bin files contain the locations of all SNPs in the data in the genome.
# fam files contain information about the individual including family pedigrees.

geno.me <- read.plink(gwas.data$bed, gwas.data$bim, gwas.data$fam)

# creating the genotype object from the geno.me matrix

genotype <- geno.me$genotype 
print(genotype)
head(genotype)

# obtain the SNP information from geno.me matrix
geno.bim <- geno.me$map
colnames(geno.bim) <- c("chr", "SNP", "gen.dist", "position", "A1", "A2")
print(head(geno.bim))


#creating SNP object from the geno.me list 
SNP.genome <- geno.me$map
colnames(SNP.genome) <- c("chromosome", "SNP", "gen.dist", "position", "A1", 
                          "A2")

print(head(SNP.genome))


# reading phenotype data from the clinical file 
# Using colClasses saves data importing time
# colClasses creates a vector of column classes used for tabular reading. 
# This is required for subsetting genotype in the next step. 

phenotype <- read.csv(clinical.data,
                     colClasses = c("character", "factor", "factor", 
                                    rep("numeric", 4)))
rownames(phenotype) <- phenotype$FamID
head(phenotype)

#subsetting genotype to only include individuals with known phenotype 

genotype.s <- genotype[phenotype$FamID, ]


### we need to save these objects somehow so they can run in future scripts

# the authors use this
# Write genotype, genoBim, clinical for future use
# does not work ####
#save(genotype, geno.bim, phenotype, file = working.data.fname(1))

# on to script 3.Data.filtering 

#___ end _______________________________________________________________________