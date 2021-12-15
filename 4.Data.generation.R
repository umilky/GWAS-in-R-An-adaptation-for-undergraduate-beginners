# === Data Generation ==========================================================

# *Note 1.Main. R is required to run first to access folder paths and to input 
# raw data required to run this script.
# All scripts must be run in numerical order 
# Load data created in previous script, 3.Data.filtering.R
load(working.data.fname(3))

# The following script generates principle components needed for genome-wide
# association analyses along with imputed SNPs which are usefull for association
# identification. The imputed SNPs will be visualized on the Manhattan plot. 

# 1) Principal Component Generation
# 2) Imputation of non-typed SNPs


# === 1) Principal Component Generation ========================================

# Principal Components (PCs) based on observed genotype data capture 
# information on substructure - genetic diversity in an apparently homogenous 
# population that is caused by population genetic history (migration, 
# selection, and/or ethnic integration). These substructures cannot be inferred 
# from self-reported race and ethnicity variables.

# Generating principal components for modeling 
# set LD threshold to 0.2
ld.cutoff <- 0.2
# Set seed for reproducibility 
set.seed(1000)
# Create vector of sample ids from genomic data
genosample.ids <- rownames(genotype.s)
# Prune SNPs
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.cutoff,
                          sample.id = genosample.ids, 
                          snp.id = colnames(genotype.s)) 
snpset.pca <- unlist(snpSUB, use.names = FALSE)
cat(length(snpset.pca), "\n") #72665 SNPs will be used in PCA analysis

# Generate principle components from SNPs within LD threshold 
pca <- snpgdsPCA(genofile, sample.id = genosample.ids,
                 snp.id = snpset.pca, num.thread = 1)

# Find and record first 10 principal components
# pcs will be a N:10 matrix. Each column is a principal component. 
pcs <- data.frame(FAMID = pca$sample.id, pca$eigenvect[,1:10],
                  stringsAsFactors = FALSE)
colnames(pcs)[2:11] <- paste("pc", 1:10, sep ="")
# Rename "FAMID" to "FamID" to match other files
names(pcs)[names(pcs) == "FAMID"] <- "FamID"
print(head(pcs))


# === 2) Imputation of non-typed SNPs ==========================================

# Genotypes of untyped SNPs may have a functional relationship to the 
# outcome and therefore provide additional power for identifying association.
# Imputation of SNPs involve estimating missing genotypes from the haplotype 
# or genotype reference panel. Here, we use the 1000 Genomes data for this. 

# Read in 1000g data for given chromosome 16 
thougeno <- read.pedfile(onethou.fn$ped, snps = onethou.fn$info, which = 1)

# Obtain genotype data for given chromosome
geno.matrix <- thougeno$genotypes

# Obtain the chromosome position for each SNP
support <- thougeno$map
colnames(support) <- c("SNP", "position", "N1", "N2")
head(support)

# Imputation of non-typed 1000g SNPs
# Subset for SNPs on given chromosome
pres.SNPS <- colnames(genotype.s)
pres.dat.chr <- geno.bim[geno.bim$SNP %in% pres.SNPS & geno.bim$chr==16,]
target.SNPS <- pres.dat.chr$SNP

# Subset 1000g data for our SNPs
# "missing" and "present" are snpMatrix objects needed for imputation rules
is.present <- colnames(geno.matrix) %in% target.SNPS
missing <- geno.matrix[,!is.present]
print(missing)

present <- geno.matrix[, is.present]
print(present)

# Obtain positions of SNPs to be used for imputation rules
pos.present <- support$position[is.present]
pos.missing <- support$position[!is.present]

# Calculate and store the imputation rules using snp.imputation()
rules <- snp.imputation(present, missing, pos.present, pos.missing)

# Remove failed imputations 
rules <- rules[can.impute(rules)]
cat("Imputation rules for", length(rules), "SNPs were estimated\n")

# Quality control for imputation certainty and MAF
# Set thresholds
r2threshold <- 0.7
minor <- 0.01

# Filter on imputation certainty and MAF
rules <- rules[imputation.r2(rules) >= r2threshold]
cat(length(rules), 
    "imputation rules remain after uncertain imputations were removed\n")
# 162565 imputation rules remain after uncertain imputation were removed

rules <- rules[imputation.maf(rules) >= minor]
cat(length(rules), "imputation rules remain after MAF filtering\n")
# 162565 imputation rules remain after MAF filtering

# Obtain posterior expectation of genotypes of imputed SNPs
target <- genotype.s[, target.SNPS]
imputed <- impute.snps(rules, target, as.numeric = FALSE)
print(imputed)


# ******************************************************************************

# Write genotype.s, geno.bim, genofile, phenotype, pcs, support & imputed 
# for future use.
# *Note this function can take a few minutes
save(genotype.s, geno.bim, genofile, phenotype, pcs, support, imputed,  
     file = working.data.fname(4))


# On to script 5.Data.analysis.R


#___ end _______________________________________________________________________
