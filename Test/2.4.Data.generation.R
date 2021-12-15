# === Data Generation ==========================================================

# *Note 2.1.Main. R is required to run first to access folder paths and to input 
# raw data required to run this script.
# All scripts must be run in numerical order 
# Load data created in previous script, 3.Data.filtering.R
load(working.data.fname(2.3))

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
pcs <- data.frame(id = pca$sample.id, pca$eigenvect[,1:10],
                  stringsAsFactors = FALSE)
colnames(pcs)[2:11] <- paste("pc", 1:10, sep ="")

print(head(pcs))


# === 2) Imputation of non-typed SNPs ==========================================

# Genotypes of untyped SNPs may have a functional relationship to the 
# outcome and therefore provide additional power for identifying association.
# Imputation of SNPs involve estimating missing genotypes from the haplotype 
# or genotype reference panel. Here, we use the 1000 Genomes data for this. 

# Unfortunately the test data is from mice such that there are no matches
# for the supplied chromosome 16 file. There will be no imputed SNPs on
# the Manhattan plot for the test data. 


# ******************************************************************************

# Write genotype.s, geno.bim, genofile, phenotype, pcs, support & imputed 
# for future use.
# *Note this function can take a few minutes
save(genotype.s, geno.bim, genofile, phenotype,
     pcs, file = working.data.fname(2.4))


# On to script 2.5.Data.analysis.R


#___ end _______________________________________________________________________

