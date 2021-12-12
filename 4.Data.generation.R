# === Data Generation ==========================================================

## Principal Components (PCs) based on observed genotype data capture 
## information on substructure - genetic diversity in an apparently homogenous 
## population that is caused by population genetic history (migration, 
## selection, and/or ethnic integration). These substructures cannot be inferred 
## from self-reported race and ethnicity variables.


## generating principal components for modeling 

# set LD threshold to 0.2
ld.cutoff <- 0.2

set.seed(1000)
genosample.ids <- rownames(genotype.s)
snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.cutoff,
                          sample.id = genosample.ids, # only analyze filtered samples
                          snp.id = colnames(genotype.s)) # only analyze filtered samples
snpset.pca <- unlist(snpSUB, use.names = FALSE)
cat(length(snpset.pca), "\n") #72665 SNPs will be used in PCA analysis

pca <- snpgdsPCA(genofile, sample.id = genosample.ids,
                 snp.id = snpset.pca, num.thread = 1)

## find and record first 10 principal components
## pcs will be a N:10 matrix. Each column is a principal component. 
pcs <- data.frame(FAMID = pca$sample.id, pca$eigenvect[,1:10],
                  stringsAsFactors = FALSE)
colnames(pcs)[2:11] <- paste("pc", 1:10, sep ="")
# rename "FAMID" to "FamID" to match other files
names(pcs)[names(pcs) == "FAMID"] <- "FamID"

print(head(pcs))

## Imputation of non-typed SNPs (honestly should get rid of this; adds work; not mandatory)
## genotypes of untyped SNPs may have a functional relationship to the 
## outcome and therefore provide additional power for identifying association.
## Imputation of SNPs involve estimating missing genotypes from the haplotype 
## or genotype reference panel. Here, we use the 1000 Genomes data for this. 


## read in 1000g data for given chromosome 16 
thougeno <- read.pedfile(onethou.fn$ped, snps = onethou.fn$info, which = 1)

## obtain genotype data for given chromosome
geno.matrix <- thougeno$genotypes

## obtain the chromosome position for each SNP
support <- thougeno$map
colnames(support) <- c("SNP", "position", "N1", "N2")
head(support)

## imputation of non-typed 1000g SNPs
## subset for SNPs on given chromosome
pres.SNPS <- colnames(genotype.s)
pres.dat.chr <- geno.bim[geno.bim$SNP %in% pres.SNPS & geno.bim$chr==16,]
target.SNPS <- pres.dat.chr$SNP

## subset 1000g data for our SNPs
## "missing" and "present are snpMatrix objects needed for imputation rules
is.present <- colnames(geno.matrix) %in% target.SNPS
missing <- geno.matrix[,!is.present]
print(missing)

present <- geno.matrix[, is.present]
print(present)

## obtain positions of SNPs to be used for imputation rules
pos.present <- support$position[is.present]
pos.missing <- support$position[!is.present]

## calculate and store the imputation rules using snp.imputation()
rules <- snp.imputation(present, missing, pos.present, pos.missing)

# remove failed imputations 
rules <- rules[can.impute(rules)]
cat("Imputation rules for", length(rules), "SNPs were estimated\n")

## quality control for imputation certainty and MAF
# set thresholds
r2threshold <- 0.7
minor <- 0.01

# filter on imputation certainty and MAF
rules <- rules[imputation.r2(rules) >= r2threshold]

cat(length(rules), "imputation rules remain after uncertain imputations were removed\n")
## 162628 imputation rules remain after uncertain imputation were removed

rules <- rules[imputation.maf(rules) >= minor]
cat(length(rules), "imputation rules remain after MAF filtering\n")

## obtain posterior expectation of genotypes of imputed SNPs
target <- genotype.s[, target.SNPS]
imputed <- impute.snps(rules, target, as.numeric = FALSE)
print(imputed)

# on to script 5.Data.analysis






