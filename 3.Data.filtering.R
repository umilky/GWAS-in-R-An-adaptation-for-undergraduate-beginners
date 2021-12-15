# === Data Quality Check & Filtering ===========================================

# *Note 1.Main. R is required to run first to access folder paths and to input 
# raw data required to run this script.
# All scripts must be run in numerical order 
# Load data created in previous script, 2.Reading.raw.data.R
load(working.data.fname(2))

# The following script performs quality checks and filtering of genomic data 
# using industry standard test such that fair analysis can be performed. 
# Refer to table 1 in the accompaning paper for a more detailed explanation. 

# 1) SNP level-filtering: Minor Allele Frequency and Call Rate
# 2) Sample level-filtering: Call Rate and Heterozygosity 
# 3) Sample level-filtering: Relatedness 
# 4) Sample level-filtering: Ancestry
# 5) SNP level-filtering: Hardy-Weinberg equilibrium


# === 1) SNP level-filtering: Minor Allele Frequency and Call Rate =============

# We want to filter out the common alleles from the SNPs. This 
# eliminates/reduces homogeneity from the sample and helps to draw statistically 
# significant relationship between the SNP and the trait (phenotype) we are 
# interested in. Another thing we use is the call rate. Call rate for a SNP 
# is the proportion of individuals in the study for which the corresponding SNP 
# information is not missing (i.e. no. of individuals with SNP information/total
# number of individuals in the study). These account for errors in genotyping. 

# Creating column summary statistics of SNP information 
snpsummary.col <- col.summary(genotype.s)
head(snpsummary.col)

# Setting thresholds for filtering 
call <- 0.95 # SNPs with 5% missing data is retained after filtering.
minor.f <- 0.01 

# Filtering SNP data by minor allele frequency(MAF) and call rate
fil.genotype <- with(snpsummary.col, (!is.na(MAF) & MAF > minor.f) 
                     & Call.rate >= call)
fil.genotype[is.na(fil.genotype)] <- FALSE #Remove NA values 

print(ncol(genotype.s)-sum(fil.genotype)) # 203287 SNPs will be removed

# Subset genotype.s for SNPs that meet the minor allele frequency(MAF)
# and call rate criteria.
genotype.s <- genotype.s[,fil.genotype]

# Subset the SNP summary data for SNPs those meet the MAF and call rate criteria
snpsummary.col <- snpsummary.col[fil.genotype,]
print(genotype.s) #658186 SNPs are remaining


# === 2) Sample level-filtering: Call Rate and Heterozygosity ==================

# Filtering to get rid of individuals with missing data and other criteria such 
# as racial, ethnic, or gender ambiguity (GWAS is population dependent). 
# We also get rid of contaminated samples (this varies)

# We will use call rate and heterozygosity as thresholds for filtering here. 
# Individuals who are missing genotype data for more than 5% of the
# typed SNPs are removed by the 95% call rate. Excess
# heterozygosity across typed SNPs within an individual may be an indication of 
# poor sample quality, while deficient heterozygosity can indicate inbreeding or 
# other substructure in that person. Thus, samples with an inbreeding 
# coefficient |F| = (1 - O/E) > 0.10 are removed, where O and E are respectively 
# the observed and expected counts of heterozygous SNPs within an individual.

# Creating row summary statistics for the sample 
snpsummary.row <- row.summary(genotype.s)

# Calculating the inbreeding coefficient for the sample and adding it to the 
# row summary statistics for the sample.
MAF <- snpsummary.col$MAF
callmatrix <- !is.na(genotype.s)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(snpsummary.row, Heterozygosity*(ncol(genotype.s))*Call.rate)
snpsummary.row$hetF <- 1-(hetObs/hetExp)

# Setting thresholds for sample filtering 
samp.call <- 0.95       # 95% call rate
het.cutoff <- 0.10      # Inbreeding coefficient cutoff

# Filtering sample based on call rate and heterozygosity
fil.sample <- with(snpsummary.row, !is.na(Call.rate) & Call.rate > samp.call &
                     abs(hetF) <= het.cutoff)
fil.sample[is.na(fil.sample)] <- FALSE   #removing NAs from the sample 

print(nrow(genotype.s)-sum(fil.sample))  # 0 subjects removed 

# Subset genotype and phenotype data for subjects meeting the heterozygosity and 
# call rate criteria.
genotype.s <- genotype.s[fil.sample,]
phenotype <- phenotype[rownames(genotype.s),]
head(genotype.s)


# === 3) Sample level-filtering: Relatedness  ==================================

# Now we will further filter the sample for relatedness. 
# We will use the linkage disequilibrium as a threshold to remove redundancy. 
# We will use the kinship threshold to remove related individuals from the 
# sample pool. 

# Setting thresholds for kinship and LD
kin.cutoff <- 0.1  # Kinship Cut-Off based on IBD coefficient
ld.cutoff <- 0.2   # LD cut-off. 0.2 seems to be the standard in GWAS. 

# Creating gds which are a combination of the PLINK files 
# SNPRelate functions require gds files to run
# Converting from PLINK to GDS
snpgdsBED2GDS(gwas.data$bed, gwas.data$fam, gwas.data$bim, gwas.data$gds)
genofile <- snpgdsOpen(gwas.data$gds, readonly = FALSE)

# Removing automatically added "-1" suffixes 
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = TRUE)

# Prune SNPs for IBD analysis 
# Set seed for reproducibility 
set.seed(1000)
# Create vector of sample ids from genomic data
genosample.ids <- rownames(genotype.s)

snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.cutoff,
                          sample.id = genosample.ids,
                          snp.id = colnames(genotype.s))
  
snpset.ibd <- unlist(snpSUB, use.names = FALSE)
cat(length(snpset.ibd), "will be used in IBD analysis\n") # expect 72890 SNP's

# Find IBD coefficients using Method of Moments procedure. 
# Method of Moments is a method of estimating population parameters 
# include pairwise kinship. 
ibd <- snpgdsIBDMoM(genofile, kinship = TRUE,
                    sample.id = genosample.ids,
                    snp.id = snpset.ibd,
                    num.thread = 1)
ibdcoeff <- snpgdsIBDSelection(ibd) # Pairwise sample comparison
head(ibdcoeff)

# Check if there are any candidates for relatedness
ibdcoeff <- ibdcoeff [ibdcoeff$kinship >= kin.cutoff, ]

# Iteratively remove samples with high kinship starting with the sample with 
# the most pairings. The statistical significance of a GWAS analysis is 
# impaired by a high degree of relatedness within the population. So, we want 
# to eliminate kinship as much as we possibly can. 
# create object for samples 
related.samples <- NULL

# Run iterative loop to add samples to related.samples from ibdcoeff
while (nrow(ibdcoeff) > 0 ) {
  # Count the number of occurrences of each and take the top one
  sample.counts <- arrange(count(c(ibdcoeff$ID1, ibdcoeff$ID2)), -freq)
  rm.sample <- sample.counts[1, "x"]
  cat("Removing sample", as.character(rm.sample), "too closely related to",
      sample.counts[1, "freq"], "other samples. \n")
  
  # Remove from ibdcoeff and add to list
  ibdcoeff <- ibdcoeff[ibdcoeff$ID1 !=rm.sample & ibdcoeff$ID2 != rm.sample,]
  related.samples <- c(as.character(rm.sample), related.samples)
}

# Subset genotype and phenotype to include only unrelated samples
genotype.s <- genotype.s[!(rownames(genotype.s) %in% related.samples),]
phenotype <- phenotype[!(phenotype$FamID %in% related.samples),]

genosample.ids <- rownames(genotype.s)

cat(length(related.samples), "similar samples removed due to 
    correlation coefficient >=", kin.cutoff, "\n")
print(genotype.s) # except all 1401 subject remain


# === 4) Sample level-filtering: Ancestry ======================================

# Checking for ancestry within sample. PCA is one approach to 
# visualizing and classifying  individuals into ancestry groups based on their 
# observed genetic makeup. This is done to eliminate sample-level errors. 
# For e.g. There could be an individual who do not fall within a racial/ethnic 
# cluster despite claiming to be of a certain racial/ethnic group.

# Find PCA matrix
pca.anc <- snpgdsPCA(genofile, sample.id = genosample.ids, snp.id = snpset.ibd, 
                 num.thread = 1)

# Create data frame of first two principal components
pctab <- data.frame(sample.id = pca.anc$sample.id,
                    PC1 = pca.anc$eigenvect[,1], # the first eigenvector
                    PC2 = pca.anc$eigenvect[,2], # the second eigenvector
                    stringsAsFactors = FALSE) 

# Plot the first two principal components
# Save plot as pdf in graphics folder 
pdf(paste(p.graphics, "pc.plot.pdf", sep = ""))
pc.plot <- plot(pctab$PC2, pctab$PC1, xlab = "Principal Component 2",
     ylab = "Principal Component 1", main = "Ancestry Plot")
dev.off()


# === 5) SNP level-filtering: Hardy-Weinberg equilibrium =======================

# Testing for violations of Hardy-Weinberg equilibrium (HWE).
# Violations of HWE can be an indication of genotyping error so it is common 
# practice to remove SNPs for which HWE is violated. 

# Setting threshold for HWE test p < 1x10^10-6
HWE.thres <- 10^-6  

# Filter based on population controls, in this case individuals without 
# coronary artery disease(CAD). 
CADcontrols <- phenotype[phenotype$CAD==0, 'FamID' ]
snpsum.col.cont <- col.summary(genotype.s[CADcontrols,])

# Filter based on HWE
HWE.use <- with(snpsum.col.cont, 
                !is.na(z.HWE) & ( abs(z.HWE) < abs( qnorm(HWE.thres/2))))

HWE.use[is.na(HWE.use)] <- FALSE          # Remove NA's as well
cat(ncol(genotype.s)-sum(HWE.use),"SNPs will be removed due to high HWE.\n") 

# Subset genotype and SNP summary data for SNPs that pass HWE criteria
genotype.s <- genotype.s[,HWE.use] #656890 SNPs remain
print(head(genotype.s)) 


# ******************************************************************************

# Write genotype.s, geno.bim, genofile & phenotype for future use.
# *Note this function can take a few minutes
save(genotype.s, geno.bim, genofile, phenotype, file = working.data.fname(3))

# on to script 4.Data.generation.R

#___ end _______________________________________________________________________
