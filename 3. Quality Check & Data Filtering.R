# === Quality Check & Filtering ================================================

## SNP level-filtering 
# We want to filter out the common alleles from the SNPs. This 
# eliminates/reduces homogeneity from the sample and helps to draw statistically 
# significant relationship between the SNP and the trait (phenotype) we are 
# interested in. Another thing we use is the call rate. Call rate for a SNP 
# is the proportion of individuals in the study for which the corresponding SNP 
# information is not missing (i.e. no. of individuals with SNP information/total
# number of individuals in the study). These account for errors in genotyping. 

# Creating summary statistics of SNP information 

snpsummary <- col.summary(genotype.s)
head(col.summary)

# Setting thresholds for the filtering 

call <- 0.95 # SNPs with 5% missing data is retained after filtering.
minor.f <- 0.01 

# Filtering the SNP data 
fil.genotype <- with(snpsummary, (!is.na(MAF) & MAF > minor.f) 
                     & Call.rate >= call)
fil.genotype[is.na(fil.genotype)] <- FALSE #Remove NA values 

print(ncol(genotype.s)-sum(fil.genotype)) # 203287 SNPs will be removed

#Subset genotype.s for SNPs those meet the MAF and call rate criteria

genotype.s <- genotype.s[,fil.genotype]

# Subset the SNP summary data for SNPs those meet the MAF and call rate criteria

snpsummary <- snpsummary[fil.genotype,]

print(genotype.s) #658186 SNPs are remaining


## Sample level filtering 

# Want to get rid of individuals with missing data and other criteria such as 
# racial, ethnic, or gender ambiguity (GWAS is population dependent)
# We also get rid of contaminated samples (this varies)

# We will use call rate and heterozygosity as thresholds for filtering here. 
# Individuals who are missing genotype data for more than 5% of the
# typed SNPs are removed by the 95% call rate. Excess
# heterozygosity across typed SNPs within an individual may be an indication of 
# poor sample quality, while deficient heterozygosity can indicate inbreeding or 
# other substructure in that person. Thus, samples with an inbreeding 
# coefficient |F| = (1 - O/E) > 0.10 are removed, where O and E are respectively 
# the observed and expected counts of heterozygous SNPs within an individual.

# Creating summary statistics for the sample 
sample.summary <- row.summary(genotype.s)

gc()
rm(clinical)

# Calculating the inbreeding coefficient for the sample 

MAF <- snpsummary$MAF 

Exphet <- (!is.na(genotype.s)) %*% (2*MAF*(1-MAF))
   # Exphet is calculated using 2*p*(1-p), where p is the dominant allele
   # frequency at that SNP


Obshet <- with(sample.summary, Heterozygosity*(ncol(genotype.s))*Call.rate)

F.value <- 1-(Obshet/Exphet)

#Including the F-value in the sample statistics 

sample.summary$hetF <- F.value 

# Setting thresholds for sample filtering 
samp.call <- 0.95       # 95% call rate
het.cutoff <- 0.10      # Inbreeding coefficient cutoff


#Filtering sample based on threshold

fil.sample <- with(sample.summary, !is.na(Call.rate) & Call.rate > samp.call &
                     abs(hetF) <= het.cutoff)

fil.sample[is.na(fil.sample)] <- FALSE   #removing NAs from the sample 

print(nrow(genotype.s)-sum(fil.sample))  # 0 subjects removed 

# Subset genotype and phenotype data for subjects meeting the heterozygosity and 
# call rate criteria 

genotype.s <- genotype.s[fil.sample,]
phenotype <- phenotype[rownames(genotype.s),]


# Now we will further filter the sample for relatedness and redundancy. 
# We will use the linkage disequilibrium as a threshold to remove redundancy. 
# We will use the kinship threshold to remove related individuals from the 
# sample pool. 

kin.cutoff <- 0.1  # Kinship Cut-Off based on IBD coefficient
ld.cutoff <- 0.2   # LD cut-off. 0.2 seems to be the standard in GWAS. 


# Creating the gds files that are required for the SNPRelate functions 
# Don't understand why we need to do this, just following the instructions in 
# the GWAS tutorials 

# Converting from PLINK to GDS
snpgdsBED2GDS(gwas.data$bed, gwas.data$bim, gwas.data$fam, gwas.data$gds)
genofile <- snpgdsOpen(gwas.data$gds, readonly = FALSE)

# Removing automatically added "-1" suffixes 
gds.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
gds.ids <- sub("-1", "", gds.ids)
add.gdsn(genofile, "sample.id", gds.ids, replace = T)

#Prune SNPs for IBD analysis 
genosample.ids <- rownames(genotype.s)

snpSUB <- snpgdsLDpruning(genofile, ld.threshold = ld.cutoff,
                          sample.id = genosample.ids,
                          snp.id = colnames(genotype.s))
  



