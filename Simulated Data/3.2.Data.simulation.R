# === Data Simulation for GWAS =================================================

# This script simulates genomic and phenotype data required for GWAS analysis.
# PLINK files are created for the genomic data and the a csv file is created
# for the phenotype file. 

# 1) Phenotype simulation
# 2) Genomic data simulation


# === 1) Phenotype simulation ==================================================

# Number of individuals in sample 
n <-  1001
# Simulate phenotype data function
simulatephenotype <- function(n, variable, m, st) {
  # Individual identification number
  id <- c(1:n)
  # Disease presence for sample control, 0 or 1 
  d.pres <- sample(c(0,1), replace=TRUE, size=n)
  # Sex, 0 or 1
  sex <- sample(c(1,2), replace=TRUE, size=n)
  # sample age between 18 and 65
  age <- sample(c(18:65), replace=TRUE, size=n)
  # Variable generation bound by input
  v <- round(rnorm(n, m, st))
  df <- cbind.data.frame(id,d.pres,sex,age,v)
  return(df)
}

# Simulation of population blood glucose levels in mg/dl. 
# Based on population data from Chengdu, China (Huang et al., 2017)
sim.phen <- simulatephenotype(n, "blood glucose (mg/dl)", 96, 21.6)
# Renaming variable column 
names(sim.phen)[names(sim.phen) == "v"] <- "blood.glucose.mg.dl)"

# Write csv file for simulated phenotype data, saved to a.simulated.data folder
write.csv(sim.phen, paste(p.simulated.data, "sim.phen.csv", sep = ""),
          row.names = FALSE)

# === 2) Genomic data simulation ===============================================

# Data Simulation Using GenIO to create a large genotype matrix for analysis. 
# We randomly generate binary genotype data and use the parameters to 
# create PLINK readable files to use in our template.

# Data dimensions.

# Number of loci.
n_loci <- 10001

# Number of individuals
n_ind <- 1001

# Overall allele frequency
p <- 0.5

# Rate of missing values
miss <- 0.1

# Total number of genotypes
n_geno <- n_ind * n_loci

# Draw random genotypes from Binomial
x <- rbinom( n_geno, 2, p)

# Add missing values
x[ sample(n_geno, n_geno * miss) ] <- NA

# Turn into matrix
x <- matrix(x, nrow = n_loci, ncol = n_ind)

# Creating bim file 
bim <- make_bim( n = n_loci )

# Making SNP IDs look like "rs" IDs
bim$id <- paste0('rs', bim$id)

# Simulating random chromosomes between 1 and 22
bim$chr <- round(runif(n_loci, min = 1, max = 22))

# Simulating random positions 
bim$pos <- round(runif(n_loci, min = 10000, max =99999 ))

# Select randomly between Cs and Gs for the reference alleles
bim$ref <- sample(c('C', 'G'), n_loci, replace = TRUE)

# Select randomly between As and Ts for the alternative alleles
bim$alt <- sample(c('A', 'T'), n_loci, replace = TRUE)

## Creating fam file 
fam <- make_fam( n = n_ind )

# Adding the sex values. They are usually 1 and 2
fam$sex <- sample(1:2, n_ind, replace = TRUE)

# Making phenotype continuous. 
fam$pheno <- rnorm(n_ind)

# Inspecting x. 
x[1:10, 1:10]

# Writing PLINK files required for GWAS analysis 
write_bed("simulated.bed",x)
write_bim("simulated.bim",bim)
write_fam("simulated.fam",fam)

# ***IMPORTANT 
# *Note these files must be manually moved to the a.simulated.data folder.
# The authors of this code did not establish a path such that manual movement of 
# the files is required by the user. 


# on to script 3.3.Reading.raw.data.R

#___ end _______________________________________________________________________






