### ========= Data Simulation Using PhenotypeSimulator =========================

# install required packages 
install.packages("PhenotypeSimulator")
library("PhenotypeSimulator")

n <-  1001
# simulate phenotype data
simulatephenotype <- function(n, variable, m, st) {
  id <- c(1:n)
  d.pres <- sample(c(0,1), replace=TRUE, size=n)
  sex <- sample(c(1,2), replace=TRUE, size=n)
  age <- sample(c(18:65), replace=TRUE, size=n)
  v <- round(rnorm(n, m, st))
  df <- cbind.data.frame(id,d.pres,sex,age,v)
  return(df)
}

# simulation of population blood glucose levels in mg/dl
sim.phen <- simulatephenotype(n, "blood glucose (mg/dl)", 96, 21.6)

# 
write.csv(sim.phen, paste(p.simulated.data, "sim.phen.csv", sep = ""), row.names = FALSE)

### ========= Data Simulation Using GenIO ======================================

## Installing Packages 

install.packages("genio")

library(genio)

## Creating large genotype matrix 

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

## making bim file 

bim <- make_bim( n = n_loci )

# Simulating random SNPs 
#bim$id <- round(runif(n_loci, min = 1, max = 1000000))

# Making SNP IDs look like "rs" IDs
bim$id <- paste0('rs', bim$id)

# Simulating random chromosomes between 1 and 22
bim$chr <- round(runif(n_loci, min = 1, max = 22))

# Simulating random positons 
bim$pos <- round(runif(n_loci, min = 10000, max =99999 ))

# Select randomly between Cs and Gs for the reference alleles
bim$ref <- sample(c('C', 'G'), n_loci, replace = TRUE)

# Select randomly between As and Ts for the alternative alleles
bim$alt <- sample(c('A', 'T'), n_loci, replace = TRUE)

## Creating the fam file 
fam <- make_fam( n = n_ind )

# Adding the sex values. They are usually 1 and 2
fam$sex <- sample(1:2, n_ind, replace = TRUE)

# Making phenotype continuous. 
fam$pheno <- rnorm(n_ind)


# Inspecting x. 
x[1:10, 1:10]


# Trying out the geno_to_char function 

geno.c <- geno_to_char(x, bim)

geno.c[1:10, 1:10]

## Writing PLINK files to match the 'template'
# This is probably not working- I cannot find the data anymore??

write_bed("simulated.bed",x)
write_bim("simulated.bim",bim)
write_fam("simulated.fam",fam)





