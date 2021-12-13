### ========= Data Simulation Using PhenotypeSimulator =========================

install.packages("PhenotypeSimulator")

library("PhenotypeSimulator")
### STILL FIGURING OUT ###
### via `runSimulation`

# simulate phenotype with population structure and observational noise effects
# only
# genetic variance
genVar <- 0.4
# random genetic variance: h2b
phenotype <- runSimulation(N = 100, P = 15, tNrSNP = 10000, SNPfrequencies = c(0.05, 0.1,0.3,0.4),
                           genVar = 0.4, h2bg = 1, phi = 1,
                           verbose = TRUE, nonlinear="exp",
                           proportionNonlinear = 0.7)
head(phenotype)
str(phenotype)
head(phenotype$phenoComponentsFinal$Y)


# step by step 
# doesnt work 
# Set parameters
genVar <- 0.4
noiseVar <- 1- genVar 
shared <- 0.6
independent <- 1 - shared
# simulate simple bi-allelic genotypes and estimate kinship
genotypes <- simulateGenotypes(N = 100, NrSNP = 10000,
                               frequencies = c(0.05, 0.1, 0.3, 0.4),
                               verbose = FALSE) genotypes_sd <- standardiseGenotypes(genotypes$genotypes)
2
kinship <- getKinship(N=100, X=genotypes_sd, verbose = FALSE)
# simulate phenotype components
genBg <- geneticBgEffects(N = 100, kinship = kinship, P = 15) noiseBg <- noiseBgEffects(N = 100, P = 15)
# rescale phenotype components
genBg_shared_scaled <- rescaleVariance(genBg$shared, shared * genVar) genBg_independent_scaled <- rescaleVariance(genBg$independent,
                                                                                                                  independent * genVar) noiseBg_shared_scaled <- rescaleVariance(noiseBg$shared, shared * noiseVar)
noiseBg_independent_scaled <- rescaleVariance(noiseBg$independent, independent *noiseVar)
# Total variance proportion shave to add up yo 1
total <- independent * noiseVar + independent * genVar + shared * noiseVar + shared * genVar
total == 1 #> [1] TRUE
# combine components into final phenotype
Y <- scale(genBg_shared_scaled$component + noiseBg_shared_scaled$component + genBg_independent_scaled$component + noiseBg_independent_scaled$component)
# transform phenotype non-linearly
Y_nl <- transformNonlinear(Y, alpha=0.7, method="exp") 

# Use exp as transformation method



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

# Adding the "chr" prefix to the chromosome values for recognition purposes
bim$chr <- paste0('chr', bim$chr)

# Making SNP IDs look like "rs" IDs
bim$id <- paste0('rs', bim$id)

# Making positions 1000 bigger
bim$pos <- bim$pos * 1000

# Select randomly between Cs and Gs for the reference alleles
bim$ref <- sample(c('C', 'G'), n_loci, replace = TRUE)

# Select randomly between As and Ts for the alternative alleles
bim$alt <- sample(c('A', 'T'), n_loci, replace = TRUE)

## Creating the fam file 
fam <- make_fam( n = n_ind )

# Add prefixes to families and IDs for recognition purposes 
fam$fam <- paste0('fam', fam$fam)
fam$id <- paste0('id', fam$id)

# Adding the sex values. They are usually 1 and 2
fam$sex <- sample(1:2, n_ind, replace = TRUE)

# Making phenotype continuous. 
fam$pheno <- rnorm(n_ind)

# Add column and row names from bim and fam tables we just created.
rownames(x) <- bim$id
colnames(x) <- fam$id

# Inspecting x. 
x[1:10, 1:10]


# Trying out the geno_to_char function 

geno.c <- geno_to_char(x, bim)

geno.c[1:10, 1:10]

## Writing PLINK files to match the 'template'
# This is probably not working- I cannot find the data anymore??

file_plink <- tempfile('vignette-random-data')

time_write_genio <- system.time(
  write_plink(file_plink,x, bim, fam)
)






