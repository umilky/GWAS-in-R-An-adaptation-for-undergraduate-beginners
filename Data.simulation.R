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



