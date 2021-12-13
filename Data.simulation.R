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
