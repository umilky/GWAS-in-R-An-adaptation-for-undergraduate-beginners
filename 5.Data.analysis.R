# === Data Analysis ============================================================

# *Note 1.Main. R is required to run first to access folder paths and to input 
# raw data required to run this script.
# All scripts must be run in numerical order 
# Load data created in previous script, 4.Data.Generation.R
load(working.data.fname(4))

# The following script performs a genome-wide association analysis(GWAA) of 
# typed SNPs this is specific to the data and here we use the HDL-cholestrol 
# as the trait. The results of the GWAA are then sorted into a data frame for 
# visualization in the next script. An association analysis is  also run for the 
# imputed SNPs which are also sorted into a data frame for further analysis. 

# 1) Manipulation of Phenotype Data
# 2) Run GWAA
# 3) GWAA Object Creation
# 4) Association analysis of imputed SNPs 
# 5) Integrating typed and imputed SNPs

# === 1) Manipulation of Phenotype Data ========================================

# Creating a suitable table to run a genome-wide association analysis(GWAA) that
# includes phenotype data and principle components. The data is also visualized 
# to check for normality and adjusted accordingly. 

# Merge phenotype data and principal components to create phenotype table
pheno.sub <- merge.data.frame(phenotype,pcs) 
# data frame => [FamID CAD sex age hdl...pc10]

# Rank-based inverse normal transformation of HDL
pheno.sub$phenotype <- rntransform(pheno.sub$hdl, family = "gaussian")

# Show that the assumptions of normality met after transformation
par(mfrow=c(1,2))
hist(pheno.sub$hdl, main = "Histogram of HDL", xlab = "HDL")
hist(pheno.sub$phenotype, main = "Histogram of Transformed HDL",
     xlab = "Transformed HDL")

# Remove columns that are unnecessary from the table 
# Age and sex are left in table as they established covariates for HDL level 
pheno.sub$hdl <- NULL
pheno.sub$ldl <- NULL
pheno.sub$tg <- NULL 
pheno.sub$CAD <- NULL 

# Rename columns to match names necessary for GWAS() function 
names(pheno.sub)[names(pheno.sub ) == "FamID"] <- "id"

# Include only subjects with HDL data 
pheno.sub <- pheno.sub[!is.na(pheno.sub$phenotype),]
# 1309 subjects identified to have phenotype data
print (head(pheno.sub))


# === 2) Run GWAA ==============================================================

# Use GWAA function to run analysis on merged data frame 

# Loading GWAS function. This function was created by Reed et al. (2015).
# Gives output as a .txt file 
source("GWAA.R")

# The GWAA function takes an hour or two to run depending on the size of 
# the data frame and your computer don't get discouraged 
# if you don't see anything for a bit.

Start <- Sys.time()

GWAA(genodata = genotype.s, phenodata = pheno.sub, filename = gwaa.out)

End <- Sys.time()

# Inform user of run time 
print(End-Start)

# The GWAA file can be found in the folder "c.GWAS" 
# under the file name "GWASout.txt".


# === 3) GWAA Object Creation ==================================================

# Reading "GWASout.txt" file to create data frame with results from GWAA.
# A negative log of the p-values is calculated to compare against the 
# Bonferroni corrected threshold and the less stringent suggestive 
# association threshold. 
# The data frame will be used to visualize a Manhattan plot in the next script. 

# Reading GWAA output 
GWAS.output <- read.table(gwaa.out, sep = "",header=TRUE, 
                         colClasses=c("character"))
head(GWAS.output)

# Calculating the -log10 of the p-values 
GWAS.output$p.value <- as.numeric(GWAS.output$p.value)
GWAS.output$neg.logp <- -log10(GWAS.output$p.value)
head(GWAS.output)

# Merge output with geno.bim by SNP name to add position and chromosome number
GWAS.output <- merge(GWAS.output, geno.bim[,c("SNP", "chr", "position")], 
                    by = "SNP")
head(GWAS.output)

# Tagging SNPs by type 
GWAS.output$type <- "typed"

# Data check to make sure plots can run correctly. 
GWAS.output$t.value <- as.numeric(GWAS.output$t.value)

# Order SNPs by significance 
GWAS.output <- arrange(GWAS.output, neg.logp)
print(head(GWAS.output))


# === 4) Association analysis of imputed SNPs ==================================

# Association analysis of imputed SNPs to be included in Manhattan plot. 
# Imputed SNPs will be compared to the same threshold as the typed SNPs 

# Setting row names for pheno.sub 
rownames(pheno.sub) <- pheno.sub$id

# Performing association testing for the imputed SNPs
imp <- snp.rhs.tests(phenotype ~ sex + age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 +
                       pc7 + pc8 + pc9 + pc10, family = "Gaussian", 
                     data = pheno.sub, snp.data = target, rules = rules)

# Obtain p values for imputed SNPs
imp.results <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), 
                      stringsAsFactors = FALSE)
imp.results <- imp.results[!is.na(imp.results$p.value),]

# Merge imputation testing results with support 
impute.out <- merge(imp.results, support[,c("SNP", "position")])

# Subsetting for chromosome 16 based on available 1000 Genomes data
impute.out$chr <- 16

# Tagging SNPs by type 
impute.out$type <- "imputed"

# Calculating the -log10 of the imputed p-values 
impute.out$neg.logp <- -log10(impute.out$p.value)

# Order by p-values 
impute.out <- arrange(impute.out, p.value)
print(head(impute.out))


# === 5) Integrating typed and imputed SNPs  ==================================

# Merging of typed and imputed SNPs for visualization in next script. 

# Combination of typed and imputed data frames 
GWAS.comb <- rbind.fill(GWAS.output, impute.out)
head(GWAS.comb)
tail(GWAS.comb)
str(GWAS.comb)


# ******************************************************************************

# Write GWAS.output & GWAS.comb for future use. 
# *Note this function can take a few minutes
save(GWAS.output, GWAS.comb, file = working.data.fname(5))
# on to script 6.Visualization 

#___ end _______________________________________________________________________
