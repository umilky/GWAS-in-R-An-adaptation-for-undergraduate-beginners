# === Data analysis ============================================================

## installing packages 
# They are in our GitHub Repository 

### we could create a folder called "packages" and add them in there?

# association analysis of typed SNPs 
# this is specific to the data and here we use the HDL-cholestrol as the trait

# merge clinical data and principal components to create phenotype table

pheno.sub <- merge.data.frame(phenotype,pcs) 
# data frame => [FamID CAD sex age hdl...pc10]

# rank-based inverse normal transformation of hdl
pheno.sub$phenotype <- rntransform(pheno.sub$hdl, family = "gaussian")

# show that the assumptions of normality met after transformation
par(mfrow=c(1,2))
hist(pheno.sub$hdl, main = "Histogram of HDL", xlab = "HDL")
hist(pheno.sub$phenotype, main = "Histogram of Transformed HDL",
     xlab = "Transformed HDL")

# remove columns that are unnecessary from the table 
pheno.sub$hdl <- NULL
pheno.sub$ldl <- NULL
pheno.sub$tg <- NULL 
pheno.sub$CAD <- NULL 


# Rename columns to match names necessary for GWAS() function 
names(pheno.sub)[names(pheno.sub ) == "FamID"] <- "id"

# Include only subjects with hdl data 
pheno.sub <- pheno.sub[!is.na(pheno.sub$phenotype),]
# 1309 subjects identified to have phenotype data

print (head(pheno.sub))


# Run GWAS analysis using parallel processing 

# loading required packages 
#### we have this in 1.main 
#library(plyr)

# Loading GWAS function. This function was created by Reed et al. (2015).
# Gives output as a .txt file 

source("GWAA.R")

# The GWAA function takes an hour or two to run depending on data and your computer
# Don't get discouraged if you don't see anything for a bit.

Start <- Sys.time()

GWAA(genodata = genotype.s, phenodata = pheno.sub, filename = gwaa.out)

End <- Sys.time()
print(End-Start)

# Reading GWAS output 
GWAS.output <- read.table(gwaa.out, sep = "",header=TRUE, 
                         colClasses=c("character"))
head(GWAS.output)

# Calculating the -log10 of the p-values 
GWAS.output$p.value <- as.numeric(GWAS.output$p.value)
GWAS.output$neg.logp <- -log10(GWAS.output$p.value)
head(GWAS.output)


# Merge output with geno.bim by SNP name to add position and chromosome number
# this does not work and kills the whole thing :(
GWAS.output <- merge(GWAS.output, geno.bim[,c("SNP", "chr", "position")], 
                    by = "SNP")
head(GWAS.output)
#rm(geno.bim)

# Order SNPs by significance 
GWAS.output <- arrange(GWAS.output, neg.logp)
print(head(GWAS.output))


## Association analysis of imputed SNPs 

# Performing association testing for the imputed SNPs

rownames(pheno.sub) <- pheno.sub$id

imp <- snp.rhs.tests(phenotype ~ sex + age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 +
                       pc7 + pc8 + pc9 + pc10, family = "Gaussian", 
                     data = pheno.sub, snp.data = target, rules = rules)

# Obtain p values for imputed SNPs
imp.results <- data.frame(SNP = imp@snp.names, p.value = p.value(imp), 
                      stringsAsFactors = FALSE)
imp.results <- imp.results[!is.na(imp.results$p.value),]

# Writing a file containing the results 
# write.csv(results, impute.out.fname, row.names = FALSE)

# Merge imputation testing results with support 
impute.out <- merge(imp.results, support[,c("SNP", "position")])
impute.out$chr <- 16

impute.out$type <- "imputed"

# Calculating the -log10 of the imputed p-values 
impute.out$neg.logp <- -log10(impute.out$p.value)

# Order by p-values 
impute.out <- arrange(impute.out, p.value)
print(head(impute.out))

# map2gene function 
# Returns the subset of SNPs that are within extend.boundary of gene
# using the coords table of gene locations
map2gene <- function(gene, coords, SNPs, extend.boundary = 5000) {
  coordsSub <- coords[coords$gene == gene,] 
  #Subset coordinate file for spcified gene
  coordsSub$start <- coordsSub$start - extend.boundary # Extend gene boundaries
  coordsSub$stop <- coordsSub$stop + extend.boundary
  SNPsub <- SNPs[SNPs$position >= coordsSub$start 
                 & SNPs$position <= coordsSub$stop &
                   SNPs$chr == coordsSub$chr,] #Subset for SNPs in gene
  return(data.frame(SNPsub, gene = gene, stringsAsFactors = FALSE))
}

#Reading file containing protein coding gene coordinates 
genes <- read.csv(protein.coding.coords.fname, stringsAsFactors = FALSE)

#Subset for CETP SNPs
impCETP <- map2gene("CETP", coords = genes, SNPs = impute.out)

# Filtering the imputed CETP SNP genotypes
impCETPgeno <- imputed[, impCETP$SNP]


## Integrating typed and imputed SNPs 

GWAS.output$type <- "typed"

GWAS.comb <- rbind.fill(GWAS.output, impute.out)
head(GWAS.comb)
tail(GWAS.comb)
str(GWAS.comb)

#Subset for CETP SNPs 
typCETP <- map2gene("CETP", coords = genes, SNPs = GWAS.output)

#Combine CETP SNPs for typed and imputed analysis 

CETP <- rbind.fill(typCETP, impCETP)[,c("SNP", "p.value", "neg.logp", "chr",
                                        "position", "type", "gene")]

print(CETP)
# on to script 6.Visualization 

#___ end _______________________________________________________________________
