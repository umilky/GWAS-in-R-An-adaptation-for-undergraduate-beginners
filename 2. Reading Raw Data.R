# === Reading Data in R ========================================================

# Read PLINK files to create a list
# PLINK is a genome association analysis toolkit
# bed files contain individual genotypes
# bin files contain the locations of all SNPs in the data in the genome
# fam files contain information about the individual including family pedigrees

geno.me <- read.plink(gwas.data$bed, gwas.data$bim, gwas.data$fam)

#creating the genotyope object from the geno.me list 

genotype <- geno.me$genotype 
print(genotype)
head(genotype)


#creating SNP object from the geno.me list 
SNP.genome <- geno.me$map
colnames(SNP.genome) <- c("chromosome", "SNP", "gen.dist", "position", "N1", 
                          "N2")
print(head(SNP.genome))

#removing raw data to create space 
rm(geno.me)

#reading phenotype data from the clinical file 
#apparently using colClasses saves data importing time
#Also, cannot subset genotype without rearranging using colClasses first
phenotype <- read.csv(clinical.data,
                     colClasses = c("character", "factor", "factor", 
                                    rep("numeric", 4)))
rownames(phenotype) <- phenotype$FamID
head(phenotype)

#subsetting genotype to only include individuals with known phenotype 

genotype.s <- genotype[phenotype$FamID, ]

#Now we move on to the next script to filter this data

