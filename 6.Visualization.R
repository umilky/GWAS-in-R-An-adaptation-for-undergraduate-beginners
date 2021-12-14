# === Visualization ============================================================
# Scripts 1-5 are required to run first to access folder paths and obtain the 
# data required for creating plots.
# This file creates plots from the analysed data.
# Plots are saved as pdfs to c.graphics folder.

## Visualizing the results of the GWAS analysis 

# Manhattan plots are used to visualize GWA significance level by chromosome
# location. Visual inspection of this plot allows for identification of SNPs 
# with relatively small p-values that are in regions with relatively large and 
# non-significant p-values, suggesting potentially false findings.
# Multiple signals in the CETP region suggest that this may be a true signal.



# creating Manhattan Plot and saving it as a .pdf. 
# we use a function created by Reed et. al (2015) to generate this plot. 

# loading Manhattan Plot function 
source("Manhattan.plot.R")

pdf(paste(p.graphics, "ManhattanPlot.pdf", sep = ""))
par(mfrow = c(1,1))
manhattan.plot.s <- GWAS.manhattan(GWAS.comb)
dev.off()

 
# Q-Q plots are used to visualize the relationship between the expected
# and observed distributions of SNP-level test statistics. We use the function 
# estlambda from the package GenABEL to generate our QQ plots. 

GWAS.output$t.value <- as.numeric(GWAS.output$t.value)

# creating QQ Plot and saving it as a .pdf
pdf(paste(p.graphics, "ManhattanPlot.pdf", sep = ""))
lambda.plot.s <- estlambda(GWAS.output$t.value^2, plot = TRUE, 
                           method = "median")
dev.off()

#___ end _______________________________________________________________________




