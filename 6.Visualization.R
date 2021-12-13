# === Visualization ============================================================

## Visualizing the results of the GWAS analysis 

# Manhattan plots are used to visualize GWA significance level by chromosome
# location. Visual inspection of this plot allows for identification of SNPs 
# with relatively small p-values that are in regions with relatively large and 
# non-significant p-values, suggesting potentially false findings.
# Multiple signals in the CETP region suggest that this may be a true signal.

source("Manhattan.plot.R")
par(mfrow = c(1,1))

# creating Manhattan Plot. 
# need to save as .pdf later
manhattan.plot.s <- GWAS.manhattan(GWAS.comb)
manhattan.plot.s


# creating QQ Plot 
# Q-Q plots are used to visualize the relationship between the expected
# and observed distributions of SNP-level test statistics. We use the function 
# estlambda from the package GenABEL to generate our QQ plots. 

#creating QQ plot. Save as pdf later
GWAS.output$t.value <- as.numeric(GWAS.output$t.value)
lambda.plot.s <- estlambda(GWAS.output$t.value^2, plot = TRUE, method = "median")
lambda.plot.s
# END? 


