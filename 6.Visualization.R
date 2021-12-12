# === Visualization ============================================================

## Visualizing the results of the GWAS analysis 

install.packages("survival")
# Manhattan plots are used to visualize GWA significance level by chromosome
# location. Visual inspection of this plot allows for identification of SNPs 
# with relatively small p-values that are in regions with relatively large and 
# non-significant p-values, suggesting potentially false findings.
# Multiple signals in the CETP region suggest that this may be a true signal.
source("Manhattan.plot.R")

#####I dont have mfrow ####### doesnt work
par(mfrow(1,1))


# creating Manhattan Plot. 
# need to save as .pdf later


### this function wont work with the data as a dataframe? it has to be a vector?
head(GWASoutput)
GWAS_Manhattan(GWASoutput)



# creating QQ Plot 
# Q-Q plots are used to visualize the relationship between the expected
# and observed distributions of SNP-level test statistics. We use the function 
# estlambda from the package GenABEL to generate our QQ plots. 

#creating QQ plot. Save as pdf later

### lambda plot shows invalid plot index? but then it shows up later, weird 
### also takes a really long time to run, but then shows an ok graph I guess
GWASoutput$t.value <- as.numeric(GWASoutput$t.value)
lambda <- estlambda(GWASoutput$t.value^2, plot = TRUE, method = "median")

# END? 


