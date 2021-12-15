# === Visualization ============================================================

# *Note 1.Main. R is required to run first to access folder paths and to input 
# raw data required to run this script.
# All scripts must be run in numerical order 

# Load data created in previous script, 5.Data.analysis.R
load(working.data.fname(2.5))

# Visualizing the results of the GWAS analysis. 
# This script creates a Manhattan plot and a QQ plot from the analysed data.
# Plots are saved as pdfs to c.graphics folder.

# 1) Manhattan plot
# 2) Q-Q plot

# === 1) Manhattan plot ========================================================

# Manhattan plots are used to visualize GWA significance level by chromosome
# location. Visual inspection of this plot allows for identification of SNPs 
# with relatively small p-values that are in regions with relatively large and 
# non-significant p-values, suggesting potentially false findings. Association 
# significance is established by p-values being above defined threshold such as 
# the Bonferroni corrected threshold or the less stringent suggestive 
# association threshold. 

# we use a function created by Reed et. al (2015) to generate this plot. 

# loading Manhattan Plot function 
source("2.Manhattan.plot.R")

# creating Manhattan Plot and saving it as a .pdf
pdf(paste(p.graphics, "ManhattanPlot.pdf", sep = ""))
par(mfrow = c(1,1))
manhattan.plot.s <- GWAS.manhattan(GWAS.comb)
dev.off()


# === 2) Q-Q plot ===============================================================

# Q-Q plots are used to visualize the relationship between the expected
# and observed distributions of SNP-level test statistics. We use the function 
# estlambda from the package GenABEL to generate our QQ plots. 

# creating QQ Plot and saving it as a .pdf
pdf(paste(p.graphics, "qq.plot.pdf", sep = ""))
lambda.plot.s <- estlambda(GWAS.output$t.value^2, plot = TRUE, 
                           method = "median")
dev.off()


#___ end _______________________________________________________________________


