################################################################################
# Scripts by Teresa Boquete, Marc Schmidt, Isolde van Riemsdijk and Lauren Walling

################################################################################
# Load packages (and common-functions script) silently
suppressPackageStartupMessages({
  library("data.table")
  library("ggplot2")
  library("tidyverse")
  library("gridExtra")
  library("hierfstat") # F statistics
  library("psych")     # Correlations among predictors
  library("vegan")     # RDA
  library("adegenet")  # Allele count 
  library("poppr")
  library("RColorBrewer")
  source("commonFunctions.R") 
})

################################################################################
# Define input files with full paths
sampleTable <- file.path("Design_with_coverage_info.txt")
infileSNP <- file.path("filteredSNPs.txt")

################################################################################
# Load input files
sampleTable <- read.table(sampleTable, sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
sampleTable$Pops1 <- paste0(sampleTable$GroupingA, "_", sampleTable$GroupingB)
sampleTable$Pops2 <- paste0(sampleTable$GroupingA)
myDataSNP <- f.load.filtered.SNPs(infileSNP)

# ##############################################################################
# GENETIC DATA ANALYSIS

# ##############################################################################
# Transform "./." into missing values (otherwise "." would be interpreted as an allele)
myDataSNP[(myDataSNP == "..")|(myDataSNP == "./.")] <- NA

# Check number of NAs in the file
sum(is.na(myDataSNP))

# Create genind object with adegenet
myDataSNP.t <- data.frame(t(myDataSNP))
colnames(myDataSNP.t) <- gsub("X", "", colnames(myDataSNP.t))
SNPgenind <- df2genind(myDataSNP.t, ploidy = 2, ind.names = rownames(sampleTable), pop = sampleTable$Pops2, sep = "/")

# Overall differentiation (Fst)
overallFst <- wc(SNPgenind)

# Pairwise genetic differentiation
pairwiseFst <- genet.dist(SNPgenind, method = "WC84")

# Write a table with the pairwise Fst output
write.table(as.matrix(pairwiseFst), "pairwiseFst_popestuary.txt", sep = "\t")

# Calculate pairwise Fst significance with a bootstrap, for the 95% CI
sign.paiw.fst <- boot.ppfst(dat = SNPgenind, nboot = 999, quant = c(0.025, 0.975),
                            diploid = TRUE)

# NOTE: Significance of the pairwise Fst: If your Confidence Interval includes ZERO,
# your Fst value is not significantly different from ZERO. If it does not include
# zero, then your Fst is different from zero.
#Check: https://www.researchgate.net/post/How_to_test_the_significance_of_FST_with_R_hierfstat

# Write lower and upper limits to files
write.table(as.matrix(sign.paiw.fst$ll), "LowerLimit_CI_pairwiseFst_popsite.txt", sep = "\t")
write.table(as.matrix(sign.paiw.fst$ul), "UpperLimit_CI_pairwiseFst_popsite.txt", sep = "\t")

# Fst tables with colors for publication
# Generate a palette with 300 shades
my_palette <- colorRampPalette(c("blue", "white"))(299)

#Defines the color breaks manually for an even color transition between the output values
col_breaks = c(seq(-0.04, 0.0, length = 100),         # for blue
               seq(0.0000000001, 0.04, length = 100), # for salmon
               seq(0.041, 0.09, length = 100))        # for reddish

#Example figure: Fst between three populations
Fst_3<-data.matrix(rev(read.table("pairwiseFst_popestuary.txt", header = T, row.names = 1)))

# The number of populations you want to display
dim_3 <- ncol(Fst_3)

# Make an image
image(1:dim_3, 1:dim_3, Fst_3, axes = FALSE, xlab = "", ylab = "",
      col = rev(my_palette), breaks = col_breaks)
axis(1, 1:dim_3, colnames(Fst_3), cex.axis = 0.7, las = 3)
axis(2, 1:dim_3, colnames(Fst_3), cex.axis = 0.7, las = 1)
text(expand.grid(1:dim_3, 1:dim_3), sprintf("%0.4f", Fst_3), cex = 0.54)

# Plot the legend
plot(seq_len(length(my_palette)), rep_len(1, length(my_palette)),
     col = my_palette, pch = 16, cex = 3, xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')

