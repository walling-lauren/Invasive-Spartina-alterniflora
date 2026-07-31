################################################################################
#                 dbRDA Spartina alterniflora, September 2021                  #
################################################################################


rm(list=ls())
## define input files:
designTable <- file.path("Design_File.csv")
allDistanceFiles <- c(file.path("epigeneticDistance_all.csv"),
                      file.path("filteredSNPs_Distances.csv"))
infileSNP <- file.path("filteredSNPs.txt")


## load packages (and common functions script) silently:
suppressPackageStartupMessages({
  library("vegan") # capscale
  library("ggplot2")
  library("cowplot")
  library("hierfstat") 
  library("adegenet")
  library("ape")
  source(file.path("commonFunctions.R"))
})

#####################################################################################
### load data
sampleTab <- f.read.sampleTable(designTable)

f.load.many.distance.files <- function(allDistanceFiles, imputeNA = TRUE) {
  require("ape")
  allDistances <- list()
  for (curFile in allDistanceFiles) {
    curFileName <- basename(curFile)
    curName <- gsub("\\.csv$", "", curFileName)
    temp <- read.csv(curFile, header = TRUE, stringsAsFactors = FALSE, row.names = 1)
    temp <- as.matrix(temp)
    diag(temp) <- 0
    temp[lower.tri(temp)] <- t(temp)[lower.tri(t(temp))]
    if (imputeNA) {
      if (mean(is.na(temp)) > 0.2) {
        cat("Skipping data set because more than 20% are NA.\n")
        next
      }
      if (sum(is.na(temp)) > 0) {
        prevNames <- dimnames(temp)
        temp <- additive(temp)
        dimnames(temp) <- prevNames
      }
      if (sum(temp < 0) > 0) {
        cat("Detected negative distances. Removing some samples, trying again\n")
        percNeg <- rowMeans(temp<0)
        toKeep <- setdiff(rownames(temp), names(percNeg)[percNeg>0.05]) # remove the samples with more than 5 % negative distances, then set neg distances to NA and remove them
        temp <- temp[toKeep, toKeep]
        temp[temp<0] <- NA
        temp <- f.remove.NA.from.distance.matrix(temp)
        if (nrow(temp) < 4) {
          cat("Skipping table with less than 4 samples after NA removal.\n")
          next
        }
      }
    } else {
      temp <- f.remove.NA.from.distance.matrix(temp)
      if (nrow(temp) < 4) {
        cat("Skipping table with less than 4 samples after NA removal.\n")
        next
      }
    }
    allDistances[[curName]] <- temp
  }
  return(allDistances)
}

allDistances <- f.load.many.distance.files(allDistanceFiles)

names(allDistances)
# [1] "epigeneticDistance_all_ALL" "geneticDistance_all" 

epiDist <- allDistances[[1]]
genDist <- allDistances[[2]]

myDataSNP <- f.load.filtered.SNPs(infileSNP)

# ##############################################################################
# GENETIC DATA PRE-PROCESSING

# ##############################################################################
# ### Transform ./. into missing values (otherwise "." would be interpreted as an allele)
myDataSNP[(myDataSNP == "..")|(myDataSNP == "./.")] <- NA

dim(myDataSNP)
# [1] 68317   211
dim(sampleTab)
# [1] 211  10
all(colnames(myDataSNP) %in% rownames(sampleTab))
# [1] TRUE
all(colnames(myDataSNP) == rownames(sampleTab))
# [1] TRUE


sum(is.na(myDataSNP)) / (dim(myDataSNP)[2] * dim(myDataSNP)[1])
# [[1] 0.09876123 => 8% of missing data

## IMPORTANT NOTE: to build a PCA with a genind object there can't be NAs in the dataset
# there are two options: (1) remove SNPs with any NA in the dataset; (2) keep
# NAs to build the genind object and then use the function scaleGen() to replalce NAs by
# mean allele frequencies (see below). The first option was problematic in the R. mangle dataset so
# I stick to the second option

################################################################################
# OPTION (1): Drop loci with missing values

#rows.with.na <- apply(myDataSNP, 1, function(x){any(is.na(x))})
#sum(rows.with.na)
# [1] 65854    by removing this I'd miss 96.4% of the SNPs (I'd only keep 2463 SNPs)


################################################################################
# OPTION (2): Build the genind object with missing values

myDataSNP.mt <- data.frame(t(myDataSNP))
colnames(myDataSNP.mt) <- gsub("X", "", colnames(myDataSNP.mt))
SNPgenind <- df2genind(myDataSNP.mt, ploidy = 2, ind.names = rownames(sampleTab), pop = sampleTab$population, sep = "/")
SNPgenind
#' /// GENIND OBJECT /////////
#' 
#' // 211 individuals; 68,317 loci; 139,617 alleles; size: 149 Mb
#' 
#' // Basic content
#' @tab:  211 x 139617 matrix of allele counts
#' @loc.n.all: number of alleles per locus (range: 2-4)
#' @loc.fac: locus factor for the 139617 columns of @tab
#' @all.names: list of allele names for each locus
#' @ploidy: ploidy of each individual  (range: 2-2)
#' @type:  codom
#' @call: df2genind(X = myDataSNP.mt, sep = "/", ind.names = rownames(sampleTab), 
#'                  pop = sampleTab$population, ploidy = 2)
#' 
#' // Optional content
#' @pop: population of each individual (group size range: 27-46)

# Obtain a matrix of scaled allele frequencies with genotypes (i.e. samples) in rows and alleles in columns
# by default it centers alleles frequencies to mean zero & scales alleles frequencies; the
# NA.method="mean" replaces NAs by the mean allele frequencies:
scaled.SNPgenind <- scaleGen(SNPgenind, NA.method="mean")
# Warning message:
#   In .local(x, ...) : Some scaling values are null.
# Corresponding alleles are removed.
dim(data.frame(scaled.SNPgenind))
# [1]   211 139455

# Do the PCA:
pca <- dudi.pca(scaled.SNPgenind, cent = FALSE, scale = FALSE, scannf = FALSE, nf = 3) # Retain 3 axes

# Create a barplot of the % of variance explained by each PC:
barplot(100*pca$eig/sum(pca$eig), main = "PCA eigenvalues",
        col = heat.colors(50), ylim = c(0,4))
title(ylab="Percent of genetic variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

# Get the % variance explained by PCs 1, 2, and 3:
(100*pca$eig/sum(pca$eig))[1:3]
# [1] 3.451675 2.567398 2.213612
sum((100*pca$eig/sum(pca$eig))[1:4])
# [1] [1] 10.40755
sum((100*pca$eig/sum(pca$eig))[1:10])  # the first 10 axes explain 19.8% of the total variance
# [1] 19.84603
sum((100*pca$eig/sum(pca$eig))[1:12])  # the first 12 axes explain 22.1% of the total variance
# [1] 22.14016

# Plot the first two PCs in a 2-dimensional space:
col <- funky(6)
s.class(pca$li, pop(SNPgenind),xax=1,yax=2, col=transp(col,.6),axesell=TRUE,
        addaxes = TRUE, cstar=0, cpoint=3, grid=FALSE)
title(ylab="PC1 (3.45%)", line = 2)
title(xlab="PC2 (2.57%)", line = 1)

col <- funky(6)
s.class(pca$li, pop(SNPgenind), xax=1, yax=3, col=transp(col,.6), axesell=TRUE,
        addaxes = TRUE, cstar=0, cpoint=3, grid=FALSE)
title(ylab="PC1 (3.45%)", line = 2)
title(xlab="PC3 (2.21%)", line = 1)

# Create a dataframe with the PCs scores that I'll use in the RDA:
pca <- dudi.pca(scaled.SNPgenind, cent = FALSE, scale = FALSE, scannf = FALSE, nf = 10) # Retain 10 axes, approx. 20% of variance
PCs_genetics <- round(pca$li, 4)


################################################################################
### dbRDA for genetic distances

### SITE MODEL

gen.dbRDA.site <- capscale(genDist ~ GroupingA, sampleTab)
anova(gen.dbRDA.site, permutations=9999) # is the model significant?
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# Model: capscale(formula = genDist ~ Site, data = sampleTab)
#           Df SumOfSqs      F Pr(>F)    
# Model      2 0.004760 7.9217  1e-04 ***
# Residual 208 0.062497                
# ---

RsquareAdj(gen.dbRDA.site)
# $r.squared
# [1] 0.07266453
# 
# $adj.r.squared
# [1] 0.06374784


### SITE MODEL

gen.dbRDA.site.hab <- capscale(genDist ~ GroupingA + GroupingB, sampleTab)
anova(gen.dbRDA.site.hab, permutations=9999) # is the model significant?
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# Model: capscale(formula = genDist ~ Site + Habitat, data = sampleTab)
#           Df SumOfSqs      F Pr(>F)    
# Model      3 0.005324 5.9319  1e-04 ***
# Residual 207 0.06193                
# ---

RsquareAdj(gen.dbRDA.site.hab)
# $r.squared
# [1] 0.08127257
# 
# $adj.r.squared
# [1] 0.06795768

gen.dbRDA.site.hab <- capscale(genDist ~ GroupingB, sampleTab)
anova(gen.dbRDA.site.hab, permutations=9999) # is the model significant?
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# Model: capscale(formula = genDist ~ Site + Habitat, data = sampleTab)
#           Df SumOfSqs      F Pr(>F)    
# Model      3 0.005324 5.9319  1e-04 ***
# Residual 207 0.06193                
# ---

RsquareAdj(gen.dbRDA.site.hab)


################################################################################
### dbRDA for epigenetic distances

epi.dbRDA <- capscale(epiDist ~ GroupingA + GroupingB + Condition(as.matrix(PCs_genetics)), sampleTab)
epi.dbRDA
anova(epi.dbRDA, permutations=9999) # is the model significant?
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# Model: capscale(formula = epiDist ~ Site + Habitat + Site:Habitat + Condition(as.matrix(PCs_genetics)), data = sampleTab)
#           Df SumOfSqs      F Pr(>F)    
# Model      5  0.04497 1.1481  1e-04 ***
# Residual 195  1.527449                
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RsquareAdj(epi.dbRDA)
# $r.squared
# [1] 0.0232559
# 
# $adj.r.squared
# [1] 0.003150587

anova(epi.dbRDA, by="terms", permutations=9999) # test the significance of each of the predictors
# Permutation test for capscale under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# Model: capscale(formula = epiDist ~ Site + Habitat + Site:Habitat + Condition(as.matrix(PCs_genetics)), data = sampleTab)
#               Df SumOfSqs      F Pr(>F)    
# Site           2  0.01777 1.1345 0.0006 ***
# Habitat        1  0.00861 1.0992 0.0166 *  
# Site:Habitat   2  0.01858 1.1863 0.0001 ***
# Residual     195  1.52744                 
# ---

### MODEL WITH GENETICS ONLY (PCs)

epi.dbRDA.onlyPCs <- capscale(epiDist ~ as.matrix(PCs_genetics), sampleTab)
epi.dbRDA.onlyPCs
anova(epi.dbRDA.onlyPCs, permutations=9999) # is the model significant?
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# Model: capscale(formula = epiDist ~ as.matrix(PCs_genetics), data = sampleTab)
#           Df SumOfSqs      F Pr(>F)    
# Model     10  0.36116 4.5938  1e-04 ***
# Residual 200  1.57241               
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RsquareAdj(epi.dbRDA.onlyPCs)
# $r.squared
# [1] 0.1867856
# 
# $adj.r.squared
# [1] 0.1461249


### MODEL WITh SITE ONLY

epi.dbRDA.site <- capscale(epiDist ~ GroupingA, sampleTab)
epi.dbRDA.site
anova(epi.dbRDA.site, permutations=9999) # is the model significant?
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# Model: capscale(formula = epiDist ~ Site, data = sampleTab)
#           Df SumOfSqs      F Pr(>F)    
# Model      2  0.09108 5.1413  1e-04 ***
# Residual 208  1.84249               
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RsquareAdj(epi.dbRDA.site)
#$r.squared
# [1] 0.04710673
# 
# $adj.r.squared
# [1] 0.03794429

### MODEL WITh HABITAT ONLY

epi.dbRDA.hab <- capscale(epiDist ~ GroupingB, sampleTab)
epi.dbRDA.hab
anova(epi.dbRDA.hab, permutations=9999) # is the model significant?
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# Model: capscale(formula = epiDist ~ Site, data = sampleTab)
#           Df SumOfSqs      F Pr(>F)    
# Model      2  0.09108 5.1413  1e-04 ***
# Residual 208  1.84249               
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RsquareAdj(epi.dbRDA.hab)
# $r.squared
# [1] 0.006136492
# 
# $adj.r.squared
# [1] 0.001381164

### MODEL WITh SITE + HABITAT

epi.dbRDA.site.hab <- capscale(epiDist ~ GroupingA + GroupingB, sampleTab)
epi.dbRDA.site.hab
anova(epi.dbRDA.site.hab, permutations=9999) # is the model significant?
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# Model: capscale(formula = epiDist ~ Site + Habitat, data = sampleTab)
#           Df SumOfSqs      F Pr(>F)    
# Model      3  0.10305 3.8843  1e-04 ***
# Residual 207  1.83052              
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RsquareAdj(epi.dbRDA.site.hab)
# $r.squared
# [1] 0.0532939
# 
# $adj.r.squared
# [1] 0.03957353

### MODEL WITh INTERACTION ONLY

epi.dbRDA.int <- capscale(epiDist ~ GroupingA:GroupingB, sampleTab)
epi.dbRDA.int
anova(epi.dbRDA.int, permutations=9999) # is the model significant?
# Permutation test for capscale under reduced model
# Permutation: free
# Number of permutations: 9999
# 
# Model: capscale(formula = epiDist ~ Site:Habitat, data = sampleTab)
#           Df SumOfSqs      F Pr(>F)    
# Model      5  0.13457 3.0669  1e-04 ***
# Residual 205  1.79900               
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
RsquareAdj(epi.dbRDA.int)
# $r.squared
# [1] 0.006136492
# 
# $adj.r.squared
# [1] 0.001381164

pvalues <- c(0.0365, .8132, .2166, .0007, .0529, .0003, .0802, .0002)
p.adjust(pvalues, method = "fdr")






