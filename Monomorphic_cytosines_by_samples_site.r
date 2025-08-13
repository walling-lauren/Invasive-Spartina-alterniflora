#####################################################################################
# Epigenetic analyses for French Spartina (Oct.2023):                               #
# create data set without NA & check number of polymorphic loci                       #
# written by Teresa Boquete, some functions from Marc Schmidt
# modified by Lauren Walling et al.                                                        #
#####################################################################################

# Define input files
designTable <- ("Design.txt")
infileName <-  ("filteredMETH.bed")
colNamesForGrouping <- c("GroupingA", "GroupingB")

# Load packages silently
suppressPackageStartupMessages({
  library("data.table") # to read input files
  library("tidyverse")  # to use the pipe (%>%)
  library("reshape2")   # to stack columns on top of each other
  library("ggplot2")    # to do most plots
  library("vioplot")    # to do violin plots
  library("jtools")     # to summ() the models
  library("naniar")     # to plot missing values
  source("commonFunctions.R")  # Marc's custom functions 
})

#####################################################################################
# Load data
sampleTab <- f.read.sampleTable(designTable, colNamesForGrouping)
myData <- f.load.methylation.bed(infileName)

# Check if the dimensions of the data are correct
dim(myData)

# Transform  data for analysis
methCov <- myData[,grep("_total$", colnames(myData))]
colnames(methCov) <- gsub("_total$", "", colnames(methCov))

commonSamples <- intersect(rownames(sampleTab), colnames(methCov))
allSamples <- union(colnames(methCov), rownames(sampleTab))
samplesToRemove <- setdiff(allSamples, commonSamples)

if (length(samplesToRemove) > 0) {
  f.print.message("Removing", length(samplesToRemove), "samples!")
  cat(paste0(samplesToRemove, collapse = '\n'), '\n')
}

sampleTab <- sampleTab[commonSamples, ]
methCov <- methCov[ ,commonSamples]
all(rownames(sampleTab) == colnames(methCov))
all(rownames(sampleTab) == commonSamples)
myData <- myData[, c("chr", "pos", "context", "samples_called",
                     paste0(rep(commonSamples, each = 2), c("_methylated", "_total")))]

colnames(methCov) <- commonSamples
colnames(myData) <- c("chr", "pos", "context", "samples_called",
                      paste0(rep(commonSamples, each = 2), c("_methylated", "_total")))


minCov <- 10
maxCov <- quantile(methCov, 0.999, na.rm = TRUE)
methCov[(methCov < minCov) | (methCov > maxCov)] <- NA

# Plot coverage per sample
melt_data <- melt(methCov) 
p <- ggplot(data = melt_data, aes(x = variable, y = value)) +
  geom_boxplot(outlier.colour = "black", outlier.shape = 16,
               outlier.size = 2, notch = FALSE) +
  geom_hline(yintercept = 10, color = "blue", size = 1) +
  geom_hline(yintercept = 5751.909 , color = "blue", size = 1) +
  ylab("Contig coverage per sample") +
  xlab("Sample") +
  theme_classic()

################################################################################
# Plot and count the amount of missing values in this data set
pdf("missing_Cs.pdf")
  vis_miss(methCov, warn_large_data = FALSE)
dev.off()

NAs.per.sample <- data.frame(apply(methCov, 2, function(x) sum(is.na(x))))
colnames(NAs.per.sample) <- "MissingLoci"
NAs.per.sample$PercentMissing <- round((NAs.per.sample$MissingLoci / nrow(methCov)) * 100, 0)
NAs.per.sample$sample <- rownames(NAs.per.sample)

pdf("missing_Cs_perSample.pdf", width = 12, height = 8)
  NAs.per.sample %>% ggplot() +
    geom_bar(aes(x = sample, y = PercentMissing), stat = "identity",
             fill = "Blue", alpha = 0.7) +
    geom_text(aes(x = sample, y = PercentMissing, label = PercentMissing),
              vjust = 1.6, color = "white", size = 5, fontface = "bold") +
    theme_bw() + theme(legend.position="none") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(text = element_text(size = 17)) +
    scale_y_continuous(name = "Percent of missing loci") +
    scale_x_discrete(name = "")
dev.off()

################################################################################
# Extract loci without missing values
lociWithNA <- apply(methCov, 1, function(x) {any(is.na(x))})
sum(lociWithNA)
nrow(methCov) - sum(lociWithNA)
(nrow(methCov) - sum(lociWithNA))/nrow(methCov)*100

#####################################################################################
# Select loci without missing values => create complete data sets (no NAs)
# Vector of sample sites
sites <- unique(sampleTab$Block)

# First have to get NA for methCov and then use that same file for lociWithNA
for (i in 1:7){
  site <-sites[i]
  samples <- rownames(sampleTab[sampleTab$Block == site, ])
  submeth <- subset(methCov[ , samples])
  methcol <- paste(samples, "_methylated", sep = "")
  totalcol <- paste(samples, "_total", sep = "")
  twocol <- c(methcol, totalcol)
  submethContext <- myData[ , twocol]
  submethContext$context <- myData$context
  submethContext$chr <- myData$chr
  submethContext$pos <- myData$pos
  lociWithNA <- apply(submeth, 1, function(x) {any(is.na(x))})
  submethNoNA <- submeth[!lociWithNA, ]
  submethContextNoNA <- submethContext[!lociWithNA, ]
  filename <- paste(site, "_filtMETH.withoutNA.bed", sep = "")
  write.table(submethContextNoNA, filename, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
}

methCov <- methCov[!lociWithNA, ]
sum(is.na(submethContextNoNA))
dim(submethContextNoNA)

myData <- myData[!lociWithNA, ]
sum(is.na(myData))
dim(myData)

write.table(myData, "filtMETH.withoutNA.bed", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

#####################################################################################
#Load infileName for each sample site (Da,Db,Dc,Ha,Hb,Hc,Pn) and run following code separatley
infileName <- "Ha_filtMETH.withoutNA.bed"

# Load methylation data as percentages
mePerc <- f.load.methylation.bed(infileName, percentages = TRUE)
myData <- f.load.methylation.bed(infileName)

# Get and plot the number of cytosines per context 
CG_Cs <- sum(mePerc$context == "CG")
CHG_Cs <- sum(mePerc$context == "CHG")
CHH_Cs <- sum(mePerc$context == "CHH")

#input each sample name for both Total_Cs and Meth_Cs
CG_CsTotal <- sum(myData$Ha1_total, myData$Ha2_total, myData$Ha3_total, myData$Ha4_total, myData$Ha8_total, myData$Ha10_total)
CG_CsMeth <- sum(myData$Ha1_methylated, myData$Ha2_methylated, myData$Ha3_methylated, myData$Ha4_methylated, myData$Ha8_methylated, myData$Ha10_methylated)
Per_Cs <- (CG_CsMeth/ CG_CsTotal) * 100

c <- c(CG_Cs, CHG_Cs, CHH_Cs)
contexts <- c("CG", "CHG", "CHH")
df <- data.frame(context = contexts, numberC = c)

#Change name of pdf based on sample site
pdf("number_Cs_perContextHA.pdf", width = 8, height = 8)
  df %>% ggplot(aes(x = context, y = numberC)) +
    geom_bar(stat="identity", fill = "#1B9E77", alpha = 0.7) +
    geom_text(aes(label = numberC), vjust = 1.6, color = "white", size = 5,
              fontface = "bold") + theme_bw() + theme(legend.position="none") +
    theme(text = element_text(size = 17)) +
    scale_y_continuous(name = "Number of cytosines") +
    scale_x_discrete(name = "")
dev.off()

#####################################################################################
# Get the number of monomorphic loci for 0% methylation (C):
monomorphic_unmC <- which(rowSums(mePerc[ ,-c(1:3)] <= 5) >= round(ncol(mePerc[ ,-c(1:3)]) * 0.95, 0))
length(monomorphic_unmC)

# Percentage of Cs that are monomorphic for unmethylation
# and breakdown by context
length(monomorphic_unmC) / nrow(mePerc) * 100

mono_C <- mePerc[monomorphic_unmC, ]
sum(mono_C$context == "CG")
sum(mono_C$context == "CG") / CG_Cs * 100

sum(mono_C$context == "CHG")
sum(mono_C$context == "CHG") / CHG_Cs * 100

sum(mono_C$context == "CHH")
sum(mono_C$context == "CHH") / CHH_Cs * 100

################################################################################
# Number of monomorphic loci for 100% methylation (mC)
monomorphic_mC <- which(rowSums(mePerc[ ,-c(1:3)] >= 95) >= round(ncol(mePerc[ ,-c(1:3)]) * 0.95,0))
length(monomorphic_mC)
length(monomorphic_mC) / nrow(mePerc) * 100

# Percentage of Cs that are monomorphic for methylation
# and breakdown by context
length(monomorphic_mC) / nrow(mePerc) * 100
mono_mC <- mePerc[monomorphic_mC, ]

sum(mono_mC$context == "CG")
sum(mono_C$context == "CG") / CG_Cs * 100

sum(mono_mC$context == "CHG")
sum(mono_mC$context == "CHG") /CHG_Cs * 100

sum(mono_mC$context == "CHH")
sum(mono_mC$context == "CHH") / CHH_Cs * 100

(length(monomorphic_unmC) + length(monomorphic_mC)) / nrow(mePerc) * 100

################################################################################
# Select polymorphic loci
monomorphicLoci <- c(monomorphic_unmC, monomorphic_mC)
length(monomorphicLoci)

lociToKeep <- setdiff(rownames(mePerc), monomorphicLoci)
length(lociToKeep)

mePerc.polymorphic <- mePerc[lociToKeep, ]
myData <- myData[lociToKeep, ]
dim(myData)

write.table(myData, "filtMETH.polymorphicCs.bed", sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
