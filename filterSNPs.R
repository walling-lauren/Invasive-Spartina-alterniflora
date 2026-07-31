
# define output file name
if (grepl("\\.vcf$", infileName)) {
  outfileName <- gsub("\\.vcf$", ".filtSNP", infileName)
} else {
  outfileName <- paste0(infileName, ".filtSNP")
}

## load packages silently
#install.packages("pegas")
suppressPackageStartupMessages({
 	library("vcfR")
	library("adegenet")
  library("pegas")
	library("plyr")
	library("RColorBrewer")
  library("dplyr")
  source(file.path(scriptDir, "commonFunctions.R"))
})

################################################################################################
### load data
sampleTab <- f.read.sampleTable(designTable, colNamesForGrouping) # see commonFunctions.R

# load the SNP data, extract GT and DP
temp <- f.read.vcf.file(infileName, sampleName)  # see commonFunctions.R
genoTypes <- temp$genoTypes
snpCoverage <- temp$snpCoverage

# match the samples
commonSamples <- sort(intersect(colnames(genoTypes), rownames(sampleTab)))
allSamples <- union(colnames(genoTypes), rownames(sampleTab))
samplesToRemove <- setdiff(allSamples, commonSamples)
if (length(samplesToRemove) > 0) {
  f.print.message("Removing", length(samplesToRemove), "samples!")
  cat(paste0(samplesToRemove, collapse = '\n'), '\n')
}
sampleTab <- sampleTab[commonSamples,]
genoTypes <- genoTypes[, commonSamples]
snpCoverage <- snpCoverage[, commonSamples]


################################################################################################
### filter
snpCoverage[(snpCoverage < minCov) | (snpCoverage > maxCov)] <- NA
tabForSummary <- data.frame(sample = rownames(sampleTab), group = sampleTab$group, stringsAsFactors = FALSE)
callsPerGroup <- f.summarize.columns(!is.na(snpCoverage), tabForSummary, sum)
maskToKeep <- rowSums(callsPerGroup >= minCountPerGroup) == ncol(callsPerGroup)

if (sum(maskToKeep) == 0) {
  f.print.message("No SNP passed the filter, no output.")
  quit("no", 0)
}

################################################################################################
### change the format of the genotype that it fits to Adegenet
genoTypes <- genoTypes[maskToKeep,]
genoTypes[is.na(genoTypes)] <- "./."
if (is.vector(genoTypes)) {
  genoTypes <- matrix(genoTypes, nrow = 1, ncol = length(genoTypes), dimnames = list(rownames(snpCoverage)[maskToKeep],names(genoTypes)))
}
write.table(genoTypes, file = outfileName, sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)

################################################################################################
### EXIT, the things below are just for testing

quit("no", 0)

################################################################################################

