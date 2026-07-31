
# define output file name
if (grepl("\\.bed$", infileName)) {
  outfileName <- gsub("\\.bed$", ".filtMETH", infileName)
} else {
  outfileName <- paste0(infileName, ".filtMETH")
}

## load packages silently
suppressPackageStartupMessages({
  library("data.table")
  source(file.path(scriptDir, "commonFunctions.R"))
})

################################################################################################
### load data
sampleTab <- f.read.sampleTable(designTable, colNamesForGrouping) # see commonFunctions.R

# load the DNA methylation data
myData <- f.load.methylation.bed(infileName)
totalCols <- grep("_total$", colnames(myData), value = TRUE)
methCov <- myData[,totalCols]
colnames(methCov) <- gsub("_total$", "", colnames(methCov))

# match the samples
commonSamples <- sort(intersect(colnames(methCov), rownames(sampleTab)))
allSamples <- union(colnames(methCov), rownames(sampleTab))
samplesToRemove <- setdiff(allSamples, commonSamples)
if (length(samplesToRemove) > 0) {
  f.print.message("Removing", length(samplesToRemove), "samples!")
  cat(paste0(samplesToRemove, collapse = '\n'), '\n')
}
sampleTab <- sampleTab[commonSamples,]
methCov <- methCov[, commonSamples]
myData <- myData[, c("chr", "pos", "context", "samples_called", paste0(rep(commonSamples, each = 2), c("_methylated", "_total")))]

################################################################################################
### filter
methCov[(methCov < minCov) | (methCov > maxCov)] <- NA
tabForSummary <- data.frame(sample = rownames(sampleTab), group = sampleTab$group, stringsAsFactors = FALSE)
callsPerGroup <- f.summarize.columns(!is.na(methCov), tabForSummary, sum)
maskToKeep <- rowSums(callsPerGroup >= minCountPerGroup) == ncol(callsPerGroup)

if (sum(maskToKeep) == 0) {
  f.print.message("No cytosine passed the filter, no output.")
  quit("no", 0)
}

################################################################################################
### subset and store
myData <- myData[maskToKeep,]
if (is.vector(myData)) {
  myData <- matrix(myData, nrow = 1, ncol = length(myData), dimnames = list(rownames(methCov)[maskToKeep],names(myData)))
}
write.table(myData, file = outfileName, sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

################################################################################################
### EXIT, the things below are just for testing

quit("no", 0)

################################################################################################
### missing values
myMiss <- round(rowMeans(is.na(methCov))*100, 1)
summary(myMiss)

