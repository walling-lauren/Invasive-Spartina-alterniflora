
designTable <- as.character("Designtable.txt")
annotationFile <- as.character("mergedAnnot.csv")
infileName <- as.character("filteredMETH.USOnly.bed")

## load packages silently
suppressPackageStartupMessages({
  library("data.table") # file reading
  library("vioplot") # plotting
  source(file.path("R_scripts/commonFunctions.R"))
})

## create the result directory if it does not exist
#if (!dir.exists(rDir)) { dir.create(rDir, FALSE) }

################################################################################################
### load data
#sampleTab <- f.read.sampleTable("forCommonFilter.design.US.txt") # see commonFunctions.R
sampleTab <- read.table(designTable, sep = '\t', header = TRUE, stringsAsFactors = FALSE, row.names = 1)
mePerc <- f.load.methylation.bed(infileName, percentages = TRUE) # see commonFunctions.R
infoColumns <- c("chr", "pos", "context")
#allSamples <- setdiff(colnames(mePerc), infoColumns)
#sampleTab <- sampleTab[allSamples,]

###Making sure the samples are the same between mePerc and sampleTab
#### need to make it so that the first three columns still stay for mePerc
commonSamples <- sort(intersect(colnames(mePerc), rownames(sampleTab)))
allSamples <- union(colnames(mePerc), rownames(sampleTab))
allSamples<- allSamples[!allSamples %in% c("chr", "pos", "context")]
samplesToRemove <- setdiff(allSamples, commonSamples)
if (length(samplesToRemove) > 0) {
  f.print.message("Removing", length(samplesToRemove), "samples!")
  cat(paste0(samplesToRemove, collapse = '\n'), '\n')
}
sampleTab <- sampleTab[commonSamples,]
mePerc <- mePerc[, c(infoColumns, commonSamples)]

################################################################################################
## average within groups (select on what to average)
#sampleTab$TS__ <- paste0(sampleTab$GroupingA, '_', sampleTab$GroupingB)
#aveData <- f.summarize.columns(mePerc, data.frame(sample = rownames(sampleTab), group = sampleTab$TS__, stringsAsFactors = FALSE), function(x) mean(x, na.rm = TRUE))
aveData <- f.summarize.columns(mePerc, data.frame(sample = rownames(sampleTab), group = sampleTab$GroupingA, stringsAsFactors = FALSE), function(x) mean(x, na.rm = TRUE))
aveDataInfo <- mePerc[,infoColumns]
rownames(aveDataInfo) <- paste0("chr", aveDataInfo$chr, "_", aveDataInfo$pos)
rownames(aveData) <- rownames(aveDataInfo)

################################################################################################
## choose a group order for the plot and set the colors
forPlotOrder <- sort(colnames(aveData))
#aveData <- aveData[, match(forPlotOrder, colnames(aveData))]
temp <- unique(sampleTab[,c("TS__", "popcol")])
#plotColors <- temp$color; names(plotColors) <- temp$TS__
plotColors <- temp$popcol; names(plotColors) <- temp$TS__
################################################################################################
## Draw violin plots
allContextsInData <- unique(aveDataInfo$context)
checkSumNormal <- !grepl("^m", allContextsInData)
checkSumMixed <- grepl("^m", allContextsInData)
allContexts <- c()
if (sum(checkSumNormal) > 0) { allContexts <- c("CG", "CHG", "CHH", allContexts) }
if (sum(checkSumMixed) > 0) { allContexts <- c(allContexts, "mCG", "mCHG", "mCHH") }
numPlotRows <- 1 # just one, there would also be the option that you do one row of plots per group and then within each figure separate groups
numPlotCols <- length(allContexts)
allMeans <- matrix(NA, nrow = length(forPlotOrder), ncol = numPlotCols, dimnames = list(forPlotOrder, allContexts))
aveData <- aveData[rownames(aveDataInfo),] # that would not be necessary, but just in case you modify things further up
pdf(file.path(rDir, "methylationLevelsViolinPlot.pdf"), height = 5*numPlotRows, width = 2+length(forPlotOrder)*numPlotCols)
par(oma = c(10, 10, 0, 0), mar = c(0, 0, 0, 0))
layout(matrix(1:(numPlotRows*numPlotCols), nrow = numPlotRows, byrow = TRUE))
for (ctxt in allContexts) {
  if (ctxt == "ALL") {
    subData <- aveData[aveDataInfo$context %in% c("CG", "CHG", "CHH"),]
  } else if (ctxt == "mALL") {
    subData <- aveData[aveDataInfo$context %in% c("mCG", "mCHG", "mCHH"),]
  } else {
    subData <- aveData[aveDataInfo$context == ctxt,]
  }
  plot(NA, main = "", bty = "n", xaxs = "r", yaxs = "r", xlab = "", ylab = "", las = 1, cex = 0.4, tck = 0.01, xlim = c(0.5, length(forPlotOrder)+0.5), ylim = c(0, 100), xaxt = "n", yaxt = "n")
  curPos <- 1
  for (curGroup in forPlotOrder) {
    toPlot <- subData[,curGroup]
    toPlot <- toPlot[!is.na(toPlot)]
    curCol <- plotColors[curGroup]
    if (sum(toPlot > 0) > 4) {
      vioplot(toPlot, names = c(curGroup), col = curCol, ylim = c(0,100), drawRect = TRUE, add = TRUE, at = curPos)
    }
    curMean <- mean(toPlot)
    lines(c(curPos-0.3,curPos+0.3), c(curMean, curMean), col = "black", lwd = 4, lty = 1)
    curPos <- curPos + 1
    allMeans[curGroup,ctxt] <- curMean # add the mean to the collection
  }
  if (ctxt == "all") { axis(2, at = seq(0, 100, by = 20), labels = seq(0, 100, by = 20), outer = TRUE, las = 1, line=2, lwd=2, cex.axis=3) }
  axis(1, at = 1:length(forPlotOrder), labels = forPlotOrder, outer = TRUE, las = 2, line=2, lwd=2, cex.axis=3)
}
invisible(dev.off())
write.csv(round(allMeans, 3), file.path(rDir, "methylationLevelsViolinPlot_means.csv"))

################################################################################################
## Get the average methylation level per group, context, feature
# it's ordered in a way that fits to the violinplot above
# NOTE: you can't assign a single feature to the contigs, the gene and repeatmasker annotation are not mutually exclusive
# that's why the code is a bit more lengthy 
# we'll only do it for gene, transposon, repeat, and nothing (unannotated), if you have a lot of data, you can also do it for the detailed annotation
allFeatures <- c("gene", "transposon", "repeat", "nothing")
forMask <- paste0("chr", aveDataInfo$chr)
listForPlot <- list()
for (ctxt in allContexts) {
  if (ctxt == "all") {
    contextMask <- rep(TRUE, nrow(aveDataInfo))
  } else {
    contextMask <- aveDataInfo$context == ctxt
  }
  featureMeans <- matrix(NA, nrow = length(allFeatures), ncol = length(forPlotOrder), dimnames = list(allFeatures, forPlotOrder))
  for (feature in allFeatures) {
    mergedAnno <- f.load.merged.annotation(annotationFile, feature)
    annoMask <- forMask %in% rownames(mergedAnno)
    subData <- aveData[annoMask & contextMask,]
    if (sum(annoMask & contextMask) > 1) {
      featureMeans[feature, colnames(aveData)] <- colMeans(subData, na.rm = TRUE)
    } else {
      featureMeans[feature, colnames(aveData)] <- NA
    }
  }
  listForPlot[[ctxt]] <- featureMeans
}

################################################################################################
## Do the plot, you could have also done it above, but it's sometimes easier to copy-paste later on if it's split :)
imageColors <- f.blackblueyellowredpinkNICE(51) # I need to send you a legend for this
pdf(file.path(rDir, "methylationLevelsPerFeature.pdf"), height = 5*numPlotRows, width = 2+length(forPlotOrder)*numPlotCols)
layout(matrix(1:(numPlotRows*numPlotCols), nrow = numPlotRows, byrow = TRUE))
for (ctxt in allContexts) {
  temp <- listForPlot[[ctxt]]
  f.image.without.text(forPlotOrder, allFeatures, t(temp), xLabel = "", yLabel = "", mainLabel = "", useLog = FALSE, col = imageColors, zlim = c(0, 100))
  write.csv(round(temp, 3), file.path(rDir, paste0("methylationLevelsPerFeature_", ctxt, ".csv")))
}
invisible(dev.off())

################################################################################################



