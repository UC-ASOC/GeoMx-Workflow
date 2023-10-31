##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Heewon Seo (Heewon.Seo@UCalgary.ca)
# Written on Oct 02, 2023
# Updated on Oct 30, 2023
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Directory settings
setwd("/home/rstudio/R")
source("userDefinedFunctions.R")

baseDir <- "/home/rstudio/analysis/"
outDir <- file.path(baseDir, "Results")
dir.create(outDir, showWarnings = FALSE)

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Load libraries
suppressPackageStartupMessages(library(NanoStringNCTools))
suppressPackageStartupMessages(library(GeoMxWorkflows))
suppressPackageStartupMessages(library(GeomxTools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsankey))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(DescTools))

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameters loading
paramFile <- file.path(baseDir, "config.yaml")

if (file.exists(paramFile)) { 
	config <- yaml::yaml.load_file(paramFile)

        slideOrder <- config$Slide$Name
        slideCols <- config$Slide$Color
        names(slideCols) <- slideOrder
        if (length(slideOrder) != config$Slide$Count) stop("Check YAML: Slide")

        patientOrder <- config$Patient$Name
        patientCols <- config$Patient$Color
        names(patientCols) <- patientOrder
        if (length(patientOrder) != config$Patient$Count) stop("Check YAML: Patient")

        regionOrder <- config$Region$Name
        regionCols <- config$Region$Color
        names(regionCols) <- regionOrder
        if (length(regionOrder) != config$Region$Count) stop("Check YAML: Region")

        segmentOrder <- config$Segment$Name
        segmentCols <- config$Segment$Color
        names(segmentCols) <- segmentOrder
        if (length(segmentOrder) != config$Segment$Count) stop("Check YAML: Segment")

        segmentQcParams <- list(
                minSegmentReads = config$QCparam$SegmentQC$minSegmentReads,
                percentTrimmed = config$QCparam$SegmentQC$percentTrimmed,
                percentStitched = config$QCparam$SegmentQC$percentStitched,
                percentAligned = config$QCparam$SegmentQC$percentAligned,
                percentSaturation = config$QCparam$SegmentQC$percentSaturation,
                minNegativeCount = config$QCparam$SegmentQC$minNegativeCount,
                maxNTCCount = config$QCparam$SegmentQC$maxNTCCount,
                minNuclei = config$QCparam$SegmentQC$minNuclei,
                minArea = config$QCparam$SegmentQC$minArea
        )

        segmentQC_colBy <- config$QCparam$Wrap$segmentQC_colBy
        segmentQC_rowBy <- config$QCparam$Wrap$segmentQC_rowBy
        segmentQC_trimmedThre <- config$SegmentQC$percentTrimmed
        segmentQC_stitchedThre <- config$SegmentQC$percentStitched
        segmentQC_alignedThre <- config$SegmentQC$percentAligned
        segmentQC_saturatedThre <- config$SegmentQC$percentSaturation
        segmentQC_nucleiThre <- config$SegmentQC$minNuclei
        segmentQC_areaThre <- config$SegmentQC$minArea

        probeQcParams <- list(
                minProbeRatio = config$QCparam$ProbeQC$minProbeRatio,
                percentFailGrubbs = config$QCparam$ProbeQC$percentFailGrubbs
        )

        loqCutoff <- config$QCparam$LOQ$loqCutoff
        loqMin <- config$QCparam$LOQ$loqMin

        geneDetectionRateThre <- config$QCparam$DetectionRate$geneDetectionRateThre
        geneDetectionRateBins <- unlist(config$QCparam$DerectionRateBins$geneDetectionRate)
        geneDetectionRateBinLabels <- config$QCparam$DerectionRateBins$geneDetectionRateLabel

        coefVariationThreshold <- config$VariableGenes$coefVariationThreshold
        howManyFeatures <- config$VariableGenes$numberOfFeatures

} else {
	stop("Provide config.yaml file with -v option")
}


##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 1 -- * --
# (1) Load DCC/PKC file into a GeoMx object, e.g., geomxSet.RDS
# (2) Export raw count matrix
# (3) Create a plot of the study design in sankey plot
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 01: Load DCC/PKC/Annotation files and export to an")
message("                 geomxObject.")
dccFiles  <- dir(file.path(baseDir, "DCC"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
pkcFile   <- dir(file.path(baseDir, "PKC"), pattern = ".pkc$", full.names = TRUE, recursive = TRUE)
annotFile <- dir(file.path(baseDir, "Annot"), pattern = ".xlsx$", full.names = TRUE, recursive = TRUE)
message(paste0("\n\t>> Files found in the working directory (", baseDir, "):"))
message(paste0("\t  - DCC       : ", length(dccFiles), " file(s) found"))
message(paste0("\t  - PKC       : ", length(pkcFile), " file(s) found"))
message(paste0("\t  - Annotation: ", length(annotFile), " file(s) found"))
if (length(dccFiles) < 1 || length(pkcFile) < 1 || length(annotFile) < 1) {
        stop("\nProvide DCC, PKC, and Annotation files")
} else {
        message("\n\t# Heads-up: Nanostring folks hard-coded in the GeomxTools package")
        message("\t         where AOISurfaceArea and AOINucleiCount column names should be")
        message("\t         area and nuclei (lower case), respectively.")
        message("\t         For more information, visit https://github.com/Nanostring-Biostats/GeomxTools/blob/master/R/NanoStringGeoMxSet-qc.R")
        message("\t")
}

# Load data
initSet <- readNanoStringGeoMxSet(dccFiles = dccFiles,
                           pkcFiles = pkcFile,
                           phenoDataFile = annotFile,
                           phenoDataSheet = "Annotation",
                           phenoDataDccColName = "FileName",
                           protocolDataColNames = c("ROI", "AOI"),
                           experimentDataColNames = c("Panel")
)
saveRDS(initSet, file.path(outDir, "01_initSet.GeoMx.RDS"))
message("\n\t>> GeoMx data exported to an GeoMxSet object:")
message(paste0("\t  - Probes : ", dim(initSet)[1]))
message(paste0("\t  - Samples: ", dim(initSet)[2]))

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 2 -- * --
# (1) Shift all zero values to one to transform in downstream analysis
gSet <- shiftCountsOne(initSet, useDALogic = TRUE)
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 02: Shift any expression counts with a value of zero (0) to")
message("                 one (1) to enable in downstream transformations.")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 3 -- * --
# (1) Quality check at the segment/sample level
# (2) Filter out segments/samples with low qual
module <- gsub(".pkc", "", annotation(gSet))

gSet <- setSegmentQCFlags(gSet, qcCutoffs = segmentQcParams)
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 03: Assess sequencing quality and adequate tissue sampling")
message("                 for every ROI/AOI segment.")

segmentQcResults <- protocolData(gSet)[["QCFlags"]]
flagColumns <- colnames(segmentQcResults)
qcSummary <- data.frame(Pass = colSums(!segmentQcResults[, flagColumns]), Warning = colSums(segmentQcResults[, flagColumns]))
segmentQcResults$QCStatus <- apply(segmentQcResults, 1L, function(x) { ifelse(sum(x) == 0L, "PASS", "WARNING") })
qcSummary["TOTAL FLAGS", ] <- c(sum(segmentQcResults[, "QCStatus"] == "PASS"), sum(segmentQcResults[, "QCStatus"] == "WARNING"))
qcSummary[,"TOTAL"] <- apply(qcSummary, 1, sum)

write.table(pasteListValues(segmentQcParams), file.path(outDir, "03_Segment_QC.txt"), row.names=F, col.names=F, quote=F, sep="\t")
colnames(qcSummary)[1] <- paste0("Filter\t", colnames(qcSummary)[1])
suppressWarnings( write.table(qcSummary, file.path(outDir, "03_Segment_QC.txt"), row.names=T, col.names=T, quote=F, sep="\t", append=T) )

qcMat <- sData(gSet)
qcMat$Segment <- factor(qcMat$Segment, levels = segmentOrder)
suppressWarnings({
        qc1 <- histQC(qcMat, "Trimmed (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_trimmedThre, cols = segmentCols)
        qc2 <- histQC(qcMat, "Stitched (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_stitchedThre, cols = segmentCols)
        qc3 <- histQC(qcMat, "Aligned (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_alignedThre, cols = segmentCols)
        qc4 <- histQC(qcMat, "Saturated (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_saturatedThre, cols = segmentCols)
        qc5 <- histQC(qcMat, "area", segmentQC_colBy, segmentQC_rowBy, segmentQC_areaThre, "log10", "AOI Area (log10)", cols = segmentCols)
        qc6 <- histQC(qcMat, "nuclei", segmentQC_colBy, segmentQC_rowBy, segmentQC_nucleiThre, "log10", "AOI nuclei count", cols = segmentCols)
})

pdf(file.path(outDir, "03_Segment_QC_before.pdf"))
print(qc1)
print(qc2)
print(qc3)
print(qc4)
print(qc5)
print(qc6)
invisible(dev.off())

gSet <- gSet[, segmentQcResults$QCStatus == "PASS"]
message("\n\t>> Removed flagged samples(segments):")
message(paste0("\t  - Probes : ", dim(gSet)[1]))
message(paste0("\t  - Samples: ", dim(gSet)[2]))

afterQcMat <- sData(gSet)
afterQcMat$Segment <- factor(afterQcMat$Segment, levels = segmentOrder)
suppressWarnings({
        qc1 <- histQC(afterQcMat, "Trimmed (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_trimmedThre, cols = segmentCols)
        qc2 <- histQC(afterQcMat, "Stitched (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_stitchedThre, cols = segmentCols)
        qc3 <- histQC(afterQcMat, "Aligned (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_alignedThre, cols = segmentCols)
        qc4 <- histQC(afterQcMat, "Saturated (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_saturatedThre, cols = segmentCols)
        qc5 <- histQC(afterQcMat, "area", segmentQC_colBy, segmentQC_rowBy, segmentQC_areaThre, "log10", "AOI Area (log10)", cols = segmentCols)
        qc6 <- histQC(afterQcMat, "nuclei", segmentQC_colBy, segmentQC_rowBy, segmentQC_nucleiThre, "log10", "AOI nuclei count", cols = segmentCols)
})

pdf(file.path(outDir, "03_Segment_QC_after.pdf"))
print(qc1)
print(qc2)
print(qc3)
print(qc4)
print(qc5)
print(qc6)
invisible(dev.off())
message("\n\t++ Created:")
message("\t  - Results/03_Segment_QC.txt")
message("\t  - Results/03_Segment_QC_before.pdf")
message("\t  - Results/03_Segment_QC_after.pdf")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 4 -- * --
# (1) Check negative probes/signal distribution
negCol <- "NegGeoMean"

negativeGeoMeans <- esBy(
        negativeControlSubset(gSet),
        GROUP = "Module", 
        FUN = function(x) { 
                assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
        }
)
protocolData(gSet)[[negCol]] <- negativeGeoMeans
pData(gSet)[, negCol] <- sData(gSet)[[negCol]]
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 04: Calculate negative count which is the geometric mean of")
message("                 the several unique negative probes in the GeoMx panel that")
message("                 do not target mRNA and establish the background count level")
message("                 per segment.")

backgrounMat <- sData(gSet)
backgrounMat <- backgrounMat[, c("NegGeoMean", segmentQC_colBy, segmentQC_rowBy)]
backgrounMat$Segment <- factor(backgrounMat$Segment, levels = segmentOrder)
qc <- histQC(backgrounMat, "NegGeoMean", segmentQC_colBy, segmentQC_rowBy, 2, "log10", "GeoMean(negative probes)", cols = segmentCols)
pdf(file.path(outDir, "04_Negative_probes.pdf"))
print(qc)
invisible(dev.off())
message("\n\t++ Created:")
message("\t  - Results/04_Negative_probes.pdf")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 5 -- * --
# (1) Quality check at the probe level
# (2) Filter out low qual probes (then, aggregate expression at the gene-level)
gSet <- setBioProbeQCFlags(gSet, qcCutoffs = probeQcParams, removeLocalOutliers = TRUE)
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 05: Remove low-performing probes. In short, this QC is an outlier")
message("                 removal process, whereby probes are either removed entirely from")
message("                 the study (global) or from specific segments (local).")
probeQcResults <- fData(gSet)[["QCFlags"]]

qcDf <- data.frame(
        Passed = sum(rowSums(probeQcResults[, -1]) == 0),
        Global = sum(probeQcResults$GlobalGrubbsOutlier),
        Local = sum(rowSums(probeQcResults[, -2:-1]) > 0 & !probeQcResults$GlobalGrubbsOutlier),
        TOTAL = nrow(probeQcResults)
)
write.table(pasteListValues(probeQcParams), file.path(outDir, "05_Probe_QC.txt"), row.names=F, col.names=F, quote=F, sep="\t")
suppressWarnings(write.table(qcDf, file.path(outDir, "05_Probe_QC.txt"), row.names=F, col.names=T, quote=F, sep="\t", append=T))
message("\n\t++ Created:")
message("\t  - Results/05_Probe_QC.txt")

probeQcPassed <- subset(
        gSet, 
        fData(gSet)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
        fData(gSet)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE
)
gSet <- probeQcPassed

message("\n\t>> Removed flagged probes:")
message(paste0("\t  - Probes : ", dim(gSet)[1]))
message(paste0("\t  - Samples: ", dim(gSet)[2]))

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 6 -- * --
# (1) Aggregate multi-probe and generate a gene expression profile
newSet <- aggregateCounts(gSet)

#Exclude negative control from the gene list
newSet <- subset(newSet, fData(newSet)$TargetName != "NegProbe-WTX")
saveRDS(as.data.frame(pData(newSet)), file.path(outDir, "00_finalStats.RDS"))

message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 06: Generate a gene-level count matrix where the count for")
message("                 any gene with multiple probes per segment is calculated as")
message("                 the geometric mean of those probes.")
message("\n\t>> Collapse to targets/genes:")
message(paste0("\t  - Genes  : ", dim(newSet)[1]))
message(paste0("\t  - Samples: ", dim(newSet)[2]))

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 7 -- * --
# (1) Calculate LOQ
# (2) Filter out low signal segments/genes compared to the background
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 07-1: Check the distribution of the number of nuclei estimated")
message("                 and area of the AOIs, and their correlations.")

lowLvqcDf <- data.frame(
        Slide = pData(newSet)$Slide,
        Sample = pData(newSet)$Library,
        Segment = pData(newSet)$Segment,
        Area = pData(newSet)$area,
        Nuclei = pData(newSet)$nuclei
)
lowLvqcDf$Slide <- factor(lowLvqcDf$Slide, levels=slideOrder)
lowLvqcDf$Segment <- factor(lowLvqcDf$Segment, levels=segmentOrder)

dodge <- position_dodge(width = 0.5)
suppressWarnings({
        violin1 <- ggplot(data = lowLvqcDf, aes(x = Segment, y = Area, fill = Segment)) +
                geom_violin(position = dodge, size = 0) +
                geom_boxplot(width = 0.1, position = dodge, fill="white") +
                scale_fill_manual(values = segmentCols) +
                facet_wrap( ~ Slide) +
                labs(
                        title = "",
                        x = "", 
                        y = "Area"
                ) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                        legend.position = "none", 
                        text = element_text(size = 12)
                ) + 
                scale_y_continuous(trans = "log10")
})

suppressWarnings({
        violin2 <- ggplot(data = lowLvqcDf, aes(x = Segment, y = Nuclei, fill = Segment)) +
                geom_violin(position = dodge, size = 0) +
                geom_boxplot(width = 0.1, position = dodge, fill="white") +
                scale_fill_manual(values = segmentCols) +
                facet_wrap( ~ Slide) +
                labs(
                        title = "",
                        x = "", 
                        y = "#Nuclei"
                ) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                        legend.position = "none", 
                        text = element_text(size = 12)
                ) + 
                scale_y_continuous(trans = "log10")
})

suppressWarnings({
        scatterPlot <- ggplot(data = lowLvqcDf, aes(x = Area, y = Nuclei, color = Segment, label = Sample, label2 = Slide)) +
                geom_point() +
                scale_color_manual(values = segmentCols) +
                labs(
                        title = "",
                        x = "Area", 
                        y = "#Nuclei"
                ) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                        text = element_text(size = 12)
                ) + 
                scale_x_continuous(trans = "log10") + 
                scale_y_continuous(trans = "log10")
})

pdf(file.path(outDir, "07_1_Low_level_QC.pdf"))
print(violin1)
print(violin2)
print(scatterPlot)
invisible(dev.off())
message("\n\t++ Created:")
message("\t  - Results/07_1_Low_level_QC.txt")

loqDf <- data.frame(row.names = colnames(newSet))
varNames <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
if (all(varNames[1:2] %in% colnames(pData(newSet)))) {
        loqDf[, module] <- pData(newSet)[, varNames[1]] * (pData(newSet)[, varNames[2]]^loqCutoff)
}
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 07-2: Determine the limit of quantification (LOQ) per segment")
message("                 where the LOQ is calculated based on the distribution of")
message("                 negative control probes and is intended to approximate")
message("                 the quantifiable limit of gene expression per segment.")

statDf <- data.frame(
        Slide = pData(newSet)$Slide,
        Sample = pData(newSet)$Sample,
        Segment = pData(newSet)$Segment,
        LOQ = loqDf[,1]
)
statShort <- dcast(statDf, Sample~Segment,  fun.aggregate = length, value.var = "LOQ")
perSample <- statShort[,c(2:ncol(statShort))]
rownames(perSample) <- as.character(statShort[,1])
statShort <- data.frame(
        Sample = rownames(perSample),
        LOQmean = as.matrix(apply(perSample, 1, mean, na.rm=TRUE))
)
statShort <- statShort[order(statShort$LOQmean),]
statDf$Sample <- factor(statDf$Sample, levels=statShort$Sample)
statDf$Segment <- factor(statDf$Segment, levels=segmentOrder)

dodge <- position_dodge(width = 0.5)
suppressWarnings({
        violin3 <- ggplot(data = statDf, aes(x = Segment, y = log10(LOQ), fill = Segment)) +
                geom_violin(position = dodge, size = 0) +
                geom_boxplot(width = 0.1, position = dodge, fill="white") +
                scale_fill_manual(values = segmentCols) +
                labs(
                        title = "",
                        subtitle = "",
                        x = "Segment", 
                        y = "LOQ, log10"
                ) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                        legend.position = "none", 
                        text = element_text(size = 12)
                ) + 
                geom_hline(aes(yintercept = log10(loqMin)), lty=2, col="grey50")
})

suppressWarnings({
        violin4 <- ggplot(data = statDf, aes(x = Segment, y = log10(LOQ), fill = Segment)) +
                geom_violin(position = dodge, size = 0) +
                geom_boxplot(width = 0.1, position = dodge, fill="white") +
                scale_fill_manual(values = segmentCols) +
                facet_grid(~Slide) + 
                labs(
                        title = "",
                        subtitle = "",
                        x = "Segment", 
                        y = "LOQ, log10"
                ) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0.5),
                        legend.position = "none", 
                        text = element_text(size = 12)
                ) + 
                geom_hline(aes(yintercept = log10(loqMin)), lty=2, col="grey50")
})
 
pdf(file.path(outDir, "07_2_LOQ_Distribution.pdf"))
print(violin3)
print(violin4)
invisible(dev.off())
message("\t  - Results/07_2_LOQ_Distribution.pdf")

suppressWarnings({
        scatterPlot2 <- ggplot(data=statDf, aes(x = Sample, y = log10(LOQ), group = Segment)) +
                geom_line(aes(color=Segment), lwd=0.5) +
                geom_point(aes(color=Segment)) +
                scale_color_manual(values = segmentCols) +
                labs(
                        title = "",
                        x = "", 
                        y = "LOQ, log10"
                ) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0.5),
                        text = element_text(size = 12)
                ) + 
                geom_hline(aes(yintercept = log10(loqMin)), lty=2, col="grey50")
})

pdf(file.path(outDir, "07_2_LOQ_perSample.pdf"), width=15)
print(scatterPlot2)
invisible(dev.off())
message("\t  - Results/07_2_LOQ_perSample.pdf")

loqDf[loqDf < loqMin] <- loqMin
pData(newSet)$LOQ <- loqDf

loqMat <- t(esApply(newSet, MARGIN = 1, FUN = function(x) { x > LOQ[, module] }) )
loqMat <- loqMat[fData(newSet)$TargetName, ] # Ordering

# Segment gene detection
pData(newSet)$GenesDetected <- colSums(loqMat, na.rm = TRUE)
pData(newSet)$GeneDetectionRate <- pData(newSet)$GenesDetected / nrow(newSet)
pData(newSet)$DetectionThreshold <- cut(pData(newSet)$GeneDetectionRate, breaks = geneDetectionRateBins, labels = geneDetectionRateBinLabels)
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 07-3: Filter out segments with exceptionally low signal that")
message("                 have a small fraction of panel genes detected above the LOQ")
message("                 relative to the other segments in the study.")

write.table(paste0("# Gene detection rate threshold > ", geneDetectionRateThre), file.path(outDir, "07_3_Segment_gene_detection.txt"), row.names=F, col.names=F, quote=F, sep="\t")
suppressWarnings(write.table(t(table(pData(newSet)$DetectionThreshold)), file.path(outDir, "07_3_Segment_gene_detection.txt"), row.names=F, col.names=T, quote=F, sep="\t", append=T))
message("\n\t++ Created:")
message("\t  - Results/07_3_Segment_gene_detection.txt")

rateMat <- pData(newSet)
rateMat$Segment <- factor(rateMat$Segment, levels = segmentOrder)
suppressWarnings({
        barplot <- ggplot(rateMat, aes(x = DetectionThreshold)) +
                geom_bar(aes(fill = Segment)) +
                geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
                scale_fill_manual(values = segmentCols) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0),
                        text = element_text(size = 12)
                ) + 
                scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
                labs(x = "Gene Detection Rate",
                        y = "Number of Segments/Samples",
                        fill = "Segment"
                ) + 
                facet_grid(as.formula(paste("~", segmentQC_rowBy)))
})

pdf(file.path(outDir, "07_3_DetectedGenes_summary.pdf"))
print(barplot)
invisible(dev.off())
message("\t  - Results/07_3_DetectedGenes_summary.txt")

statDf <- data.frame(
        Sample = protocolData(newSet)$AOI, # pData(newSet)$Sample,
        Segment = pData(newSet)$Segment,
        GenesDetected = colSums(loqMat, na.rm = TRUE),
        Nuclei = pData(newSet)$nuclei
)
statDf$Segment <- factor(statDf$Segment, segmentOrder)

pdf(file.path(outDir, "07_3_DetectedGenes.pdf"))
for (segment in segmentOrder) {
        tmpDf <- statDf[which(statDf$Segment == segment),]
        tmpDf <- tmpDf[order(tmpDf$GenesDetected, decreasing = F),]
        tmpDf$Sample <- factor(tmpDf$Sample, levels=tmpDf$Sample)
        suppressWarnings({
                barplot <- ggplot(tmpDf, aes(x = Sample, y = GenesDetected, fill = Segment)) +
                        geom_bar(stat = "identity") +
                        scale_fill_manual(values = segmentCols[segment]) +
                        theme_minimal() +
                        labs(
                                title = "",
                                x = "Sample"
                        ) +
                        coord_flip() +
                        theme_bw() +
                        theme(
                                axis.line = element_line(colour = "black"),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.border = element_blank(),
                                panel.background = element_blank(),
                                axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0),
                                legend.position = "none"
                        )
        })
        print(barplot)
}
invisible(dev.off())
message("\t  - Results/07_3_DetectedGenes.txt")

newSet <- newSet[, pData(newSet)$GeneDetectionRate >= geneDetectionRateThre]
message("\n\t>> After filtering out low signal sample/segment based on LOQ:")
message(paste0("\t  - Genes  : ", dim(newSet)[1]))
message(paste0("\t  - Samples: ", dim(newSet)[2]))

# Gene detection rate
loqMat <- loqMat[, colnames(newSet)]
fData(newSet)$DetectedSegments <- rowSums(loqMat, na.rm = TRUE)
fData(newSet)$DetectionRate <- fData(newSet)$DetectedSegments / nrow(pData(newSet))
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 07-4: Determine the detection rate for genes across the study")
message("                 where individual genes are detected to varying degrees in")
message("                 the segments. In other words, we will calculate the total number")
message("                 of genes detected in different percentages of segments to filter")
message("                 out low detected genes.")

geneDetectionRateDf <- data.frame(
        Gene = rownames(newSet),
        Number = fData(newSet)$DetectedSegments,
        DetectionRate = percent(fData(newSet)$DetectionRate)
)
geneDetectionRateDf <- geneDetectionRateDf[order(geneDetectionRateDf$Number, geneDetectionRateDf$DetectionRate, geneDetectionRateDf$Gene),]
write.table(geneDetectionRateDf, file.path(outDir, "07_4_Gene_detection_rate.txt"), row.names=F, col.names=T, quote=F, sep="\t")

negativeProbefData <- subset(fData(newSet), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
finalSet <- newSet[fData(newSet)$DetectionRate >= geneDetectionRateThre | fData(newSet)$TargetName %in% neg_probes, ]
saveRDS(finalSet, file.path(outDir, "07_finalSet.GeoMx.RDS"))
message("\n\t>> After filtering out low signal gene:")
message(paste0("\t  - Genes  : ", dim(finalSet)[1]))
message(paste0("\t  - Samples: ", dim(finalSet)[2]))

countDf <- pData(finalSet) %>% make_long(Region, Segment, Patient)
suppressWarnings({
        studyDesign <- ggplot(countDf, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
                geom_sankey(flow.alpha = .6, node.color = "gray30") +
                geom_sankey_label(size = 3, color = "black", fill = "white") +
                scale_fill_viridis_d(option = "A", alpha = 0.95) +
                theme_sankey(base_size = 18) +
                labs(
                        title = "Study Design",
                        x = NULL, 
                        y = NULL
                ) +
                theme_bw() +
                theme(
                        axis.line = element_blank(),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                        axis.ticks.x = element_blank(),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank(),
                        legend.position = "none", 
                        text = element_text(size = 12)
                )
})

pdf(file.path(outDir, "07_4_Study_design.pdf"))
suppressWarnings(print(studyDesign))
invisible(dev.off())
message("\n\t++ Created:")
message("\t  - Results/07_4_Study_design.pdf")
message("\t  - Results/07_4_Gene_detection_rate.txt")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 8 -- * --
# (1) Normalization
finalSet <- normalize(finalSet, norm_method = "quant", desiredQuantile = .75, toElt = "q_norm")
# finalSet <- normalize(finalSet, norm_method = "neg", fromElt = "exprs", toElt = "neg_norm") # Background normalization for WTA without custom spike-in
annot <- pData(finalSet)
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 08: Upper-quartile (Q3) normalization method estimates")
message("                 a normalization factor per segment to bring the segment")
message("                 data distributions together.")

segmentIdx <- c()
for (segment in segmentOrder) segmentIdx <- c(segmentIdx, which(annot$Segment == segment))

rawMat <- assayData(finalSet)$exprs
colnames(rawMat) <- annot$Library
rawDf <- melt(rawMat)
colnames(rawDf) <- c("Gene", "Sample", "Expression")
rawDf$Segment <- rep(annot$Segment, each=nrow(rawMat))
rawDf$Segment <- factor(rawDf$Segment, levels = segmentOrder)

normMat <- assayData(finalSet)$q_norm
colnames(normMat) <- annot$Library
normDf <- melt(normMat)
colnames(normDf) <- c("Gene", "Sample", "Expression")
normDf$Segment <- rep(annot$Segment, each=nrow(rawMat))
normDf$Segment <- factor(normDf$Segment, levels = segmentOrder)

pdf(file.path(outDir, "08_Before_norm_perSegment.pdf"))
suppressWarnings(boxplot(log10(Expression) ~ Segment, data = rawDf, col = segmentCols, pch = 20, ylab = "Expression, log10", xlab = "Raw count"))
invisible(dev.off())

rawMat <- rawMat[, segmentIdx]
pdf(file.path(outDir, "08_Before_norm_perSample.pdf"), width=12)
suppressWarnings(boxplot(log10(rawMat), col = segmentCols[annot$Segment[segmentIdx]], pch = 20, names = NA, xlab = "Sample", ylab = "Raw count, log10"))
invisible(dev.off())

pdf(file.path(outDir, "08_After_norm_perSegment.pdf"))
suppressWarnings(boxplot(log10(Expression) ~ Segment, data = normDf, col = segmentCols, pch = 20, ylab = "", xlab = "Upper-quartile norm"))
invisible(dev.off())

normMat <- normMat[, segmentIdx]
pdf(file.path(outDir, "08_After_norm_perSample.pdf"), width=12)
suppressWarnings(boxplot(log10(normMat), col = segmentCols[annot$Segment[segmentIdx]], pch = 20, names = NA, xlab = "Sample", ylab = "Upper-quartile norm, log10"))
invisible(dev.off())

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 9 -- * --
# (1) PCA
normMat <- assayData(finalSet)$q_norm
pcaRes <- prcomp(t(normMat), scale = TRUE)
pcaDf <- data.frame(
        Sample = annot$Library,
        X = pcaRes$x[, 1],
        Y = pcaRes$x[, 2],
        Slide = annot$Slide,
        Segment = annot$Segment,
        Region = annot$Region
)
pcaDf$Slide <- factor(pcaDf$Slide, levels = slideOrder)
pcaDf$Segment <- factor(pcaDf$Segment, levels = segmentOrder)
pcaDf$Region <- factor(pcaDf$Region)

suppressWarnings({
        pcaPlot1 <- ggplot(data = pcaDf, aes(x = X, y = Y, color = Slide, label = Sample)) +
                geom_point(size = 2, shape = 20) +
                scale_color_manual(values = slideCols) +
                labs(
                        x = paste0("PC1 (", round(summary(pcaRes)$importance[2, c(1)] * 100, 1), "%)"),
                        y = paste0("PC2 (", round(summary(pcaRes)$importance[2, c(2)] * 100, 1), "%)")
                ) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                        text = element_text(size = 12)
                )
})
pdf(file.path(outDir, "09_PCA_perSlide.pdf"))
print(pcaPlot1)
invisible(dev.off())

suppressWarnings({
        pcaPlot2 <- ggplot(data = pcaDf, aes(x = X, y = Y, color = Segment, label = Sample)) +
                geom_point(size = 2, shape = 20) +
                scale_color_manual(values = segmentCols) +
                labs(
                        x = paste0("PC1 (", round(summary(pcaRes)$importance[2, c(1)] * 100, 1), "%)"),
                        y = paste0("PC2 (", round(summary(pcaRes)$importance[2, c(2)] * 100, 1), "%)")
                ) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                        text = element_text(size = 12)
                )
})
pdf(file.path(outDir, "09_PCA_perSegment.pdf"))
print(pcaPlot2)
invisible(dev.off())

suppressWarnings({
        pcaPlot <- ggplot(data = pcaDf, aes(x = X, y = Y, color = Region, label = Sample)) +
                geom_point(size = 2, shape = 20) +
                scale_color_manual(values = regionCols) +
                labs(
                        x = paste0("PC1 (", round(summary(pcaRes)$importance[2, c(1)] * 100, 1), "%)"),
                        y = paste0("PC2 (", round(summary(pcaRes)$importance[2, c(2)] * 100, 1), "%)")
                ) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank(),
                        axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5),
                        text = element_text(size = 12)
                )
})
pdf(file.path(outDir, "09_PCA_perRegion.pdf"))
print(pcaPlot)
invisible(dev.off())

message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
file.copy(getHistoryFile(), file.path(outDir, "00_Rhistory.R"))
message("Fin")

q("no")
