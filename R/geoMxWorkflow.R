##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Heewon Seo (Heewon.Seo@UCalgary.ca)
# Written on Oct 02, 2023
# Updated on Oct 23, 2023
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Settings
setwd("/home/rstudio/R")
source("userDefinedFunctions.R")

baseDir <- "/home/rstudio/analysis/"
outDir <- file.path(baseDir, "Results")
dir.create(outDir, showWarnings = FALSE)

paramFile <- file.path(baseDir, "Settings", "parameterSettings.R")
if (file.exists(paramFile)) { 
	source(paramFile)
} else {
	stop("Provide parameterSettings.R file with -v option")
}

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Load libraries
suppressPackageStartupMessages(library(NanoStringNCTools))
suppressPackageStartupMessages(library(GeoMxWorkflows))
suppressPackageStartupMessages(library(GeomxTools))

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggsankey))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(RColorBrewer))

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
gSet <- readNanoStringGeoMxSet(dccFiles = dccFiles,
                           pkcFiles = pkcFile,
                           phenoDataFile = annotFile,
                           phenoDataSheet = "Annotation",
                           phenoDataDccColName = "FileName",
                           protocolDataColNames = c("ROI", "AOI"),
                           experimentDataColNames = c("Panel")
)
saveRDS(gSet, file.path(outDir, "01_geomxSet.RDS"))
message("\n\t>> GeoMx data exported to an GeoMxSet object:")
message(paste0("\t  - Probes : ", dim(gSet)[1]))
message(paste0("\t  - Samples: ", dim(gSet)[2]))

# Export raw matrix
outMat <- assayData(gSet)$exprs
rownames(outMat) <- fData(gSet)$TargetName
colnames(outMat) <- pData(gSet)$Library
colnames(outMat)[1] <- paste0("Gene\t", colnames(outMat)[1])
write.table(outMat, file.path(outDir, "01_RAW_count.txt"), row.names=T, col.names=T, quote=F, sep="\t")
message("\n\t++ Created:")
message("\t  - Results/01_RAW_count.txt")

# Study design in a glance
countDf <- pData(gSet) %>% make_long(Region, Segment, Patient)
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

pdf(file.path(outDir, "01_Study_design.pdf"))
suppressWarnings(print(studyDesign))
invisible(dev.off())
message("\t  - Results/01_Study_design.txt")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 2 -- * --
# (1) Shift all zero values to one to transform in downstream analysis
gSet <- shiftCountsOne(gSet, useDALogic = TRUE)
saveRDS(gSet, file.path(outDir, "02_geomxSet.RDS"))
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 02: Shift any expression counts with a value of zero (0) to")
message("                 one (1) to enable in downstream transformations.")

outMat <- assayData(gSet)$exprs
rownames(outMat) <- fData(gSet)$TargetName
colnames(outMat) <- pData(gSet)$Library
colnames(outMat)[1] <- paste0("Gene\t", colnames(outMat)[1])

write.table(outMat, file.path(outDir, "02_RAW_noZero.txt"), row.names=T, col.names=T, quote=F, sep="\t")
message("\n\t++ Created:")
message("\t  - Results/02_RAW_noZero.txt")

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

tmpMat <- sData(gSet)
tmpMat$Segment <- factor(tmpMat$Segment, levels = segmentOrder)
suppressWarnings({
        qc1 <- histQC(tmpMat, "Trimmed (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_trimmedThre, cols = segmentCols)
        qc2 <- histQC(tmpMat, "Stitched (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_stitchedThre, cols = segmentCols)
        qc3 <- histQC(tmpMat, "Aligned (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_alignedThre, cols = segmentCols)
        qc4 <- histQC(tmpMat, "Saturated (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_saturatedThre, cols = segmentCols)
        qc5 <- histQC(tmpMat, "area", segmentQC_colBy, segmentQC_rowBy, segmentQC_areaThre, "log10", "AOI Area (log10)", cols = segmentCols)
        qc6 <- histQC(tmpMat, "nuclei", segmentQC_colBy, segmentQC_rowBy, segmentQC_nucleiThre, "log10", "AOI nuclei count", cols = segmentCols)
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
saveRDS(gSet, file.path(outDir, "03_geomxSet.RDS"))
message("\n\t>> Removed flagged samples(segments):")
message(paste0("\t  - Probes : ", dim(gSet)[1]))
message(paste0("\t  - Samples: ", dim(gSet)[2]))

tmpMat <- sData(gSet)
tmpMat$Segment <- factor(tmpMat$Segment, levels = segmentOrder)
suppressWarnings({
        qc1 <- histQC(tmpMat, "Trimmed (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_trimmedThre, cols = segmentCols)
        qc2 <- histQC(tmpMat, "Stitched (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_stitchedThre, cols = segmentCols)
        qc3 <- histQC(tmpMat, "Aligned (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_alignedThre, cols = segmentCols)
        qc4 <- histQC(tmpMat, "Saturated (%)", segmentQC_colBy, segmentQC_rowBy, segmentQC_saturatedThre, cols = segmentCols)
        qc5 <- histQC(tmpMat, "area", segmentQC_colBy, segmentQC_rowBy, segmentQC_areaThre, "log10", "AOI Area (log10)", cols = segmentCols)
        qc6 <- histQC(tmpMat, "nuclei", segmentQC_colBy, segmentQC_rowBy, segmentQC_nucleiThre, "log10", "AOI nuclei count", cols = segmentCols)
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
saveRDS(gSet, file.path(outDir, "04_geomxSet.RDS"))
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 04: Calculate negative count which is the geometric mean of")
message("                 the several unique negative probes in the GeoMx panel that")
message("                 do not target mRNA and establish the background count level")
message("                 per segment.")

tmpMat <- sData(gSet)
tmpMat <- tmpMat[,c(negCol, segmentQC_colBy, segmentQC_rowBy)]
tmpMat$Segment <- factor(tmpMat$Segment, levels = segmentOrder)
qc <- histQC(tmpMat, negCol, segmentQC_colBy, segmentQC_rowBy, 2, "log10", "GeoMean(negative probes)", cols = segmentCols)
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
saveRDS(gSet, file.path(outDir, "05_geomxSet.RDS"))

message("\n\t>> Removed flagged probes:")
message(paste0("\t  - Probes : ", dim(gSet)[1]))
message(paste0("\t  - Samples: ", dim(gSet)[2]))

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 6 -- * --
# (1) Aggregate multi-probe and generate a gene expression profile
newSet <- aggregateCounts(gSet)

#Exclude negative control from the gene list
newSet <- subset(newSet, fData(newSet)$TargetName != "NegProbe-WTX")

message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 06: Generate a gene-level count matrix where the count for")
message("                 any gene with multiple probes per segment is calculated as")
message("                 the geometric mean of those probes.")
message("\n\t>> Collapse to targets/genes:")
message(paste0("\t  - Genes  : ", dim(newSet)[1]))
message(paste0("\t  - Samples: ", dim(newSet)[2]))

saveRDS(newSet, file.path(outDir, "06_geomxSet.RDS"))
saveRDS(as.data.frame(pData(newSet)), file.path(outDir, "00_finalStats.RDS"))

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 7 -- * --
# (1) Calculate LOQ
# (2) Filter out low signal segments/genes compared to the background
loqDf <- data.frame(row.names = colnames(newSet))
tmpVar <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
if(all(tmpVar[1:2] %in% colnames(pData(newSet)))) {
        loqDf[, module] <- pData(newSet)[, tmpVar[1]] * (pData(newSet)[, tmpVar[2]] ^ loqCutoff)
}
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 07-1: Determine the limit of quantification (LOQ) per segment")
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
        violin <- ggplot(data = statDf, aes(x = Segment, y = log10(LOQ), fill = Segment)) +
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
                geom_hline(aes(yintercept = log10(2)), lty=2, col="grey50")
})

pdf(file.path(outDir, "07_LOQ_perSegment.pdf"))
print(violin)
invisible(dev.off())
message("\n\t++ Created:")
message("\t  - Results/07_LOQ_perSegment.pdf")

violin <- ggplot(data = statDf, aes(x = Segment, y = log10(LOQ), fill = Segment)) +
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
        geom_hline(aes(yintercept = log10(2)), lty=2, col="grey50")
 
pdf(file.path(outDir, "07_LOQ_perSlideSegment.pdf"))
print(violin)
invisible(dev.off())
message("\t  - Results/07_LOQ_perSlideSegment.pdf")

qc <- ggplot(data=statDf, aes(x=Sample, y=log10(LOQ), group=Segment)) +
        geom_line(aes(color=Segment), lwd=0.5) +
        geom_point(aes(color=Segment)) +
        scale_color_brewer(palette="Dark2") + 
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
        geom_hline(aes(yintercept = log10(2)), lty=2, col="grey50")
 
pdf(file.path(outDir, "07_LOQ_perSample.pdf"), width=15)
print(qc)
invisible(dev.off())
message("\t  - Results/07_LOQ_perSample.pdf")

loqDf[loqDf < 2] <- 2
pData(newSet)$LOQ <- loqDf

loqMat <- t(esApply(newSet, MARGIN = 1, FUN = function(x) { x > LOQ[, module] }) )
loqMat <- loqMat[fData(newSet)$TargetName, ] # Ordering

# Segment gene detection
pData(newSet)$GenesDetected <- colSums(loqMat, na.rm = TRUE)
pData(newSet)$GeneDetectionRate <- pData(newSet)$GenesDetected / nrow(newSet)
pData(newSet)$DetectionThreshold <- cut(pData(newSet)$GeneDetectionRate, breaks = geneDetectionRate, labels = geneDetectionRateLabel)
saveRDS(newSet, file.path(outDir, "07_geomxSet_beforeFiltering.RDS"))
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 07-2: Filter out segments with exceptionally low signal that")
message("                 have a small fraction of panel genes detected above the LOQ")
message("                 relative to the other segments in the study.")

write.table(paste0("# Gene detection rate threshold > ", geneDetectionRateThre), file.path(outDir, "07_Segment_gene_detection.txt"), row.names=F, col.names=F, quote=F, sep="\t")
suppressWarnings(write.table(t(table(pData(newSet)$DetectionThreshold)), file.path(outDir, "07_Segment_gene_detection.txt"), row.names=F, col.names=T, quote=F, sep="\t", append=T))
message("\n\t++ Created:")
message("\t  - Results/07_Segment_gene_detection.txt")

tmpMat <- pData(newSet)
tmpMat$Segment <- factor(tmpMat$Segment, levels = segmentOrder)
suppressWarnings({
        barplot <- ggplot(tmpMat, aes(x = DetectionThreshold)) +
                geom_bar(aes(fill = Segment)) +
                geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
                scale_fill_manual(values = segmentCols) +
                theme_bw() +
                theme(
                        axis.line = element_line(colour = "black"),
                        # panel.grid.major = element_blank(),
                        # panel.grid.minor = element_blank(),
                        # panel.border = element_blank(),
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

pdf(file.path(outDir, "07_detectedGenes_summary.pdf"))
print(barplot)
invisible(dev.off())
message("\t  - Results/07_detectedGenes_summary.txt")

statDf <- data.frame(
        Sample = protocolData(newSet)$AOI, # pData(newSet)$Sample,
        Segment = pData(newSet)$Segment,
        GenesDetected = colSums(loqMat, na.rm = TRUE),
        Nuclei = pData(newSet)$nuclei
)
statDf$Segment <- factor(statDf$Segment, segmentOrder)

pdf(file.path(outDir, "07_detectedGenes.pdf"))
for (segment in segmentOrder) {
        tmpDf <- statDf[which(statDf$Segment == segment),]
        tmpDf <- tmpDf[order(tmpDf$GenesDetected, decreasing = F),]
        tmpDf$Sample <- factor(tmpDf$Sample, levels=tmpDf$Sample)
        barplot <- ggplot(tmpDf, aes(x = Sample, y = GenesDetected, fill = Segment)) +
                geom_bar(stat = "identity") +
                # geom_line(aes(x = Sample, y = Nuclei), size = 1, group = 1) +
                # scale_y_continuous(
                #         name = "Number of Genes detected",
                #         sec.axis = sec_axis(trans = ~ ., name = "Number of Nuclei")
                # ) +
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
        print(barplot)
}
invisible(dev.off())
message("\t  - Results/07_detectedGenes.txt")

newSet <- newSet[, pData(newSet)$GeneDetectionRate >= geneDetectionRateThre]
message("\n\t>> After filtering out low signal sample/segment based on LOQ:")
message(paste0("\t  - Genes  : ", dim(newSet)[1]))
message(paste0("\t  - Samples: ", dim(newSet)[2]))
saveRDS(newSet, file.path(outDir, "07_geomxSet_afterFiltering.RDS"))

# Gene detection rate
loqMat <- loqMat[, colnames(newSet)]
fData(newSet)$DetectedSegments <- rowSums(loqMat, na.rm = TRUE)
fData(newSet)$DetectionRate <- fData(newSet)$DetectedSegments / nrow(pData(newSet))
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 07-3: Determine the detection rate for genes across the study")
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
write.table(geneDetectionRateDf, file.path(outDir, "07_Gene_detection_rate.txt"), row.names=F, col.names=T, quote=F, sep="\t")

negativeProbefData <- subset(fData(newSet), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
finalSet <- newSet[fData(newSet)$DetectionRate >= geneDetectionRateThre | fData(newSet)$TargetName %in% neg_probes, ]
saveRDS(finalSet, file.path(outDir, "07_geomxSet.RDS"))
message("\n\t>> After filtering out low signal gene:")
message(paste0("\t  - Genes  : ", dim(finalSet)[1]))
message(paste0("\t  - Samples: ", dim(finalSet)[2]))

countDf <- pData(gSet) %>% make_long(Region, Segment, Patient)
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

pdf(file.path(outDir, "07_Study_design.pdf"))
suppressWarnings(print(studyDesign))
invisible(dev.off())
message("\n\t++ Created:")
message("\t  - Results/07_Gene_detection_rate.txt")
message("\t  - Results/07_Study_design.pdf")

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 8 -- * --
# (1) Normalization
finalSet <- normalize(finalSet, norm_method = "quant", desiredQuantile = .75, toElt = "q_norm")
# finalSet <- normalize(finalSet, norm_method = "neg", fromElt = "exprs", toElt = "neg_norm") # Background normalization for WTA without custom spike-in
saveRDS(finalSet, file.path(outDir, "08_geomxSet.RDS"))
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 08: Upper-quartile (Q3) normalization method estimates")
message("                 a normalization factor per segment to bring the segment")
message("                 data distributions together.")

outMat <- assayData(finalSet)$exprs
colnames(outMat) <- pData(finalSet)$Library
write.table(outMat, file.path(outDir, "08_exprMat.txt"), row.names=T, col.names=T, quote=F, sep="\t")
message("\n\t++ Created:")
message("\t  - Results/08_exprMat.txt")

outMat <- assayData(finalSet)$q_norm
colnames(outMat) <- pData(finalSet)$Library
write.table(outMat, file.path(outDir, "08_exprMat_UQ.txt"), row.names=T, col.names=T, quote=F, sep="\t")
message("\t  - Results/08_exprMat_UQ.txt")

# library(vsn)
# vsnNorm <- vsn2(outMat)@hx
# write.table(vsnNorm, file.path(outDir, "08_exprMat_UQ_VSN.txt"), row.names=T, col.names=T, quote=F, sep="\t")
# message("  - Output/08_exprMat_UQ.txt")

# To compare before & after
noFilterSet <- readRDS(file.path(outDir, "07_geomxSet_beforeFiltering.RDS"))
noFilterSet <- normalize(noFilterSet, norm_method = "quant", desiredQuantile = .75, toElt = "q_norm")
saveRDS(noFilterSet, file.path(outDir, "08_geomxSet_noFiltering.RDS"))

##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# -- * -- Step 9 -- * --
# (1) Before & after diagnostic plots
finalSetExpr <- apply(assayData(finalSet)$q_norm, 2, median)
noFilterExpr <- apply(assayData(noFilterSet)$q_norm, 2, median)
message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Preprocessing 09: Generate diagnostic plots to investigate before and")
message("                 after QC.")

segmentNames <- list(All = rownames(pData(noFilterSet)), FilterIn = rownames(pData(finalSet)), FilterOut = setdiff(rownames(pData(noFilterSet)), rownames(pData(finalSet))))
write.table(pasteListValues(sapply(segmentNames, length), FALSE), file.path(outDir, "09_Stats.txt"), row.names=F, col.names=F, quote=F)
message("\n\t++ Created:")
message("\t  - Results/09_Stats.txt")

filteredOut <- setdiff(rownames(pData(noFilterSet)), rownames(pData(finalSet)))
if (!identical(filteredOut, character(0))) {
        stats <- pData(noFilterSet)
        df <- data.frame(
                Slide = stats$Slide,
                Sample = stats$Sample,
                Segment = stats$Segment,
                Area = stats$area,
                Nuclei = stats$nuclei,
                LOQ = stats$LOQ,
                GenesDetected = stats$GenesDetected,
                GenesDetRate = stats$GeneDetectionRate,
                Filter = 0
        )
        colnames(df)[6] <- "LOQ"
        df$Filter[which(rownames(df) %in% colnames(finalSet))] <- 1
        df$Slide <- factor(df$Slide, levels=slideOrder)
        df$Segment <- factor(df$Segment, levels=segmentOrder)

        for (segment in segmentOrder) {
                subDf <- df[which(df$Segment == segment),]
                table <- xtabs( ~ Filter + Slide, data = subDf)
                fisherObj <- try(fisher.test(table), silent=TRUE)
                if(is(fisherObj, "try-error")) { fisherP <- NA } else { fisherP <- formatC(fisherObj$p.value, format="e", digits=2) }
                cattObj <- try(CochranArmitageTest(table, alternative = c("two.sided")), silent=TRUE)
                if(is(cattObj, "try-error")) { cattP <- NA } else { cattP <- formatC(cattObj$p.value, format="e", digits=2) }
                rowSum <- apply(table, 1, sum)
                table <- cbind(table, rowSum)
                colSum <- apply(table, 2, sum)
                table <- rbind(table, colSum)
                colnames(table)[1] <- paste0("Filter\t", colnames(table)[1])
                colnames(table)[ncol(table)] <- "TOTAL"
                rownames(table) <- c("Out", "In", "TOTAL")

                write.table(paste0("# Fisher Exact Test, p-value = ", fisherP, "\n# Cochran-Armitage Trend Test, p-value = ", cattP), file.path(outDir, paste0("09_Stats_", segment, ".txt")), row.names=F, col.names=F, quote=F, sep="\t")
                suppressWarnings(write.table(table, file.path(outDir, paste0("09_Stats_", segment, ".txt")), row.names=T, col.names=T, quote=F, sep="\t", append=T))
                message(paste0("\t  - Results/09_Stats_", segment, ".txt"))
        }

        dummies <- sapply(seq_along(segmentNames), function(idx) {        
                prefix <- names(segmentNames)[idx]
                selectSegments <- segmentNames[[idx]]

                if (length(selectSegments) > 0) {
                        subDf <- df[which(rownames(df) %in% selectSegments),]
                        subDf <- data.frame(subDf, MedianExpression = noFilterExpr[selectSegments])
                        
                        b1 <- boxQC(subDf, "Segment", "Area", "Slide", "log10", cols = segmentCols)
                        b2 <- boxQC(subDf, "Segment", "Nuclei", "Slide", "log10", cols = segmentCols)
                        b3 <- boxQC(subDf, "Segment", "LOQ", "Slide", "log10", cols = segmentCols)
                        b4 <- boxQC(subDf, "Segment", "MedianExpression", "Slide", "log10", cols = segmentCols, yLabel = "Median expression")
                        
                        pdf(file.path(outDir, paste0("09_Boxplot_", prefix, ".pdf")))
                        print(b1)
                        print(b2)
                        print(b3)
                        print(b4)
                        invisible(dev.off())
                        message(paste0("\t  - Results/09_Boxplot_", prefix, ".pdf"))

                        pdf(file.path(outDir, paste0("09_Scatterplot_", prefix, ".pdf")))
                        for (segment in segmentOrder) {
                                suppressWarnings({
                                        segDf <- subDf[which(subDf$Segment == segment), c("Area", "Nuclei", "LOQ", "MedianExpression")]
					if (nrow(segDf) > 3) {	
                                        	pairs(segDf, lower.panel=panel.cor, upper.panel=panel.lm, pch=20, cex=2, main=segment, col=segmentCols[segment])
					}
                                })
                        }
                        invisible(dev.off())
                        message(paste0("\t  - Results/09_Scatterplot_", segment, ".pdf"))
                }

                return(idx) # dummy
        })

        prefix <- "FilterIn_UQ3"
        selectSegments <- segmentNames[[2]]

        subDf <- df[which(rownames(df) %in% selectSegments),]
        subDf <- data.frame(subDf, MedianExpression = finalSetExpr)

        b1 <- boxQC(subDf, "Segment", "Area", "Slide", "log10", cols = segmentCols)
        b2 <- boxQC(subDf, "Segment", "Nuclei", "Slide", "log10", cols = segmentCols)
        b3 <- boxQC(subDf, "Segment", "LOQ", "Slide", "log10", cols = segmentCols)
        b4 <- boxQC(subDf, "Segment", "MedianExpression", "Slide", "log10", cols = segmentCols, yLabel = "Median expression")

        pdf(file.path(outDir, paste0("09_Boxplot_", prefix, ".pdf")))
        print(b1)
        print(b2)
        print(b3)
        print(b4)
        invisible(dev.off())
        message(paste0("\t  - Results/09_Boxplot_", prefix, ".pdf"))

        pdf(file.path(outDir, paste0("09_Scatterplot_", prefix, ".pdf")))
        for (segment in segmentOrder) {
                suppressWarnings({
                        segDf <- subDf[which(subDf$Segment == segment), c("Area", "Nuclei", "LOQ", "MedianExpression")]
			if (nrow(segDf) > 3) {
                        	pairs(segDf, lower.panel=panel.cor, upper.panel=panel.lm, pch=20, cex=2, main=segment, col=segmentCols[segment])
			}
                })
        }
        invisible(dev.off())
        message(paste0("\t  - Results/09_Scatterplot_", prefix, ".pdf"))
}

message("\n##### ##### ##### ##### ##### ##### ##### ##### ##### #####")
message("Fin")

q("no")
