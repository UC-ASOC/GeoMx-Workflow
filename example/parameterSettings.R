##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Heewon Seo (Heewon.Seo@UCalgary.ca)
# Updated on Oct 23, 2023
# SOA Brain
##### ##### ##### ##### ##### ##### ##### ##### ##### #####
# Parameter settings
slideOrder <- c("74", "77", "73-1", "hu_brain_004a", "hu_brain_004b")

segmentOrder <- c("Neuron", "Iba1", "GFAP", "Neuropil", "Full")
segmentCols <- c("magenta", "gold", "cyan", "indianred2", "chartreuse1")
names(segmentCols) <- segmentOrder

segmentQcParams <- list(
        minSegmentReads = 1000, # Minimum number of reads (1000)
        percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
        percentStitched = 80,   # Minimum % of reads stitched (80%)
        percentAligned = 75,    # Minimum % of reads aligned (80%)
        percentSaturation = 50, # Minimum sequencing saturation (50%)
        minNegativeCount = 1,   # Minimum negative control counts (10)
        maxNTCCount = 9000,     # Maximum counts observed in NTC well (1000)
        minNuclei = 20,         # Minimum # of nuclei estimated (100)
        minArea = 1000          # Minimum segment area (5000)
)

segmentQC_colBy <- "Segment"
segmentQC_rowBy <- "Slide"
segmentQC_trimmedThre <- 80     # per cent
segmentQC_stitchedThre <- 80    # per cent
segmentQC_alignedThre <- 75     # per cent
segmentQC_saturatedThre <- 50   # per cent
segmentQC_areaThre <- 1000      # micro meter square
segmentQC_nucleiThre <- 20      # nuclei count

probeQcParams <- list(
        minProbeRatio = 0.1,    # geometric mean of that probe’s counts from all segments divided by the geometric mean of all probe counts representing the target from all segments is less than 0.1
        percentFailGrubbs = 20  # the probe is an outlier according to the Grubb’s test in at least 20% of the segments
)

loqCutoff <- 2
loqMin <- 2

geneDetectionRate <- c(0, 0.01, 0.05, 0.1, 0.15, 1) # 1%, 5%, 10%, 15%, >15%
geneDetectionRateLabel <- c("<1%", "1-5%", "5-10%", "10-15%", ">15%")
geneDetectionRateThre <- 0.05 # Nanostring recommended 5% or 10%

message("\nParameters loaded.")