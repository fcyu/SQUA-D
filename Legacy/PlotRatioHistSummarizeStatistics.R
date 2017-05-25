# Plot some histogram.
# Summarize mean and variance of log2 ratios

library(stats)
library(graphics)

rm(list = ls())

inputPath <- "D:\\Results\\Mascot\\Acetylation\\summary.UPS.FR.tsv"
resultDir <- "D:\\Results\\Mascot\\Acetylation\\summary.UPS.FR.tsv.ratios.histograms\\"
dataPointNumT <- 0

inputTable <- read.table(inputPath, header = TRUE, sep = "\t", fill = TRUE, comment.char = "", quote = "")

f1Means <- numeric()
f1Sds <- numeric()
f2Means <- numeric()
f2Sds <- numeric()
r1Means <- numeric()
r1Sds <- numeric()
r2Means <- numeric()
r2Sds <- numeric()

f1AllRatios <- numeric()
f2AllRatios <- numeric()
r1AllRatios <- numeric()
r2AllRatios <- numeric()

cat("Summarizing...\n")
i <- 0
lastProgress <- 0
for (rowIdx in which(inputTable$Pass.XIC.criteria)) {
    i <- i + 1
    progress = floor(i / 100)
    if (progress > lastProgress) {
        cat(sprintf("Summarizing %d of %d...\n", progress * 100, sum(inputTable$Pass.XIC.criteria)))
        lastProgress <- progress
    }

    myRow <- inputTable[rowIdx,]
    ratioTable <- read.table(paste(inputPath, ".ratios\\", myRow$UPS, ".tsv", sep = ""), header = TRUE, sep = "\t", comment.char = "", quote = "")
    f1Ratios <- as.numeric(as.character(ratioTable[grepl("Batch1_.F.", ratioTable$fraction), 2]))
    f2Ratios <- as.numeric(as.character(ratioTable[grepl("Batch2_.F.", ratioTable$fraction), 2]))
    r1Ratios <- as.numeric(as.character(ratioTable[grepl("Batch1_.R.", ratioTable$fraction), 2]))
    r2Ratios <- as.numeric(as.character(ratioTable[grepl("Batch2_.R.", ratioTable$fraction), 2]))

    if ((length(f1Ratios) >= dataPointNumT) & (length(f2Ratios) >= dataPointNumT) & (length(r1Ratios) >= dataPointNumT) & (length(r2Ratios) >= dataPointNumT)) {
        f1Means <- c(f1Means, mean(f1Ratios))
        f1Sds <- c(f1Sds, sd(f1Ratios))
        f2Means <- c(f2Means, mean(f2Ratios))
        f2Sds <- c(f2Sds, sd(f2Ratios))
        r1Means <- c(r1Means, mean(r1Ratios))
        r1Sds <- c(r1Sds, sd(r1Ratios))
        r2Means <- c(r2Means, mean(r2Ratios))
        r2Sds <- c(r2Sds, sd(r2Ratios))

        f1AllRatios <- c(f1AllRatios, f1Ratios)
        f2AllRatios <- c(f2AllRatios, f2Ratios)
        r1AllRatios <- c(r1AllRatios, r1Ratios)
        r2AllRatios <- c(r2AllRatios, r2Ratios)
    }
}

xLim <- range(c(f1AllRatios, f2AllRatios, r1AllRatios, r2AllRatios))
layout(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2, byrow = TRUE))
f1q <- quantile(f1AllRatios, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
f2q <- quantile(f2AllRatios, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
r1q <- quantile(r1AllRatios, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
r2q <- quantile(r2AllRatios, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
hist(f1AllRatios, breaks = 100, xlim = xLim, main = sprintf("mean: %f; std: %f; \n quantile: %f, %f, %f, %f, %f", mean(f1AllRatios), sd(f1AllRatios), f1q[1], f1q[2], f1q[3], f1q[4], f1q[5]))
hist(f2AllRatios, breaks = 100, xlim = xLim, main = sprintf("mean: %f; std: %f; \n quantile: %f, %f, %f, %f, %f", mean(f2AllRatios), sd(f2AllRatios), f2q[1], f2q[2], f2q[3], f2q[4], f2q[5]))
hist(r1AllRatios, breaks = 100, xlim = xLim, main = sprintf("mean: %f; std: %f; \n quantile: %f, %f, %f, %f, %f", mean(r1AllRatios), sd(r1AllRatios), r1q[1], r1q[2], r1q[3], r1q[4], r1q[5]))
hist(r2AllRatios, breaks = 100, xlim = xLim, main = sprintf("mean: %f; std: %f; \n quantile: %f, %f, %f, %f, %f", mean(r2AllRatios), sd(r2AllRatios), r2q[1], r2q[2], r2q[3], r2q[4], r2q[5]))

cat("Done!")
