# extract ion chromatogram

rm(list = ls())

options(java.parameters = "-Xmx32g")

library(xcms)
library(stats)
library(caTools)
library(signal)
library(xlsx)
library(stringr)


digitizeMzIntensity <- function(mzIntensity, binSize, minMz, maxMz) {
    outputVector <- rep(0, ceiling((maxMz - minMz) / binSize) + 1)
    for (i in seq(1, nrow(mzIntensity))) {
        outputVector[floor((mzIntensity[i, 1] - minMz) / binSize) + 1] <- outputVector[floor((mzIntensity[i, 1] - minMz) / binSize) + 1] + mzIntensity[i, 2]
    }
    return(outputVector)
}


## A function to plot XIC
plotXicMs1 <- function(xcmsRawObj, ppmTol, mzC, scanNumC, rtC, rtCLeft, rtCRight, mzT, scanNumT, rtT, rtTLeft, rtTRight, plotDir, runName, charge, adjustRatio, firstUpsp) {
    C13Step <- 1.00727646688 / charge

    # get control XIC
    rtIntensityC <- NA
    if (!is.na(rtC)) {
        mzRangeC <- matrix(c(mzC * (1 - ppmTol * 1e-6), mzC * (1 + ppmTol * 1e-6)), nrow = 1, ncol = 2)
        rtRangeC <- matrix(c(rtCLeft, rtCRight), nrow = 1, ncol = 2)
        eicObj <- getEIC(xcmsRawObj, mzrange = mzRangeC, rtrange = rtRangeC)
        rtIntensityC <- eicObj@eic$xcmsRaw[[1]]
    }

    # get treatment XIC
    rtIntensityT <- NA
    if (!is.na(rtT)) {
        mzRangeT <- matrix(c(mzT * (1 - ppmTol * 1e-6), mzT * (1 + ppmTol * 1e-6)), nrow = 1, ncol = 2)
        rtRangeT <- matrix(c(rtTLeft, rtTRight), nrow = 1, ncol = 2)
        eicObj <- getEIC(xcmsRawObj, mzrange = mzRangeT, rtrange = rtRangeT)
        rtIntensityT <- eicObj@eic$xcmsRaw[[1]]
        # adjust treatment intensity based on fixing error
        rtIntensityT[, 2] <- rtIntensityT[, 2] * adjustRatio
    }

    minRt <- min(rtCLeft, rtTLeft, na.rm = TRUE)
    maxRt <- max(rtCRight, rtTRight, na.rm = TRUE)

    # get MS1 peaks
    mzIntensityC <- NA
    mzIntensityT <- NA
    if (!is.na(rtC)) {
        idx <- which.min(abs(rtC - xcmsRawObj@scantime)) # the index is not scan num
        mzIntensityC <- getScan(xcmsRawObj, idx, c(mzC - 0.3, mzC + 3 * C13Step + 0.3))
    }
    if (!is.na(rtT)) {
        idx <- which.min(abs(rtT - xcmsRawObj@scantime)) # the index is not scan num
        mzIntensityT <- getScan(xcmsRawObj, idx, c(mzT - 0.3, mzT + 3 * C13Step + 0.3))
        # adjust treatment intensity based on fixing error
        mzIntensityT[, 2] <- mzIntensityT[, 2] * adjustRatio
    }

    # plot and save the figure
    jpeg(paste(plotDir, firstUpsp, "-", runName, "-", scanNumC, "-", scanNumT, ".jpg", sep = ""), width = 18.75, height = 22.5, units = "cm", res = 300)
    par(mfrow = c(2, 1), yaxs = "i")

    ColorC <- "blue"
    ColorT <- "red"

    if (!is.na(rtC) && !is.na(rtT)) {
        plot(rtIntensityC[, 1], rtIntensityC[, 2], main = "Extracted Ion Chromatogram", xlab = "Retention Time (second)", ylab = "Intensity", xlim = c(minRt, maxRt), ylim = c(0, max(rtIntensityC[, 2], rtIntensityT[, 2]) * 1.1), type = "l", col = ColorC, lwd = 2)
        par(new = TRUE)
        plot(rtIntensityT[, 1], rtIntensityT[, 2], main = "Extracted Ion Chromatogram", xlab = "Retention Time (second)", ylab = "Intensity", xlim = c(minRt, maxRt), ylim = c(0, max(rtIntensityC[, 2], rtIntensityT[, 2]) * 1.1), type = "l", col = ColorT, lwd = 2)
        par(new = FALSE)
        myTitle <- sprintf("MS1, RT(control) = %ds, RT(treatment) = %ds", round(rtC), round(rtT))
        myXlim <- c(min(mzC, mzT) - 0.3, max(mzC, mzT) + 2 * C13Step + 0.3)
        myYlim <- c(0, max(mzIntensityC[, 2], mzIntensityT[, 2]) * 1.1)
        plot(mzIntensityC[, 1], mzIntensityC[, 2], main = myTitle, xlab = "MZ (Th)", ylab = "Intensity", xlim = myXlim, ylim = myYlim, type = "h", col = ColorC, lwd = 2)
        par(new = TRUE)
        plot(mzIntensityT[, 1], mzIntensityT[, 2], main = myTitle, xlab = "MZ (Th)", ylab = "Intensity", xlim = myXlim, ylim = myYlim, type = "h", col = ColorT, lwd = 2)
        par(new = FALSE)
    } else if (!is.na(rtC)) {
        plot(rtIntensityC[, 1], rtIntensityC[, 2], main = "Extracted Ion Chromatogram", xlab = "Retention Time (second)", ylab = "Intensity", xlim = c(minRt, maxRt), ylim = c(0, max(rtIntensityC[, 2]) * 1.1), type = "l", col = ColorC, lwd = 2)
        myTitle <- sprintf("MS1, RT(control) = %ds, RT(treatment) = NAs", round(rtC))
        myXlim <- c(min(mzC, mzT) - 0.3, max(mzC, mzT) + 2 * C13Step + 0.3)
        myYlim <- c(0, max(mzIntensityC[, 2]) * 1.1)
        plot(mzIntensityC[, 1], mzIntensityC[, 2], main = myTitle, xlab = "MZ (Th)", ylab = "Intensity", xlim = myXlim, ylim = myYlim, type = "h", col = ColorC, lwd = 2)
    } else if (!is.na(rtT)) {
        plot(rtIntensityT[, 1], rtIntensityT[, 2], main = "Extracted Ion Chromatogram", xlab = "Retention Time (second)", ylab = "Intensity", xlim = c(minRt, maxRt), ylim = c(0, max(rtIntensityT[, 2]) * 1.1), type = "l", col = ColorT, lwd = 2)
        myTitle <- sprintf("MS1, RT(control) = NAs, RT(treatment) = %ds", round(rtT))
        myXlim <- c(min(mzC, mzT) - 0.3, max(mzC, mzT) + 2 * C13Step + 0.3)
        myYlim <- c(0, max(mzIntensityT[, 2]) * 1.1)
        plot(mzIntensityT[, 1], mzIntensityT[, 2], main = myTitle, xlab = "MZ (Th)", ylab = "Intensity", xlim = myXlim, ylim = myYlim, type = "h", col = ColorT, lwd = 2)
    } else {
        stop(sprintf("Something wrong in runName = %s, scanNumC = %d, scanNumT = %d, UPSP = %s", runName, scanNumC, scanNumT, firstUpsp))
        dev.off()
    }
    dev.off()
}


# main part
cat(str_c(Sys.time(), ": Setting parameters...\n"))
ppmTol <- 20 # twice of the MS1 tolerance from including the whole peak
fdrT <- 0.05
upspPath <- "D:\\Results\\Quantification\\Acetylation\\12batches\\Pretable3.xlsx"
psmPath <- "D:\\Results\\Quantification\\Acetylation\\12batches\\filtered_psm_table.tsv"
spectraDir <- "D:\\Li_group\\Acetylation_raw_data\\"

cat(str_c(Sys.time(), ": Creating a dir holding dirs of plots...\n"))
plotsDir <- paste(upspPath, ".plot\\", sep = "")
if (!dir.exists(plotsDir)) {
    dir.create(plotsDir)
}

cat(str_c(Sys.time(), ": Reading the table...\n"))
upspTable <- read.xlsx2(upspPath, 1, header = TRUE, colClasses = c(rep("character", 2), rep("integer", 32), rep("numeric", 6)))
psmTable <- read.table(psmPath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, na.strings = c("NA", ""), comment.char = "", quote = "", colClasses = c(rep("character", 3), "integer", rep("numeric", 4), "integer", rep("character", 2), "integer", rep("character", 3), rep("numeric", 3), "character", "character", rep("numeric", 4), "integer", rep("numeric", 6)))
logRatioT <- sd(upspTable$adjusted_log_ratio_mean)

# read spectra
cat(str_c(Sys.time(), ": Reading and ploting...\n"))
spectraFileNames <- sort(list.files(spectraDir, pattern = ".+\\.mzXML$", full.names = FALSE, no.. = TRUE))
lastProgress <- 0
for (i in seq(1, length(spectraFileNames))) {
    progress = floor(i * 100 / length(spectraFileNames))
    if (progress > lastProgress) {
        cat(sprintf("\r%d%%", progress))
        lastProgress <- progress
    }

    # read the mzXML file
    xcmsRawObj <- xcmsRaw(paste(spectraDir, spectraFileNames[i], sep = ""), profstep = 0)

    # read ratio table and plot
    for (j in which((upspTable$BH_FDR <= fdrT) & (abs(upspTable$adjusted_log_ratio_mean) >= logRatioT))) {
        subPsmTable <- psmTable[psmTable$UPSP == upspTable[j,]$UPSP,]

        if (nrow(subPsmTable) > 0) {
            firstUpsp <- str_split(upspTable[j,]$UPSP, ";")[[1]][1]

            # create a dir holding plots
            plotDir <- paste(plotsDir, firstUpsp, "\\", sep = "")
            if (!dir.exists(plotDir)) {
                dir.create(plotDir)
            }

            for (k in seq(1, nrow(subPsmTable))) {
                runName <- subPsmTable[k,]$run_name
                if (grepl(runName, spectraFileNames[i]) && !is.na(subPsmTable[k, ]$normalized_log_ratio)) {
                    charge <- subPsmTable[k,]$charge
                    group_type <- subPsmTable[k,]$group_type
                    scanNumC <- NA
                    rtC <- NA
                    rtCLeft <- NA
                    rtCRight <- NA
                    mzC <- NA
                    scanNumT <- NA
                    rtT <- NA
                    rtTLeft <- NA
                    rtTRight <- NA
                    mzT <- NA
                    if (group_type == "C") {
                        scanNumC <- subPsmTable[k,]$scan_num
                        rtC <- subPsmTable[k,]$rt_corrected
                        rtCLeft <- subPsmTable[k,]$rt_left
                        rtCRight <- subPsmTable[k,]$rt_right
                        mzC <- subPsmTable[k,]$spectrum_mz
                        scanNumT <- subPsmTable[k,]$another_scan_num
                        rtT <- subPsmTable[k,]$another_rt
                        rtTLeft <- subPsmTable[k,]$another_rt_left
                        rtTRight <- subPsmTable[k,]$another_rt_right
                        mzT <- subPsmTable[k,]$another_mz
                    } else {
                        scanNumT <- subPsmTable[k,]$scan_num
                        rtT <- subPsmTable[k,]$rt_corrected
                        rtTLeft <- subPsmTable[k,]$rt_left
                        rtTRight <- subPsmTable[k,]$rt_right
                        mzT <- subPsmTable[k,]$spectrum_mz
                        scanNumC <- subPsmTable[k,]$another_scan_num
                        rtC <- subPsmTable[k,]$another_rt
                        rtCLeft <- subPsmTable[k,]$another_rt_left
                        rtCRight <- subPsmTable[k,]$another_rt_right
                        mzC <- subPsmTable[k,]$another_mz
                    }

                    adjustRatio <- 2 ^ (subPsmTable[k,]$normalized_log_ratio - subPsmTable[k,]$log_ratio)

                    plotXicMs1(xcmsRawObj, ppmTol, mzC, scanNumC, rtC, rtCLeft, rtCRight, mzT, scanNumT, rtT, rtTLeft, rtTRight, plotDir, runName, charge, adjustRatio, firstUpsp)
                }
            }
        }
    }
}

cat(str_c("\n", Sys.time(), ": Done!\n"))
