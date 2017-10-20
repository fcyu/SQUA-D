# Adjust batch effect for log ratios
# It is based on ComBat. The major difference is that we don't have X term in equation (2.1) because all runs are under the same "experimental condition"--half from control and half from treatment. Y_{ijg} is log-ratios.
# Only batches that have more than two log ratios are adjusted.

suppressMessages(suppressWarnings(library(stats)))
suppressMessages(suppressWarnings(library(MASS)))
suppressMessages(suppressWarnings(library(plyr)))
suppressMessages(suppressWarnings(library(stringr)))


# A function of adjusting batch effect
adjust <- function(expData, batchNum, meanOnly = TRUE, plotPrior = FALSE) {
    iterPrecision <- 0.0001

    # Step 1: standardize the data
    stdData <- lapply(expData, step1)


    # Step 2
    gammaHat <- matrix(NA, length(stdData), batchNum)
    deltaHatSquare <- matrix(NA, length(stdData), batchNum)
    if (meanOnly) {
        deltaHatSquare <- matrix(1, length(stdData), batchNum)
    }
    ni <- rep(0, batchNum) # this is for Step 3
    for (g in seq(1, length(stdData))) {
        localStdData <- stdData[[g]]
        if ((length(localStdData) > 0) && (sum(localStdData$adjustableBatches > 0))) {
            Y <- localStdData$Y
            batchStructure <- localStdData$batchStructure
            adjustableBatches <- localStdData$adjustableBatches
            if (sum(adjustableBatches) == 1) {
                gammaHat[g, adjustableBatches] <- mean(Y[batchStructure[, adjustableBatches]], na.rm = TRUE)
                if (!meanOnly) {
                    deltaHatSquare[g, adjustableBatches] <- var(Y[batchStructure[, adjustableBatches]], na.rm = TRUE)
                }
                ni[adjustableBatches] <- ni[adjustableBatches] + 1
            } else {
                gammaHat[g, adjustableBatches] <- apply(batchStructure[, adjustableBatches], 2, function(x, y) { mean(y[x], na.rm = TRUE) }, y = Y)
                if (!meanOnly) {
                    deltaHatSquare[g, adjustableBatches] <- apply(batchStructure[, adjustableBatches], 2, function(x, y) { var(y[x], na.rm = TRUE) }, y = Y)
                }
                ni[adjustableBatches] <- ni[adjustableBatches] + colSums(batchStructure[, adjustableBatches])
            }
        }
    }
    gammaBar <- apply(gammaHat, 2, mean, na.rm = TRUE)
    tauBarSquare <- apply(gammaHat, 2, var, na.rm = TRUE)
    VBar <- apply(deltaHatSquare, 2, mean, na.rm = TRUE)
    SBarSquare <- apply(deltaHatSquare, 2, var, na.rm = TRUE)
    lammdaBar = (VBar + 2 * SBarSquare) / SBarSquare
    thetaBar = (VBar ^ 3 + VBar * SBarSquare) / SBarSquare

    if (plotPrior) {
        if (!meanOnly) {
            par(mfrow = c(2, 2))
        } else {
            par(mfrow = c(1, 2))
        }

        tmp <- density(gammaHat[, 1], na.rm = TRUE)
        plot(tmp, type = "l", main = "Density Plot")
        xx <- seq(min(tmp$x), max(tmp$x), length = 100)
        lines(xx, dnorm(xx, gammaBar[1], sqrt(tauBarSquare[1])), col = 2)
        qqnorm(gammaHat[, 1])
        qqline(gammaHat[, 1], col = 2)

        if (!meanOnly) {
            tmp <- density(deltaHatSquare[, 1], na.rm = TRUE)
            invgam <- 1 / rgamma(nrow(deltaHatSquare), lammdaBar[1], thetaBar[1])
            tmp1 <- density(invgam, na.rm = TRUE)
            plot(tmp, type = "l", main = "Density Plot", ylim = c(0, max(tmp$y, tmp1$y)))
            lines(tmp1, col = 2)
            qqplot(deltaHatSquare[, 1], invgam, xlab = "Sample Quantiles", ylab = "Theoretical Quantiles")
            lines(c(0, max(invgam)), c(0, max(invgam)), col = 2)
            title("Q-Q Plot")
        }
    }

    # Step 3
    gammaOld <- gammaHat
    deltaSquareOld <- deltaHatSquare
    change <- 1
    while (change > iterPrecision) {
        temp1 <- matrix(rep(ni * tauBarSquare, nrow(gammaHat)), nrow = nrow(gammaHat), byrow = TRUE)
        gammaNew <- (temp1 * gammaHat + deltaSquareOld * matrix(rep(gammaBar, nrow(gammaHat)), nrow = nrow(gammaHat), byrow = TRUE)) / (temp1 + deltaSquareOld)
        if (meanOnly) {
            deltaSquareNew <- matrix(1, nrow(deltaSquareOld), ncol(deltaSquareOld))
        } else {
            mySumRes <- matrix(unlist(lapply(seq(1, length(stdData)), mySum, stdData = stdData, gamma = gammaNew, batchNum = batchNum)), length(stdData), byrow = TRUE)
            deltaSquareNew <- (matrix(rep(thetaBar, length(stdData)), nrow = length(stdData), byrow = TRUE) + mySumRes) / matrix(rep(ni / 2 + lammdaBar - 1, length(stdData)), nrow = length(stdData), byrow = TRUE)
        }
        change <- max(apply(abs(gammaNew - gammaOld) / (gammaOld + iterPrecision), 2, max, na.rm = TRUE), apply(abs(deltaSquareNew - deltaSquareOld) / (deltaSquareOld + iterPrecision), 2, max, na.rm = TRUE), na.rm = TRUE)
        gammaOld <- gammaNew
        deltaSquareOld <- deltaSquareNew
    }


    # Step 4
    adjustedData <- lapply(seq(1, length(stdData)), step4, stdData = stdData, deltaStarSquare = deltaSquareOld, gammaStar = gammaOld)

    return(adjustedData)
}


# The function of step 1
step1 <- function(localData) {
    if (length(localData) > 0) {
        if (sum(localData$adjustableBatches) > 0) {
            # Since we don't have the matrix X, least square is unnecessary.
            Y <- localData$Y
            adjustableLogRatios <- NA
            if (sum(localData$adjustableBatches) == 1) {
                adjustableLogRatios = localData$batchStructure[, localData$adjustableBatches]
            } else {
                adjustableLogRatios = rowSums(localData$batchStructure[, localData$adjustableBatches]) > 0
            }
            alphaHat <- mean(Y[adjustableLogRatios], na.rm = TRUE)
            sigmaHat <- sd(Y[adjustableLogRatios], na.rm = TRUE)
            stdLogRatio <- Y
            stdLogRatio[adjustableLogRatios] <- (Y[adjustableLogRatios] - alphaHat) / sigmaHat
            return(list(Y = stdLogRatio, batchStructure = localData$batchStructure, adjustableBatches = localData$adjustableBatches, alphaHat = alphaHat, sigmaHat = sigmaHat))
        } else {
            return(list(Y = localData$Y, adjustableBatches = localData$adjustableBatches))
        }
    } else {
        return(list())
    }
}


mySum <- function(g, stdData, gamma, batchNum) {
    localStdData <- stdData[[g]]
    outputRow <- rep(NA, batchNum)
    if ((length(localStdData) > 0) && (sum(localStdData$adjustableBatches) > 0)) {
        Y <- localStdData$Y
        batchStructure <- localStdData$batchStructure
        adjustableBatches <- localStdData$adjustableBatches
        if (sum(adjustableBatches) == 1) {
            outputRow[adjustableBatches] <- 0.5 * sum((Y[batchStructure[, adjustableBatches]] - gamma[g, adjustableBatches]) ^ 2)
        } else {
            for (i in seq(1, ncol(batchStructure))) {
                if (adjustableBatches[i]) {
                    outputRow[i] <- 0.5 * sum((Y[batchStructure[, i]] - gamma[g, i]) ^ 2)
                }
            }
        }
    }
    return(outputRow)
}


step4 <- function(g, stdData, deltaStarSquare, gammaStar) {
    localStdData <- stdData[[g]]
    if (length(localStdData) > 0) {
        if (sum(localStdData$adjustableBatches) > 0) {
            Y <- localStdData$Y
            adjustedLocalData <- unlist(lapply(seq(1, length(Y)), myFun, sigmaHat = localStdData$sigmaHat, deltaStarSquareRow = deltaStarSquare[g,], Y = Y, batchStructure = localStdData$batchStructure, adjustableBatches = localStdData$adjustableBatches, gammaStarRow = gammaStar[g,])) + localStdData$alphaHat
            return(list(Y = adjustedLocalData))
        } else {
            return(list(Y = localStdData$Y))
        }
    } else {
        return(list())
    }
}


myFun <- function(rowIdx, sigmaHat, deltaStarSquareRow, Y, batchStructure, adjustableBatches, gammaStarRow) {
    colIdx <- which(batchStructure[rowIdx, adjustableBatches])
    if (length(colIdx) == 1) {
        return((sigmaHat / sqrt(deltaStarSquareRow[colIdx])) * (Y[rowIdx] - gammaStarRow[colIdx]))
    } else if (length(colIdx) == 0) {
        return(Y[rowIdx])
    } else {
        stop("Something wrong in step 4\n")
    }
}


createBatchStructure <- function(batchType, allBatchType) {
    bacthStructureRow <- rep(FALSE, length(allBatchType))
    bacthStructureRow[allBatchType == batchType] <- TRUE
    return(bacthStructureRow)
}


calPValue <- function(expDataEntry) {
    x <- expDataEntry$Y
    x <- x[is.finite(x)]
    if (length(x) > 2) {
        if (sd(x) > 0) {
            return(t.test(x, mu = 0)$p.value)
        } else {
            return(NA)
        }
    } else {
        return(NA)
    }
}


adjustBatchEffect <- function(upspTable, psmTable, debug = FALSE) {
    allBatchType <- str_sort(unique(psmTable$batch_type))
    expData <- list()
    for (i in seq(1, nrow(upspTable))) {
        myRow <- upspTable[i,]
        ratioTable <- psmTable[(psmTable$UPSP == upspTable[i,]$UPSP) & !is.na(psmTable$normalized_log_ratio),]
        logRatio <- ratioTable$normalized_log_ratio
        if (length(logRatio) > 1) {
            batchStructure <- matrix(unlist(lapply(ratioTable$batch_type, createBatchStructure, allBatchType)), length(logRatio), byrow = TRUE)
            colnames(batchStructure) <- allBatchType
            adjustableBatches <- colSums(batchStructure) > 3
            expData[[i]] <- list(Y = logRatio, batchStructure = batchStructure, adjustableBatches = adjustableBatches)
        } else {
            expData[[i]] <- list()
        }
    }

    adjustedData <- adjust(expData, length(allBatchType), meanOnly = TRUE, plotPrior = debug)

    upspTable["Adjusted_log_ratio_mean"] <- unlist(lapply(adjustedData, function(x) { mean(x$Y, na.rm = TRUE) }))
    upspTable["Adjusted_log_ratio_sd"] <- unlist(lapply(adjustedData, function(x) { sd(x$Y, na.rm = TRUE) }))
    upspTable["p_value"] <- unlist(lapply(adjustedData, calPValue))
    upspTable <- upspTable[is.finite(upspTable$p_value),]
    upspTable["BH_FDR"] <- p.adjust(upspTable$p_value, method = "BH")

    return(upspTable)
}
