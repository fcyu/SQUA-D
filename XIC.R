# extract ion chromatogram

library(xcms)
library(caTools)
library(stats)
library(parallel)
library(stringr)
library(plyr)

## A function to find the turning point
# Based on MaxQuant criteria.
findLiftRightLocalMinima <- function(x, startIdx, correctedStartIdx) {
  startIdx1 <- startIdx
  startIdx2 <- correctedStartIdx
  if (startIdx > correctedStartIdx) {
    startIdx1 <- correctedStartIdx
    startIdx2 <- startIdx
  }
  
  # get the 1st order derivation
  xDev <- x[seq(2, length(x))] - x[seq(1, length(x) - 1)]
  xDev[is.na(xDev)] <- 0
  
  rightLocalMinimumIdx <- NA
  leftLocalMinimumIdx <- NA
  
  # find right turning point
  if (startIdx2 > length(xDev) - 4) {
    rightLocalMinimumIdx <- length(x)
  } else {
    for (i in seq(startIdx2 + 1, length(xDev) - 1)) {
      if ((xDev[i] <= 0) && (xDev[i + 1] >= 0)) {
        # this is a local minimum
        if (i < length(xDev) - 1) {
          for (j in seq(i + 1, length(xDev) - 1)) {
            if ((xDev[j] >= 0) && (xDev[j + 1] <= 0)) {
              # this is a local maxima
              if (x[i] <= (1 / 1.3) * min(x[startIdx2], x[j])) {
                rightLocalMinimumIdx <- i + 1
                break
              }
            }
          }
        } else {
          rightLocalMinimumIdx <- length(x)
        }
        if (!is.na(rightLocalMinimumIdx)) {
          break
        }
      }
    }
    if (is.na(rightLocalMinimumIdx)) {
      rightLocalMinimumIdx <- length(x)
    }
  }
  
  
  # find left turning point
  if (startIdx1 < 4) {
    leftLocalMinimumIdx <- 1
  } else {
    for (i in seq(startIdx1 - 1, 2, -1)) {
      if ((xDev[i - 1] <= 0) && (xDev[i] >= 0)) {
        # this is a local minimum
        if (i > 2) {
          for (j in seq(i - 1, 2, -1)) {
            if ((xDev[j - 1] >= 0) && (xDev[j] <= 0)) {
              # this is a local maxima
              if (x[i] <= (1 / 1.3) * min(x[startIdx1], x[j])) {
                leftLocalMinimumIdx <- i - 1
                break
              }
            }
          }
        } else {
          leftLocalMinimumIdx <- 1
        }
        if (!is.na(leftLocalMinimumIdx)) {
          break
        }
      }
    }
    if (is.na(leftLocalMinimumIdx)) {
      leftLocalMinimumIdx <- 1
    }
  }
  
  return(c(leftLocalMinimumIdx, rightLocalMinimumIdx))
}


## A function to locate RT, RT range, and calculate the maximum intensity
getRtIntensityInfo <- function(xcmsRawObj, mz, rt, ppmTol, rtTol) {
  leftMz <- min(max(mz * (1 - ppmTol * 1e-6), xcmsRawObj@mzrange[1]), xcmsRawObj@mzrange[2])
  rightMz <- max(min(mz * (1 + ppmTol * 1e-6), xcmsRawObj@mzrange[2]), xcmsRawObj@mzrange[1])
  rtLeft <- rt - rtTol
  rtRight <- rt + rtTol
  
  # get XIC
  temp <- NULL
  eicObj <- getEIC(xcmsRawObj, matrix(c(leftMz, rightMz), nrow = 1, ncol = 2))
  temp <- eicObj@eic$xcmsRaw[[1]]
  
  rtArray <- temp[, 1]
  intensityArray <- temp[, 2]
  
  # moving average smooth
  smoothedIntensityArray <- filter(x = intensityArray, filter = rep(1 / 3, 3), method = "convolution", sides = 2, circular = FALSE)
  smoothedIntensityArray[is.na(smoothedIntensityArray)] <- 0
  
  rtIdx <- which.min(abs(rtArray - rep(rt, length(rtArray))))
  
  # using local maxima to locate the correct RT
  rtLeftIdx <- which.min(abs(rtArray - rep(rtLeft, length(rtArray))))
  rtRightIdx <- which.min(abs(rtArray - rep(rtRight, length(rtArray))))
  intensityLocalMaxima <- max(smoothedIntensityArray[rtLeftIdx:rtRightIdx])
  correctedRtIdx <- which(smoothedIntensityArray[rtLeftIdx:rtRightIdx] == intensityLocalMaxima)[1] + rtLeftIdx - 1;
  # caution: There may be multiple point having intensities equals intensityLocalMaxima, we only care about the one in the [rtLeftIdx, rtRightIdx] range.
  correctedRt <- rtArray[correctedRtIdx]
  
  # find left and right end turning points
  leftRightTurningPoints <- findLiftRightLocalMinima(smoothedIntensityArray, rtIdx, correctedRtIdx)
  leftRt <- rtArray[leftRightTurningPoints[1]]
  rightRt <- rtArray[leftRightTurningPoints[2]]
  
  # get MS1 intensities of the RT range
  # tempX <- rtArray[leftRightTurningPoints[1] : leftRightTurningPoints[2]]
  # tempY <- smoothedIntensityArray[leftRightTurningPoints[1] : leftRightTurningPoints[2]]
  # intensity <- trapz(tempX, abs(tempY)) # caution: the smoothing may cause negative value
  
  # return the local mixima as the intensity as what MaxQuant does
  return(c(correctedRt, leftRt, rightRt, intensityLocalMaxima))
}


extractRtIntensity <- function(i, xcmsRawObj, psmTable, ppmTol, rtTol) {
  mz <- psmTable[i,]$spectrum_mz
  rt <- psmTable[i,]$rt
  charge <- psmTable[i,]$charge
  res <- getRtIntensityInfo(xcmsRawObj, mz, rt, ppmTol, rtTol)
  psmTable[i,]$rt_corrected <- res[1]
  psmTable[i,]$rt_left <- res[2]
  psmTable[i,]$rt_right <- res[3]
  psmTable[i,]$intensity <- res[4]
  return(psmTable[i,])
}


extractAnotherRtIntensity <- function(i, xcmsRawObj, psmTable, ppmTol, rtTol, minCharge, maxCharge) {
  if (is.na(psmTable[i, ]$another_scan_num)) {
    if ((psmTable[i,]$another_mz > xcmsRawObj@mzrange[1]) && (psmTable[i,]$another_mz < xcmsRawObj@mzrange[2])) {
      anotherRes <- getRtIntensityInfo(xcmsRawObj, psmTable[i,]$another_mz, psmTable[i,]$rt, ppmTol, rtTol)
      # check charge state
      charges <- inferCharge(xcmsRawObj, anotherRes[1], psmTable[i,]$another_mz, ppmTol, minCharge, maxCharge)
      if ((length(charges) > 0) && (psmTable[i,]$charge %in% charges)) {
        psmTable[i,]$another_rt <- anotherRes[1]
        psmTable[i,]$another_rt_left <- anotherRes[2]
        psmTable[i,]$another_rt_right <- anotherRes[3]
        psmTable[i,]$another_intensity <- anotherRes[4]
      }
    }
  }
  return(psmTable[i,])
}


inferCharge <- function(xcmsRawObj, rt, mz, ppmTol, minCharge, maxCharge) {
  # This method is very simple and naive. But we cannot do much about centroid data.
  scanIdx <- which.min(abs(xcmsRawObj@scantime - rt)) # the index is not scan num
  specMatrix <- getScan(xcmsRawObj, scanIdx, c(max(mz - 0.5, xcmsRawObj@mzrange[1]), min(mz + 3, xcmsRawObj@mzrange[2])))
  charges <- c()
  if (nrow(specMatrix) > 1) {
    mzIdx <- which.min(abs(specMatrix[, 1] - mz))
    if (mzIdx < nrow(specMatrix)) {
      daTol <- mz * ppmTol * 1e-6
      for (charge in seq(minCharge, maxCharge)) {
        mzDiff <- 1.00335483 / charge
        if (min(abs(specMatrix[seq(mzIdx + 1, nrow(specMatrix)), 1] - mz - mzDiff)) <= daTol) {
          charges <- unique(c(charges, charge))
        }
      }
    }
  }
  return(charges)
}


xic <- function(spectraPath, psmTable, ppmTol, rtTol, minCharge, maxCharge, threadNum) {
  # add new columns to the table to store results
  psmTable["rt_corrected"] <- NA
  psmTable["rt_left"] <- NA
  psmTable["rt_right"] <- NA
  psmTable["intensity"] <- NA
  psmTable["another_scan_num"] <- NA
  psmTable["another_rt"] <- NA
  psmTable["another_rt_left"] <- NA
  psmTable["another_rt_right"] <- NA
  psmTable["another_intensity"] <- NA
  
  xcmsRawObj <- xcmsRaw(spectraPath, profstep = 0)
  
  cl <- makeCluster(threadNum)
  clusterExport(cl, c("getRtIntensityInfo", "findLiftRightLocalMinima", "getEIC", "filter", "str_detect", "inferCharge", "getScan"))
  
  # extract rt and intensity information
  psmRowList <- parLapply(cl, seq(1, nrow(psmTable)), extractRtIntensity, xcmsRawObj, psmTable, ppmTol, rtTol)
  # psmRowList <- lapply(seq(1, nrow(psmTable)), extractRtIntensity, xcmsRawObj, psmTable, ppmTol, rtTol)
  psmTable <- ldply(psmRowList, data.frame)
  
  # pair light and heavy PSMs. Here, one heavy PSM is allowed to pair multiple light PSM; for a light PSM, only the closed and overlapped heavy one is paired.
  for (i in which(psmTable$label_type == "L")) {
    tempIdx <- which((psmTable$label_type == "H") & (psmTable[i,]$label_free_peptide == psmTable$label_free_peptide))
    if (length(tempIdx) > 0) {
      j <- tempIdx[which.min(abs(psmTable[i,]$rt - psmTable[tempIdx,]$rt))]
      if (abs(psmTable[i,]$rt_left < psmTable[j, ]$rt_right) || abs(psmTable[i,]$rt_right < psmTable[j,]$rt_left)) {
        psmTable[i,]$another_scan_num <- psmTable[j,]$scan_num
        psmTable[i,]$another_rt <- psmTable[j,]$rt_corrected
        psmTable[i,]$another_rt_left <- psmTable[j,]$rt_left
        psmTable[i,]$another_rt_right <- psmTable[j,]$rt_right
        psmTable[i,]$another_intensity <- psmTable[j,]$intensity
        psmTable[j,]$another_scan_num <- psmTable[i,]$scan_num
        psmTable[j,]$another_rt <- psmTable[i,]$rt_corrected
        psmTable[j,]$another_rt_left <- psmTable[i,]$rt_left
        psmTable[j,]$another_rt_right <- psmTable[i,]$rt_right
        psmTable[j,]$another_intensity <- psmTable[i,]$intensity
      }
    }
  }
  
  # extract another rt and intensity information for unpaired PSMs
  psmRowList <- parLapply(cl, seq(1, nrow(psmTable)), extractAnotherRtIntensity, xcmsRawObj, psmTable, ppmTol, rtTol, minCharge, maxCharge)
  # psmRowList <- lapply(seq(1, nrow(psmTable)), extractAnotherRtIntensity, xcmsRawObj, psmTable, ppmTol, rtTol, minCharge, maxCharge)
  psmTable <- ldply(psmRowList, data.frame)
  
  stopCluster(cl)
  
  return(psmTable)
}
