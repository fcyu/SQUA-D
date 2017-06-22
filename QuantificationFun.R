# functions for Quantification

formatPeptide <- function(peptide, modTable) {
  peptide <- str_replace(peptide, "^[A-Z]\\.", "n")
  peptide <- str_replace(peptide, "\\.[A-Z]", "c")
  for (i in seq(1, nrow(modTable))) {
    peptide <- str_replace_all(peptide, modTable[i,]$mod, modTable[i,]$mass)
  }
  return(peptide)
}


calMass <- function(peptide) {
  temp <- str_sub(peptide, 2, str_length(peptide) - 1)
  temp2 <- str_match_all(temp, "\\(([0-9.-]+)\\)")[[1]]
  ptmFreePeptide <- str_replace_all(temp, "\\([0-9.-]+\\)", "")
  return(mw(ptmFreePeptide, monoisotopic = TRUE) + sum(as.numeric(temp2[, 2])))
}


getLabelInfo <- function(peptide, lightLabel, heavyLabel) {
  tempLight <- str_detect(peptide, lightLabel)
  tempHeavy <- str_detect(peptide, heavyLabel)
  if (!(tempLight && tempHeavy)) {
    if (tempLight) {
      labelFreePeptide <- str_replace_all(peptide, lightLabel, "")
      return(data.frame(labelType = "L", num = str_count(peptide, lightLabel), labelFreePeptide = labelFreePeptide, stringsAsFactors = FALSE))
    } else if (tempHeavy) {
      labelFreePeptide <- str_replace_all(peptide, heavyLabel, "")
      return(data.frame(labelType = "H", num = str_count(peptide, heavyLabel), labelFreePeptide = labelFreePeptide, stringsAsFactors = FALSE))
    }
  }
  return(data.frame(labelType = NA, num = NA, labelFreePeptide = NA, stringsAsFactors = FALSE))
}


getGroupType <- function(labelType, experimentType) {
  if (str_detect(experimentType, "f")) {
    if (labelType == "L") {
      return("C")
    } else {
      return("T")
    }
  } else {
    if (labelType == "L") {
      return("T")
    } else {
      return("C")
    }
  }
}


getAnotherMz <- function(i, labelType, mz, labelNum, charge, labelDiff) {
  if (labelType[i] == "L") {
    return(mz[i] + labelNum[i] * labelDiff / charge[i])
  } else {
    return(mz[i] - labelNum[i] * labelDiff / charge[i])
  }
}


getProteinInfo <- function(labelFreePeptide, db, kingMod) {
  temp <- unlist(vmatchPattern(str_replace_all(labelFreePeptide, "[0-9.()nc-]+", ""), db, max.mismatch = 0, min.mismatch = 0, with.indels = FALSE, fixed = TRUE))
  if (length(temp) > 0) {
    # get local ptm site and index
    temp2 <- str_match_all(labelFreePeptide, "([A-Znc])(\\([0-9.-]+\\))?")[[1]]
    ptmLocalIdx <- which(temp2[, 3] == kingMod) - 1
    ptmSite <- temp2[ptmLocalIdx + 1, 2]
    # generate all possible UPSP
    upsp <- c()
    for (i in seq(1, length(temp@start))) {
      ptmGlobalIdx <- ptmLocalIdx + temp@start[i] - 1
      upsp <- c(upsp, str_c(str_match(temp@NAMES[i], "^[^ ]+"), str_c(str_c(ptmSite, ptmGlobalIdx, sep = ""), collapse = "-"), sep = "-"))
    }
    # eliminate redunant UPSP with .2 .3...
    temp3 <- str_split(upsp, "\\.", simplify = TRUE)
    temp3 <- cbind(temp3, 0)
    for (i in seq(1, nrow(temp3))) {
      if (as.numeric(str_sub(temp3[i, 2], 1, 1)) == min(as.numeric(str_sub(temp3[temp3[i, 1] == temp3[, 1], 2], 1, 1)))) {
        temp3[i, 3] <- 1
      }
    }
    keptIdx <- which(temp3[, 3] == 1)
    return(c(str_c(upsp[keptIdx], collapse = ";"), str_c(temp@NAMES[keptIdx], collapse = "&&")))
  } else {
    return(c(NA, NA))
  }
}


calculateMs2Intensity <- function(scanNum, spectra) {
  peak <- peaks(spectra, scanNum)
  return(sum(peak[, 2]))
}


calculateDeltaScore <- function(idx, scanNums, peptides, mascotDatTable) {
  ptmFreePeptide = str_replace_all(peptides[idx], "[0-9\\.\\(\\)nc\\-]+", "")
  subTable <- mascotDatTable[mascotDatTable$scanNum == scanNums[idx] & mascotDatTable$ptmFreePeptide == ptmFreePeptide,]
  if (nrow(subTable) == 1) {
    return(unlist(subTable$score))
  } else if (nrow(subTable) > 1) {
    scores <- sort(unlist(subTable$score), decreasing = TRUE)
    return(scores[1] - scores[2])
  } else {
    cat(str_c(Sys.time(), ": Cannot find the dat results for ", scanNums[idx], " ", peptides[idx], "\n", sep = ""))
    return(1)
  }
}


analyzeEachRun <- function(runName, mascotResultDir, spectraDir, experimentStructure, db, ppmTol, rtTol, minCharge, maxCharge, threadNum, modTable, lightLabel, heavyLabel, kingMod) {
  cat(str_c(Sys.time(), ": Analyzing ", runName, "...\n", sep = ""))
  mascotTable <- read_tsv(str_c(mascotResultDir, runName, ".percolator.psms", sep = ""), col_types = "cdddc", quote = "", comment = "")
  mascotTable <- mascotTable[mascotTable$`q-value` <= 0.01,]
  
  mascotDatTable <- read_tsv(str_c(mascotResultDir, runName, "_target.dat.tsv", sep = ""), col_types = "ciccd", quote = "", comment = "")
  temp <- str_match(mascotDatTable$scanId, "scan=([0-9]+)\"")
  mascotDatTable$scanNum <- as.integer(temp[, 2])
  
  if (nrow(mascotTable) > 0) {
    scanInfoMatrix <- str_match(mascotTable$PSMId, "scan=([0-9]+)\";rt:([0-9.]+);mz:([0-9.]+);charge:([0-9])")
    scanNums <- as.integer(scanInfoMatrix[, 2])
    
    peptide <- unlist(lapply(mascotTable$peptide, formatPeptide, modTable))
    labelInfoTable <- lapply(peptide, getLabelInfo, lightLabel, heavyLabel)
    labelInfoMatrix <- matrix(unlist(labelInfoTable), nrow = length(labelInfoTable), byrow = TRUE)
    
    deltaScore <- unlist(lapply(seq(1, length(scanNums)), calculateDeltaScore, scanNums, peptide, mascotDatTable))
    
    psmTable <- data.frame(
      run_name = runName,
      batch_type = experimentStructure[experimentStructure$run_name == runName,]$batch_type,
      experiment_type = experimentStructure[experimentStructure$run_name == runName,]$experiment_type,
      scan_num = scanNums,
      spectrum_mz = as.numeric(scanInfoMatrix[, 4]),
      spectrum_mass = (as.numeric(scanInfoMatrix[, 4]) - 1.00727646688) * as.numeric(scanInfoMatrix[, 5]),
      theo_mass = unlist(lapply(peptide, calMass)),
      rt = as.numeric(scanInfoMatrix[, 3]),
      charge = as.integer(scanInfoMatrix[, 5]),
      ptm_free_peptide = str_replace_all(peptide, "[0-9.()nc-]+", ""),
      label_free_peptide = labelInfoMatrix[, 3],
      label_count = as.integer(labelInfoMatrix[, 2]),
      protein_anno = NA,
      UPSP = NA,
      label_type = labelInfoMatrix[, 1], # heavy or light labelling
      score = as.numeric(mascotTable$score),
      delta_score = deltaScore,
      q = as.numeric(mascotTable$`q-value`),
      another_mz = NA,
      original_peptide = peptide,
      group_type = NA, # control or treatment group
      stringsAsFactors = FALSE
    )
    
    psmTable <- psmTable[!is.na(psmTable$label_type),]
    psmTable <- psmTable[str_detect(psmTable$label_free_peptide, kingMod),]
    
    if (nrow(psmTable) > 0) {
      psmTable$group_type <- unlist(lapply(psmTable$label_type, getGroupType, experimentStructure[experimentStructure$run_name == runName,]$experiment_type))
      psmTable$another_mz <- unlist(lapply(seq(1, nrow(psmTable)), getAnotherMz, labelType = psmTable$label_type, mz = psmTable$spectrum_mz, labelNum = psmTable$label_count, charge = psmTable$charge, labelDiff = labelDiff))
      
      cl <- makeCluster(threadNum)
      clusterExport(cl, c("vmatchPattern", "str_replace_all", "str_c", "str_match_all", "str_match", "str_sort", "str_split", "str_sub"))
      
      cat(str_c(Sys.time(), ": - Getting protein information...\n"), sep = "")
      proteinInfoList <- parLapply(cl, psmTable$label_free_peptide, getProteinInfo, db, kingMod)
      
      proteinInfoMatrix <- matrix(unlist(proteinInfoList), nrow = length(proteinInfoList), byrow = TRUE)
      psmTable$UPSP <- proteinInfoMatrix[, 1]
      psmTable$protein_anno <- proteinInfoMatrix[, 2]
      psmTable <- psmTable[!is.na(psmTable$UPSP),]
      
      # calculate XIC
      cat(str_c(Sys.time(), ": - Extracting ion chromatography...\n"), sep = "")
      spectraPath <- str_c(spectraDir, runName, ".mzXML", sep = "")
      psmTable <- xic(spectraPath, psmTable, ppmTol, rtTol, minCharge, maxCharge, threadNum)
      
      # calculate log-ratio
      # After the calculation, if log_ratio == NA, it means that it is paired to a light one.
      # If only another_scan_num == NA, it means that it's another pair doesn't have the identification result.
      # If all "another_*" == NA, it means that it doesn't have another pair, we calculate the ratio based on min intensity.
      minIntensity <- min(psmTable$intensity)
      psmTable["log_ratio"] <- NA
      lightIdx <- (psmTable$label_type == "L")
      heavyIdx <- ((psmTable$label_type == "H") & is.na(psmTable$another_scan_num)) # paired heavy version has been condisered in the light version
      noNaIdx <- (!is.na(psmTable$another_intensity) & (psmTable$another_intensity > 0))
      naIdx <- (is.na(psmTable$another_intensity) | (psmTable$another_intensity <= 0))
      if (str_detect(experimentStructure[experimentStructure$run_name == runName, ]$experiment_type, "f")) {
        psmTable[lightIdx & noNaIdx,]$log_ratio <- log2(psmTable[lightIdx & noNaIdx,]$another_intensity / psmTable[lightIdx & noNaIdx,]$intensity)
        psmTable[lightIdx & naIdx,]$log_ratio <- log2(0.5 * minIntensity / psmTable[lightIdx & naIdx,]$intensity)
        psmTable[heavyIdx & noNaIdx,]$log_ratio <- log2(psmTable[heavyIdx & noNaIdx,]$intensity / psmTable[heavyIdx & noNaIdx,]$another_intensity)
        psmTable[heavyIdx & naIdx,]$log_ratio <- log2(psmTable[heavyIdx & naIdx,]$intensity / (0.5 * minIntensity))
      } else {
        psmTable[lightIdx & noNaIdx,]$log_ratio <- log2(psmTable[lightIdx & noNaIdx,]$intensity / psmTable[lightIdx & noNaIdx,]$another_intensity)
        psmTable[lightIdx & naIdx,]$log_ratio <- log2(psmTable[lightIdx & naIdx,]$intensity / (0.5 * minIntensity))
        psmTable[heavyIdx & noNaIdx,]$log_ratio <- log2(psmTable[heavyIdx & noNaIdx,]$another_intensity / psmTable[heavyIdx & noNaIdx,]$intensity)
        psmTable[heavyIdx & naIdx,]$log_ratio <- log2(0.5 * minIntensity / psmTable[heavyIdx & naIdx,]$intensity)
      }
      
      stopCluster(cl)
      
      # calculate MS2 peak intensity for SIn later on
      cat(str_c(Sys.time(), ": - Extracting MS2 peak intensity...\n"), sep = "")
      spectra <- openMSfile(spectraPath)
      psmTable$MS2_intensity <- unlist(lapply(psmTable$scan_num, calculateMs2Intensity, spectra))
      
      return(psmTable)
    } else {
      return(data.frame())
    }
  } else {
    return(data.frame())
  }
}


generatePeptideTable <- function(i, allPeptide, psmTable) { # TODO: check
  subTable <- psmTable[psmTable$label_free_peptide == allPeptide[i],]
  ptmFreePeptide <- subTable[1,]$ptm_free_peptide
  labelFreePeptide <- subTable[1,]$label_free_peptide
  upsp <- unlist(str_split(subTable[1,]$UPSP, ";"))
  temp <- str_match_all(upsp, "([^-]+)-([A-Zn])([0-9]+).*")
  temp <- matrix(unlist(temp), length(temp), byrow = TRUE)
  localLocation <- 1
  if (temp[1, 3] != "n") {
    localLocation <- str_locate(ptmFreePeptide, temp[1, 3])[1]
  }
  peptideLocation <- as.integer(temp[, 4]) - localLocation + 1
  lightCount <- sum(subTable$label_type == "L")
  ms2Light <- ""
  if (lightCount > 0) {
    ms2Light <- str_c(str_c(subTable[subTable$label_type == "L",]$run_name, subTable[subTable$label_type == "L",]$scan_num, sep = "-"), collapse = ";")
  }
  heavyCount <- sum(subTable$label_type == "H")
  ms2Heavy <- ""
  if (heavyCount > 0) {
    ms2Heavy <- str_c(str_c(subTable[subTable$label_type == "H",]$run_name, subTable[subTable$label_type == "H",]$scan_num, sep = "-"), collapse = ";")
  }
  peptideTable <- data.frame(
    Accession_No. = str_c(temp[, 2], collapse = ";"),
    Description = subTable[1,]$protein_anno,
    Peptide_Location = str_c(peptideLocation, collapse = ";"),
    Sequence = ptmFreePeptide,
    Modification_Peptide = labelFreePeptide,
    IDs_Light = lightCount,
    MS2_Light = ms2Light,
    IDs_Heavy = heavyCount,
    MS2_Heavy = ms2Heavy,
    total_IDs = nrow(subTable),
    replicates = length(unique(subTable$experiment_type))
  )
  return(peptideTable)
}


filterBasedOnPsmReplication <- function(x, psmTableFinal) {
  outputIdx <- which(psmTableFinal$original_peptide == x)
  if (length(outputIdx) >= psmReplicateT) {
    return(outputIdx)
  } else {
    return(c())
  }
}


filterBasedOnSpectralCountCriteria <- function(x, psmTableFinal, labellingReplicateT, experimentReplicateT, frT) {
  outputIdx <- which(psmTableFinal$label_free_peptide == x)
  fCount <- sum(str_detect(psmTableFinal[outputIdx,]$experiment_type, "f"))
  rCount <- sum(str_detect(psmTableFinal[outputIdx,]$experiment_type, "r"))
  if ((sum(psmTableFinal[outputIdx, ]$label_type == "L") >= labellingReplicateT) && (sum(psmTableFinal[outputIdx, ]$label_type == "H") >= labellingReplicateT) && (length(unique(psmTableFinal[outputIdx, ]$experiment_type)) >= experimentReplicateT) && (fCount / (fCount + rCount) >= frT) && (rCount / (fCount + rCount) >= frT)) {
    return(outputIdx)
  } else {
    return(c())
  }
}


generateUpspTable <- function(i, allUpsp, psmTableFinal) {
  psmTable <- psmTableFinal[psmTableFinal$UPSP == allUpsp[i],]
  upspTable <- data.frame(
    UPSP = allUpsp[i],
    Annotation = psmTable[1,]$protein_anno,
    Peptide = str_c(unique(psmTable$ptm_free_peptide), collapse = ";"),
    Modified_peptide = str_c(unique(psmTable$original_peptide), collapse = ";"),
    Mascot_delta_score = max(psmTable$delta_score),
    Forward_1 = mean(psmTable[psmTable$experiment_type == "f1",]$normalized_log_ratio, na.rm = TRUE),
    Forward_2 = mean(psmTable[psmTable$experiment_type == "f2",]$normalized_log_ratio, na.rm = TRUE),
    Forward_3 = mean(psmTable[psmTable$experiment_type == "f3",]$normalized_log_ratio, na.rm = TRUE),
    Reciprocal_1 = mean(psmTable[psmTable$experiment_type == "r1",]$normalized_log_ratio, na.rm = TRUE),
    Reciprocal_2 = mean(psmTable[psmTable$experiment_type == "r2",]$normalized_log_ratio, na.rm = TRUE),
    Reciprocal_3 = mean(psmTable[psmTable$experiment_type == "r3",]$normalized_log_ratio, na.rm = TRUE),
    Log_ratio_mean = mean(psmTable$normalized_log_ratio, na.rm = TRUE),
    Log_ratio_sd = sd(psmTable$normalized_log_ratio, na.rm = TRUE),
    Number_of_log_ratio_from_forward = sum(str_detect(psmTable$experiment_type, "f") & is.finite(psmTable$normalized_log_ratio)),
    Number_of_log_ratio_from_reciprocal = sum(str_detect(psmTable$experiment_type, "r") & is.finite(psmTable$normalized_log_ratio)),
    stringsAsFactors = FALSE
  )
}


calculateNsaf <- function(upsp, psmTable, experimentStructure) {
  allExpType <- sort(unique(experimentStructure$experiment_type))
  nsafVector <- rep(0, length(allExpType) * 2)
  for (i in seq(1, nrow(experimentStructure))) {
    tempIdx <- which(allExpType == experimentStructure[i,]$experiment_type)
    temp1 <- psmTable$run_name == experimentStructure[i,]$run_name
    temp2 <- psmTable$group_type == "C"
    temp3 <- psmTable$group_type == "T"
    temp4 <- psmTable$UPSP == upsp
    nsafVector[2 * tempIdx - 1] <- nsafVector[2 * tempIdx - 1] + sum(temp1 & temp2 & temp4) / sum(temp1 & temp2)
    nsafVector[2 * tempIdx] <- nsafVector[2 * tempIdx] + sum(temp1 & temp3 & temp4) / sum(temp1 & temp3)
  }
  nsafVector[nsafVector == 0] <- NA
  return(nsafVector)
}


calculateSin <- function(upsp, psmTable, experimentStructure) {
  allExpType <- sort(unique(experimentStructure$experiment_type))
  sinVector <- rep(0, length(allExpType) * 2)
  for (i in seq(1, nrow(experimentStructure))) {
    tempIdx <- which(allExpType == experimentStructure[i,]$experiment_type)
    temp1 <- psmTable$run_name == experimentStructure[i,]$run_name
    temp2 <- psmTable$group_type == "C"
    temp3 <- psmTable$group_type == "T"
    temp4 <- psmTable$UPSP == upsp
    sinVector[2 * tempIdx - 1] <- sinVector[2 * tempIdx - 1] + sum(psmTable[temp1 & temp2 & temp4,]$MS2_intensity) / sum(psmTable[temp1 & temp2,]$MS2_intensity)
    sinVector[2 * tempIdx] <- sinVector[2 * tempIdx] + sum(psmTable[temp1 & temp3 & temp4,]$MS2_intensity) / sum(psmTable[temp1 & temp3,]$MS2_intensity)
  }
  sinVector[sinVector == 0] <- NA
  return(sinVector)
}