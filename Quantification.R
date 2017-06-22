# Quantification

rm(list = ls())

options(java.parameters = "-Xmx32g")

library(stringr)
library(mzR)
library(xlsx)
library(Peptides)
library(Biostrings)
library(plyr)
library(parallel)
library(NOISeq)
library(stringr)
library(readr)


source("QuantificationFun.R")
source("AdjustLogRatioBatchEffect.R")
source("XIC.R")


# parameters. Will be put to command arguement later on
experimentStructurePath <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\experiment_structure.tsv"
mascotResultDir <- "D:\\Dropbox\\Results\\Mascot\\Acetylation\\"
spectraDir <- "D:\\Li_group\\Acetylation_raw_data\\"
fastaPath <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\TAIR10_pep_20101214_updated_with_AT3G53420.3.fasta"
outputDir <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\"
allPsmTablePath <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\all_psm_table.tsv"
# allPsmTablePath <- ""
nsaf <- FALSE
sin <- FALSE
ppmTol <- 20
rtTol <- 10
threadNum <- detectCores()
psmReplicateT <- 2
labellingReplicateT <- 1
experimentReplicateT <- 5
frT <- 0.15
performBatchEffectAdjustment <- TRUE
minCharge <- 2
maxCharge <- 5
deltaScoreT <- 10

# constants
modTable <- data.frame(mod = character(), mass = character(), stringsAsFactors = FALSE)
modTable <- rbind(modTable, data.frame(mod = "\\(Dimethyl\\)", mass = "(28.03)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Dimethyl:2H\\(4\\)13C\\(2\\)\\)", mass = "(34.06)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Acetyl\\)", mass = "(42.01)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Oxidation\\)", mass = "(15.99)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Carbamidomethyl\\)", mass = "(57.02)", stringsAsFactors = FALSE))
lightLabel <- "\\(28.03\\)"
heavyLabel <- "\\(34.06\\)"
labelDiff <- 6.03
kingMod <- "(42.01)"

# read experiment type
experimentStructure <- read.table(experimentStructurePath, sep = "\t", quote = "", header = TRUE, comment.char = "", colClasses = c("character", "factor", "factor"))

psmTableFinal <- NULL
if (str_length(allPsmTablePath) > 0) {
  cat(str_c(Sys.time(), ": Already have the PSM table. Reading it...\n", sep = ""))
  psmTableFinal <- read.table(allPsmTablePath, sep = "\t", quote = "", header = TRUE, comment.char = "", na.strings = c("NA", ""), colClasses = c("character", rep("factor", 2), "integer", rep("numeric", 4), "integer", rep("character", 2), "integer", rep("character", 2), "factor", rep("numeric", 4), "character", "factor", rep("numeric", 4), "integer", rep("numeric", 6)))
} else {
  cat(str_c(Sys.time(), ": Analyzing each Mascot result...\n", sep = ""))
  # read protein database
  db <- readAAStringSet(fastaPath)
  # read Mascot results, calculate spectral counting meassure, and generate a PSM table
  psmTableList <- lapply(experimentStructure$run_name, analyzeEachRun, mascotResultDir, spectraDir, experimentStructure, db, ppmTol, rtTol, minCharge, maxCharge, threadNum, modTable, lightLabel, heavyLabel, kingMod)
  psmTableFinal <- ldply(psmTableList, data.frame)
  write.table(psmTableFinal, str_c(outputDir, "all_psm_table.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)
}

cl <- makeCluster(threadNum)
clusterExport(cl, c("str_detect", "str_c", "str_match_all", "str_locate", "str_split", "psmReplicateT"))

cat(str_c(Sys.time(), ": Filtering PSMs based on delta score threshold ", deltaScoreT, "...\n", sep = ""))
psmTableFinal <- psmTableFinal[psmTableFinal$delta_score >= deltaScoreT,]

cat(str_c(Sys.time(), ": Generating Pretable1...\n", sep = ""))
allPeptide <- unique(psmTableFinal$label_free_peptide)
peptideTableList <- parLapply(cl, seq(1, length(allPeptide)), generatePeptideTable, allPeptide, psmTableFinal)
peptideTable <- ldply(peptideTableList, data.frame)
write.xlsx2(peptideTable, str_c(outputDir, "Pretable1.xlsx"), showNA = FALSE, row.names = FALSE, col.names = TRUE)

# summarize number of king modification site in each protein from all PSMs
cat(str_c(Sys.time(), ": Summarizing protein-modification count from all PSMs...\n", sep = ""))
allUpsp <- unlist(str_split(unlist(psmTableFinal$UPSP), ";"))
allProId <- unique(unlist(str_extract(allUpsp, "[^-]+")))
proteinModTable <- data.frame(proteinId = allProId, site = NA, stringsAsFactors = FALSE)
for (upsp in allUpsp) {
  tempVector <- str_extract_all(upsp, "[^-]+")[[1]]
  if (is.na(proteinModTable[proteinModTable$proteinId == tempVector[1], ]$site)) {
    proteinModTable[proteinModTable$proteinId == tempVector[1],]$site <- list(tempVector[seq(2, length(tempVector))])
  } else {
    proteinModTable[proteinModTable$proteinId == tempVector[1],]$site <- list(unique(c(unlist(proteinModTable[proteinModTable$proteinId == tempVector[1],]$site), tempVector[seq(2, length(tempVector))])))
  }
}
proteinModNumTable <- data.frame(proteinId = allProId, count = unlist(lapply(proteinModTable$site, length)))
write.table(proteinModNumTable, str_c(outputDir, "all_protein_mod_count.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

cat(str_c(Sys.time(), ": Filtering PSMs with more than ", psmReplicateT, " replicates...\n", sep = ""))
allPeptide <- unique(psmTableFinal$original_peptide)
filteredIdx <- unlist(parLapply(cl, allPeptide, filterBasedOnPsmReplication, psmTableFinal))
psmTableFinal <- psmTableFinal[filteredIdx,]

write.table(psmTableFinal, str_c(outputDir, "psm_table.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

# generate Pretable2
cat(str_c(Sys.time(), ": Generating Pretable2...\n", sep = ""))
allPeptide <- unique(psmTableFinal$label_free_peptide)
peptideTableList <- parLapply(cl, seq(1, length(allPeptide)), generatePeptideTable, allPeptide, psmTableFinal)
peptideTable <- ldply(peptideTableList, data.frame)
write.xlsx2(peptideTable, str_c(outputDir, "Pretable2.xlsx"), showNA = FALSE, row.names = FALSE, col.names = TRUE)

# summarize number of king modification site in each protein
cat(str_c(Sys.time(), ": Summarizing protein-modification count from filtered PSMs...\n", sep = ""))
allUpsp <- unlist(str_split(unlist(psmTableFinal$UPSP), ";"))
allProId <- unique(unlist(str_extract(allUpsp, "[^-]+")))
proteinModTable <- data.frame(proteinId = allProId, site = NA, stringsAsFactors = FALSE)
for (upsp in allUpsp) {
  tempVector <- str_extract_all(upsp, "[^-]+")[[1]]
  if (is.na(proteinModTable[proteinModTable$proteinId == tempVector[1], ]$site)) {
    proteinModTable[proteinModTable$proteinId == tempVector[1],]$site <- list(tempVector[seq(2, length(tempVector))])
  } else {
    proteinModTable[proteinModTable$proteinId == tempVector[1],]$site <- list(unique(c(unlist(proteinModTable[proteinModTable$proteinId == tempVector[1],]$site), tempVector[seq(2, length(tempVector))])))
  }
}
proteinModNumTable <- data.frame(proteinId = allProId, count = unlist(lapply(proteinModTable$site, length)))
write.table(proteinModNumTable, str_c(outputDir, "protein_mod_count.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

# filtering based on spectral counting criteria
cat(str_c(Sys.time(), ": Filtering PSMs based on spectral counting criteria...\n", sep = ""))
allPeptide <- unique(psmTableFinal$label_free_peptide)
filteredIdx <- unlist(parLapply(cl, allPeptide, filterBasedOnSpectralCountCriteria, psmTableFinal, labellingReplicateT, experimentReplicateT, frT))
psmTableFinal <- psmTableFinal[filteredIdx,]

# normalized log-ratio
psmTableFinal["normalized_log_ratio"] <- NA
allBatchType <- unique(experimentStructure$batch_type)
for (batchType in allBatchType) {
  logRatioMedian <- median(psmTableFinal[psmTableFinal$batch_type == batchType,]$log_ratio, na.rm = TRUE)
  psmTableFinal[psmTableFinal$batch_type == batchType,]$normalized_log_ratio <- psmTableFinal[psmTableFinal$batch_type == batchType,]$log_ratio - logRatioMedian
}
write.table(psmTableFinal, str_c(outputDir, "filtered_psm_table.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

# convert PSM table to UPSP table
cat(str_c(Sys.time(), ": Converting PSMs to UPSPs...\n", sep = ""))
allUpsp <- unique(psmTableFinal$UPSP)
upspTableList <- parLapply(cl, seq(1, length(allUpsp)), generateUpspTable, allUpsp, psmTableFinal)
upspTableFinal <- ldply(upspTableList, data.frame)

if (nsaf) {
  # NSAF
  cat(str_c(Sys.time(), ": Calculating NSAF...\n", sep = ""))
  allExpType <- sort(unique(experimentStructure$experiment_type))
  nsafMatrix <- matrix(unlist(parLapply(cl, upspTableFinal$UPSP, calculateNsaf, psmTableFinal, experimentStructure)), nrow = nrow(upspTableFinal), byrow = TRUE)
  rownames(nsafMatrix) <- unlist(upspTableFinal$UPSP)
  colnames(nsafMatrix) <- rep(NA, ncol(nsafMatrix))
  colnames(nsafMatrix)[seq(1, ncol(nsafMatrix), 2)] <- str_c(allExpType, "-C", sep = "")
  colnames(nsafMatrix)[seq(2, ncol(nsafMatrix), 2)] <- str_c(allExpType, "-T", sep = "")
  write.table(nsafMatrix, str_c(outputDir, "NSAF.tsv"), quote = FALSE, sep = "\t", na = "NA", row.names = TRUE, col.names = NA)
  
  #NOISeq
  noiseqFactors <- data.frame(experiment = rep(c("C", "T"), ncol(nsafMatrix) / 2))
  nsafNoiseqData <- readData(data = nsafMatrix, factors = noiseqFactors)
  nsafNoiseq <- noiseqbio(nsafNoiseqData, k = 0.5 * min(nsafMatrix, na.rm = TRUE), norm = "n", factor = "experiment", filter = 0, plot = TRUE)
  temp <- log2(nsafMatrix[, seq(2, ncol(nsafMatrix), 2)] / nsafMatrix[, seq(1, ncol(nsafMatrix), 2)])
  upspTableFinal["NSAF_log_ratio"] <- apply(temp, 1, mean, na.rm = TRUE)
  upspTableFinal["NSAF_NOISeq_Prob"] <- nsafNoiseq@results[[1]]$prob
}

if (sin) {
  # SIn
  cat(str_c(Sys.time(), ": Calculating SIn\n", sep = ""))
  sinMatrix <- matrix(unlist(parLapply(cl, upspTableFinal$UPSP, calculateSin, psmTableFinal, experimentStructure)), nrow = nrow(upspTableFinal), byrow = TRUE)
  rownames(sinMatrix) <- unlist(upspTableFinal$UPSP)
  colnames(sinMatrix) <- rep(NA, ncol(sinMatrix))
  colnames(sinMatrix)[seq(1, ncol(sinMatrix), 2)] <- str_c(allExpType, "-C", sep = "")
  colnames(sinMatrix)[seq(2, ncol(sinMatrix), 2)] <- str_c(allExpType, "-T", sep = "")
  write.table(sinMatrix, str_c(outputDir, "SIn"), quote = FALSE, sep = "\t", na = "NA", row.names = TRUE, col.names = NA)
  
  # NOISeq
  noiseqFactors <- data.frame(experiment = rep(c("C", "T"), ncol(nsafMatrix) / 2))
  sinNoiseqData <- readData(data = sinMatrix, factors = noiseqFactors)
  sinNoiseq <- noiseqbio(sinNoiseqData, k = 0.5 * min(sinMatrix, na.rm = TRUE), norm = "n", factor = "experiment", filter = 0, plot = TRUE)
  temp <- log2(sinMatrix[, seq(2, ncol(sinMatrix), 2)] / sinMatrix[, seq(1, ncol(sinMatrix), 2)])
  upspTableFinal["SIn_log_ratio"] <- apply(temp, 1, mean, na.rm = TRUE)
  upspTableFinal["SIn_NOISeq_Prob"] <- sinNoiseq@results[[1]]$prob
}

stopCluster(cl)

# Adjust batch effect, t-test, and BH FDR
if (performBatchEffectAdjustment) {
  cat(str_c(Sys.time(), ": Adjusting batch effects, estimating p-values and FDRs...\n", sep = ""))
  upspTableFinal <- adjustBatchEffect(upspTableFinal, psmTableFinal, debug = FALSE)
} else {
  cat(str_c(Sys.time(), ": Estimating p-values and FDRs...\n", sep = ""))
  upspTableFinal["p_value"] <- unlist(lapply(seq(1, nrow(upspTableFinal)), function(i, upspTableFinal, psmTableFinal) { t.test(psmTableFinal[(upspTableFinal[i,]$UPSP == psmTableFinal$UPSP) & !is.na(psmTableFinal$normalized_log_ratio),]$normalized_log_ratio)$p.value }, upspTableFinal, psmTableFinal))
  upspTableFinal <- upspTableFinal[is.finite(upspTableFinal$p_value),]
  upspTableFinal["BH_FDR"] <- p.adjust(upspTableFinal$p_value, method = "BH")
}
write.xlsx2(upspTableFinal, str_c(outputDir, "Pretable3.xlsx"), showNA = FALSE, row.names = FALSE, col.names = TRUE)

closeAllConnections()
cat(str_c(Sys.time(), ": Done!\n"), sep = "")