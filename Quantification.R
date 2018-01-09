# Quantification

rm(list = ls())


suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(mzR)))
suppressMessages(suppressWarnings(library(Peptides)))
suppressMessages(suppressWarnings(library(Biostrings)))
suppressMessages(suppressWarnings(library(plyr)))
suppressMessages(suppressWarnings(library(parallel)))
suppressMessages(suppressWarnings(library(NOISeq)))
suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(readr)))


source("QuantificationFun.R")
source("AdjustLogRatioBatchEffect.R")
source("XIC.R")

threadNum <- detectCores()

commandLine <- FALSE

# parameters. Only for internal use.
mascotResultDir <- "D:\\Dropbox\\Results\\Mascot\\Acetylation\\"
spectraDir <- "D:\\Li_group\\Acetylation_raw_data\\"
outputDir <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\"
experimentStructurePath <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\experiment_structure.tsv"
fastaPath <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\TAIR10_pep_20101214_updated_with_AT3G53420.3.fasta"
allPsmTablePath <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\Table_S1a.tsv"
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

sdFactor <- 0.5
bhFdrT <- 0.1

if (commandLine) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) != 15) {
    for (a in args) {
      cat(str_c(a, "\n"))
    }
    stop("Incorrect args.")
  }
  
  mascotResultDir <- args[2]
  spectraDir <- args[2]
  outputDir <- args[2]
  experimentStructurePath <- args[3]
  fastaPath <- args[4]
  nsaf <- FALSE
  sin <- FALSE
  ppmTol <- as.integer(args[5])
  rtTol <- as.integer(args[6])
  psmReplicateT <- as.integer(args[7])
  labellingReplicateT <- as.integer(args[8])
  experimentReplicateT <- as.integer(args[9])
  frT <- as.numeric(args[10])
  performBatchEffectAdjustment <- args[11] == "1"
  minCharge <- as.integer(args[12])
  maxCharge <- as.integer(args[13])
  deltaScoreT <- as.numeric(args[14])
}

# print the arguments
cat("Parameters:\n")
cat(sprintf("mascotResultDir: %s\n", mascotResultDir))
cat(sprintf("spectraDir: %s\n", spectraDir))
cat(sprintf("outputDir: %s\n", outputDir))
cat(sprintf("experimentStructurePath: %s\n", experimentStructurePath))
cat(sprintf("fastaPath: %s\n", fastaPath))
cat(sprintf("nsaf: %d\n", nsaf))
cat(sprintf("sin: %d\n", sin))
cat(sprintf("ppmTol: %d\n", ppmTol))
cat(sprintf("rtTol: %d\n", rtTol))
cat(sprintf("psmReplicateT: %d\n", psmReplicateT))
cat(sprintf("labellingReplicateT: %d\n", labellingReplicateT))
cat(sprintf("experimentReplicateT: %d\n", experimentReplicateT))
cat(sprintf("frT: %f\n", frT))
cat(sprintf("performBatchEffectAdjustment: %d\n", performBatchEffectAdjustment))
cat(sprintf("minCharge: %d\n", minCharge))
cat(sprintf("maxCharge: %d\n", maxCharge))
cat(sprintf("deltaScoreT: %f\n", deltaScoreT))

# constants
modTable <- data.frame(mod = character(), mass = character(), stringsAsFactors = FALSE)
modTable <- rbind(modTable, data.frame(mod = "\\(Dimethyl\\)", mass = "(28.03)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Dimethyl:2H\\(4\\)13C\\(2\\)\\)", mass = "(34.06)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Acetyl\\)", mass = "(42.01)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Oxidation\\)", mass = "(15.99)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Carbamidomethyl\\)", mass = "(57.02)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Deamidated\\)", mass = "(0.98)", stringsAsFactors = FALSE))

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
  db <- readAAStringSet(fastaPath)
  psmTableList <- lapply(experimentStructure$run_name, analyzeEachRun, mascotResultDir, spectraDir, experimentStructure, db, ppmTol, rtTol, minCharge, maxCharge, threadNum, modTable, lightLabel, heavyLabel, kingMod)
  psmTableFinal <- ldply(psmTableList, data.frame)
  write.table(psmTableFinal, str_c(outputDir, "Table_S1a.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)
}

cl <- makeCluster(threadNum)
clusterExport(cl, c("str_detect", "str_c", "str_match_all", "str_locate", "str_split", "psmReplicateT"))

cat(str_c(Sys.time(), ": Filtering PSMs based on delta score threshold ", deltaScoreT, " (Table_S1b)...\n", sep = ""))
psmTableFinal <- psmTableFinal[psmTableFinal$delta_score >= deltaScoreT,]
write.table(psmTableFinal, str_c(outputDir, "Table_S1b.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

cat(str_c(Sys.time(), ": Summarizing protein-modification count from all PSMs (all_protein_mod_count)...\n", sep = ""))
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

cat(str_c(Sys.time(), ": Filtering PSMs with more than ", psmReplicateT, " replicates (Table_S1c)...\n", sep = ""))
allPeptide <- unique(psmTableFinal$original_peptide)
filteredIdx <- unlist(parLapply(cl, allPeptide, filterBasedOnPsmReplication, psmTableFinal, psmReplicateT))
psmTableFinal <- psmTableFinal[filteredIdx,]
write.table(psmTableFinal, str_c(outputDir, "Table_S1c.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

cat(str_c(Sys.time(), ": Generating Table_S1d...\n", sep = ""))
allPeptide <- unique(psmTableFinal$label_free_peptide)
peptideTableList <- parLapply(cl, seq(1, length(allPeptide)), generatePeptideTable, allPeptide, psmTableFinal)
peptideTable <- ldply(peptideTableList, data.frame)
write.table(peptideTable, str_c(outputDir, "Table_S1d.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

cat(str_c(Sys.time(), ": Summarizing protein-modification count from filtered PSMs (protein_mod_count)...\n", sep = ""))
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

cat(str_c(Sys.time(), ": Filtering PSMs based on spectral counting criteria (Table_S2a)...\n", sep = ""))
allPeptide <- unique(psmTableFinal$label_free_peptide)
filteredIdx <- unlist(parLapply(cl, allPeptide, filterBasedOnSpectralCountCriteria, psmTableFinal, labellingReplicateT, experimentReplicateT, frT))
psmTableFinal <- psmTableFinal[filteredIdx,]
allPeptide <- unique(psmTableFinal$label_free_peptide)
peptideTableList <- parLapply(cl, seq(1, length(allPeptide)), generatePeptideTable, allPeptide, psmTableFinal)
peptideTable <- ldply(peptideTableList, data.frame)
write.table(peptideTable, str_c(outputDir, "Table_S2a.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

cat(str_c(Sys.time(), ": Normalizing log-ratios...\n", sep = ""))
psmTableFinal["normalized_log_ratio"] <- NA
allBatchType <- unique(experimentStructure$batch_type)
for (batchType in allBatchType) {
  logRatioMedian <- median(psmTableFinal[psmTableFinal$batch_type == batchType,]$log_ratio, na.rm = TRUE)
  psmTableFinal[psmTableFinal$batch_type == batchType,]$normalized_log_ratio <- psmTableFinal[psmTableFinal$batch_type == batchType,]$log_ratio - logRatioMedian
}

cat(str_c(Sys.time(), ": Generating filtered_psm_table...\n", sep = ""))
write.table(psmTableFinal, str_c(outputDir, "filtered_psm_table.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

cat(str_c(Sys.time(), ": Converting peptides to UPSPs...\n", sep = ""))
allUpsp <- unique(psmTableFinal$UPSP)
upspTableList <- parLapply(cl, seq(1, length(allUpsp)), generateUpspTable, allUpsp, psmTableFinal)
upspTableFinal <- ldply(upspTableList, data.frame)

if (nsaf) {
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

cat(str_c(Sys.time(), ": Generating Table_S2b...\n", sep = ""))
write.table(upspTableFinal, str_c(outputDir, "Table_S2b.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

cat(str_c(Sys.time(), ": Generating Table_3a...\n", sep = ""))
if (performBatchEffectAdjustment) {
  logRatioT <- sdFactor * sd(upspTableFinal$Adjusted_log_ratio_mean)
  upspTableFinal <- upspTableFinal[upspTableFinal$BH_FDR <= bhFdrT & abs(upspTableFinal$Adjusted_log_ratio_mean) >= logRatioT,]
} else {
  logRatioT <- sdFactor * sd(upspTableFinal$Log_ratio_mean)
  upspTableFinal <- upspTableFinal[upspTableFinal$BH_FDR <= bhFdrT & abs(upspTableFinal$Log_ratio_mean) >= logRatioT,]
}
write.table(upspTableFinal, str_c(outputDir, "Table_S3a.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

cat(str_c(Sys.time(), ": Generating Table_3b...\n", sep = ""))
temp <- sum(upspTableFinal$Number_of_log_ratio_from_forward) + sum(upspTableFinal$Number_of_log_ratio_from_reciprocal)
if (performBatchEffectAdjustment) {
  upspTableFinal$x <- (upspTableFinal$Number_of_log_ratio_from_forward + upspTableFinal$Number_of_log_ratio_from_reciprocal) * (2 ^ upspTableFinal$Adjusted_log_ratio_mean) / temp
} else {
  upspTableFinal$x <- (upspTableFinal$Number_of_log_ratio_from_forward + upspTableFinal$Number_of_log_ratio_from_reciprocal) * (2 ^ upspTableFinal$Log_ratio_mean) / temp
}
write.table(upspTableFinal, str_c(outputDir, "Table_S3b.tsv"), quote = FALSE, sep = "\t", na = "", row.names = FALSE)

closeAllConnections()
cat(str_c(Sys.time(), ": Done!\n"), sep = "")