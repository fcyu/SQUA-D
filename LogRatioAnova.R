# ANOVA test with Log2-intensity ratio

rm(list = ls())

suppressMessages(suppressWarnings(library(stringr)))
suppressMessages(suppressWarnings(library(readr)))


unadjustedPath <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\filtered_psm_table.tsv"
adjustedPath <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\adjusted_table.tsv"

unadjustedTable <- read_tsv(unadjustedPath, col_types = "ccciddddiccicccddddccdddddddddddd", quote = "")

unadjustedTable <- unadjustedTable[is.finite(unadjustedTable$normalized_log_ratio), ]
allDF <- data.frame(logRatio = unadjustedTable$normalized_log_ratio, FR = str_sub(unadjustedTable$run_name, 9, 9), biologicalReplicates = as.integer(str_sub(unadjustedTable$run_name, 6, 6)), tmb = str_sub(unadjustedTable$run_name, 8, 8))

cat("ANOVA...")
unadjustAovRes <- aov(logRatio ~ FR * biologicalReplicates * tmb, allDF)
summary(unadjustAovRes)
boxplot(logRatio ~ FR * biologicalReplicates * tmb, allDF)
title("Boxplot of Unadjusted Data")

adjustedTable <- read_tsv(adjustedPath, col_types = "dcic", quote = "")
adjustedAovRes <- aov(logRatio ~ FR * biologicalReplicates * tmb, adjustedTable)
summary(adjustedAovRes)
boxplot(logRatio ~ FR * biologicalReplicates * tmb, adjustedTable)
title("Boxplot of Adjusted Data")