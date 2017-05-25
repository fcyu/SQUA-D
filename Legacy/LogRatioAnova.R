# ANOVA test with Log2-intensity ratio

library(stats)

rm(list = ls())

inputPath <- inputPath <- "D:\\Results\\Mascot\\Acetylation\\summary.UPS.FR.tsv"

organelleName <- "all"
organellePattern <- paste("Batch._", organelleName, "..", sep = "")

inputTable <- read.table(inputPath, header = TRUE, sep = "\t", fill = TRUE, comment.char = "", quote = "")

cat("Summarizing log-ratios...\n")
allDF <- data.frame(logRatio = numeric(), type = numeric(), batch = numeric(), organelle = character())
lastProgress <- 0
for (rowIdx in seq(1, nrow(inputTable))) {
  progress = floor(rowIdx * 100 / nrow(inputTable))
  if (progress > lastProgress) {
    cat(sprintf("\r%d%%", progress))
    lastProgress <- progress
  }
  
  myRow <- inputTable[rowIdx,]
  ratioTablePath <- paste(inputPath, ".ratios\\", strsplit(as.character(myRow$UPS), split = ";")[[1]][1], ".tsv", sep = "")
  if (file.exists(ratioTablePath)) {
    ratioTable <- read.table(ratioTablePath, header = TRUE, sep = "\t", comment.char = "", quote = "")
    
    if (!grepl("all", organelleName)) {
      ratioTable <- ratioTable[grepl(organellePattern, ratioTable$fraction),]
    }
    
    logRatio <- ratioTable$refinedLogRatio
    
    if (length(logRatio) > 0) {
      type <- as.factor(ifelse(grepl("Batch._.F.", ratioTable$fraction), "F", "R"))
      batch <- as.factor(ifelse(grepl("Batch1_...", ratioTable$fraction), 1, 2))
      
      if (grepl("all", organelleName)) {
        organelle <- as.factor(substr(ratioTable$fraction, 8, 8))
        localDF <- data.frame(logRatio, type, batch, organelle)
        allDF <- rbind(allDF, localDF)
      } else {
        localDF <- data.frame(logRatio, type, batch)
        allDF <- rbind(allDF, localDF)
      }
    }
  }
}

cat("\nANOVA...")
if (grepl("all", organelleName)) {
  allAovRes <- aov(logRatio ~ type * batch * organelle, allDF)
  summary(allAovRes)
  boxplot(logRatio ~ type : batch : organelle, allDF)
  title("Boxplot of Unadjusted Data")
} else {
  allAovRes <- aov(logRatio ~ type * batch, allDF)
  summary(allAovRes)
  boxplot(logRatio ~ type : batch, allDF)
  title("Boxplot of Unadjusted Data")
}
