# Convert Percolator output file to TSV file

rm(list = ls())

library(stringr)
library(readr)

# constants
modTable <- data.frame(mod = character(), mass = character(), stringsAsFactors = FALSE)
modTable <- rbind(modTable, data.frame(mod = "\\(Dimethyl\\)", mass = "(28.03)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Dimethyl:2H\\(4\\)13C\\(2\\)\\)", mass = "(34.06)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Acetyl\\)", mass = "(42.01)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Oxidation\\)", mass = "(15.99)", stringsAsFactors = FALSE))
modTable <- rbind(modTable, data.frame(mod = "\\(Carbamidomethyl\\)", mass = "(57.02)", stringsAsFactors = FALSE))

formatPeptide <- function(peptide, modTable) {
  peptide <- str_replace(peptide, "^[A-Z]\\.", "n")
  peptide <- str_replace(peptide, "\\.[A-Z]", "")
  for (i in seq(1, nrow(modTable))) {
    peptide <- str_replace_all(peptide, modTable[i,]$mod, modTable[i,]$mass)
  }
  temp <- str_extract_all(peptide, "[A-Zn](\\(([0-9\\.\\-]+)\\))?", simplify = TRUE)
  modStr <- ""
  for (i in seq(1, length(temp))) {
    if (str_detect(temp[i], "\\(")) {
      massStr <- str_sub(temp[i], 3, str_length(temp[i]) - 1)
      if (i == 1) {
        modStr <- str_c(modStr, massStr, "@N-term;", sep = "")
      } else {
        modStr <- str_c(modStr, massStr, "@", (i - 1), ";", sep = "")
      }
    }
  }
  return(c(str_replace_all(peptide, "[\\(\\)0-9\\-\\.nc]+", ""), modStr))
}

## main part
mascotResultDir <- "D:\\Dropbox\\Results\\Mascot\\Acetylation\\"
outputPath <- "D:\\Dropbox\\Results\\Mascot\\Acetylation\\ProteinProspector.tsv"

fileNames <- list.files(mascotResultDir, pattern = ".+\\.percolator\\.psms$")
outputTable <- data.frame(scanNum = numeric(), charge = integer(), fraction = character(), peptide = character())
for (fileName in fileNames) {
  mascotTable <- read_tsv(str_c(mascotResultDir, fileName, sep = ""), col_types = "cdddc", quote = "", comment = "")
  mascotTable <- mascotTable[mascotTable$`q-value` <= 0.01,]
  scanInfoMatrix <- str_match(mascotTable$PSMId, "scan=([0-9]+)\";rt:([0-9.]+);mz:([0-9.]+);charge:([0-9])")
  scanNums <- as.integer(scanInfoMatrix[, 2])
  charges <- as.integer(scanInfoMatrix[, 5])
  tempList <- lapply(mascotTable$peptide, formatPeptide, modTable)
  tempMatrix <- matrix(unlist(tempList), nrow = length(tempList), byrow = TRUE)
  tempTable <- data.frame(scanNum = scanNums, charge = charges, fraction = str_replace(fileName, "\\.percolator\\.psms", ""), peptide = tempMatrix[, 1], mod = tempMatrix[, 2])
  outputTable <- rbind(outputTable, tempTable)
}
write.table(outputTable, outputPath, quote = FALSE, sep = "\t", row.names = FALSE)
