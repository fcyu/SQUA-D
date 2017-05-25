# Plot spectral count V.S. UPSP index

# TODO: change color manuall, finish g2 and g3

library(ggplot2)
library(cowplot)

rm(list = ls())

# parameters
table1Path <- "D:\\Results\\Mascot\\Acetylation\\summary.UPS.FR.tsv"
table2Path <- "D:\\Results\\Mascot\\Acetylation\\summary.UPS.FR.tsv.adjusted.tsv"
fdrT <- 0.15
probT <- 0.8

# read tables
table1 <- read.table(table1Path, header = TRUE, sep = "\t", quote = "", comment.char = "", fill = TRUE)
table2 <- read.table(table2Path, header = TRUE, sep = "\t", quote = "", comment.char = "", fill = TRUE)

upspCountDF <- data.frame(upsp = table1$UPSP, lightCount = table1$IDs..Light., heavyCount = table1$IDs..Heavy., totalCount = table1$IDs..Light. + table1$IDs..Heavy., sinProb = table1$SIn.Probability, spectralCountProb = table1$Spectral.Count.Probability)
upspCountDF$fdr <- 1

for (upsp in table2$UPSP) {
    fdr <- table2[table2$UPSP == upsp,]$BH.FDR
    if (fdr <= fdrT) {
        temp <- upsp == upspCountDF$upsp
        if (sum(temp) == 1) {
            upspCountDF[temp,]$fdr <- fdr
        } else if (sum(temp) > 1) {
            stop(sprintf("Something is wrong: %s\n.", upsp))
        }
    }
}

upspCountDF <- upspCountDF[order(upspCountDF$totalCount),]
upspCountDF$index <- seq(1, nrow(upspCountDF))
upspCountDF$Source <- "None"

temp <- upspCountDF$fdr <= fdrT & upspCountDF$sinProb < probT & upspCountDF$spectralCountProb < probT
if (any(temp)) {
    upspCountDF[temp,]$Source <- "XIC"
}

temp <- upspCountDF$fdr > fdrT & upspCountDF$sinProb >= probT & upspCountDF$spectralCountProb < probT
if (any(temp)) {
    upspCountDF[temp,]$Source <- "SIn"
}

temp <- upspCountDF$fdr > fdrT & upspCountDF$sinProb < probT & upspCountDF$spectralCountProb >= probT
if (any(temp)) {
    upspCountDF[temp,]$Source <- "Rs"
}

temp <- upspCountDF$fdr > fdrT & upspCountDF$sinProb >= probT & upspCountDF$spectralCountProb >= probT
if (any(temp)) {
    upspCountDF[temp,]$Source <- "Rs and SIn"
}

temp <- upspCountDF$fdr <= fdrT & upspCountDF$sinProb >= probT & upspCountDF$spectralCountProb < probT
if (any(temp)) {
    upspCountDF[temp,]$Source <- "XIC and SIn"
}

temp <- upspCountDF$fdr <= fdrT & upspCountDF$sinProb < probT & upspCountDF$spectralCountProb >= probT
if (any(temp)) {
    upspCountDF[temp,]$Source <- "XIC and Rs"
}

temp <- upspCountDF$fdr <= fdrT & upspCountDF$sinProb >= probT & upspCountDF$spectralCountProb >= probT
if (any(temp)) {
    upspCountDF[temp,]$Source <- "All"
}

upspCountDF$Source <- factor(upspCountDF$Source, levels = c("None", "XIC", "SIn", "Rs", "Rs and SIn", "XIC and SIn", "XIC and Rs", "All"))

# plot parameters
myAlpha <- 1
pointSize <- 1

g1 <- ggplot(data = subset(upspCountDF, Source != "None"), aes(order = Source)) +
  geom_point(aes(x = index, y = lightCount, color = Source), alpha = myAlpha, size = pointSize, na.rm = TRUE) +
  scale_color_manual(values = c("XIC" = "red", "SIn" = "blue", "Rs" = "green", "Rs and SIn" = "cyan", "XIC and SIn" = "magenta", "XIC and Rs" = "yellow", "All" = "black")) +
  xlab("UPSP Index") +
  ylab("PSM Count") +
  theme(legend.position = c(.1, .75))

g2 <- ggplot(data = subset(upspCountDF, Source != "None"), aes(order = Source)) +
  geom_point(aes(x = index, y = heavyCount, color = Source), alpha = myAlpha, size = pointSize) +
  scale_color_manual(values = c("XIC" = "red", "SIn" = "blue", "Rs" = "green", "Rs and SIn" = "cyan", "XIC and SIn" = "magenta", "XIC and Rs" = "yellow", "All" = "black")) +
  xlab("UPSP Index") +
  ylab("PSM Count") +
  theme(legend.position = c(.1, .75))

g3 <- ggplot(data = subset(upspCountDF, Source != "None"), aes(order = Source)) +
  geom_point(aes(x = index, y = totalCount, color = Source), alpha = myAlpha, size = pointSize) +
  scale_color_manual(values = c("XIC" = "red", "SIn" = "blue", "Rs" = "green", "Rs and SIn" = "cyan", "XIC and SIn" = "magenta", "XIC and Rs" = "yellow", "All" = "black")) +
  xlab("UPSP Index") +
  ylab("PSM Count") +
  theme(legend.position = c(.1, .75))

p <- plot_grid(g1, g2, g3, align = "v", nrow = 3, labels = c("A", "B", "C"))
title <- ggdraw() + draw_label(sprintf("BH FDR <= %.2f, Prob >= %.2f", fdrT, probT), fontface = "bold")
plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
