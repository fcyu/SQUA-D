# Plot volcano plot for Pretable3 and Pretable4

rm(list = ls())

library(cowplot)
library(ggplot2)
library(xlsx)


# parameters
upspPath <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\Pretable3.xlsx"
outputPlotPath <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\Pretable3.volcano.pdf"
fdrT <- 0.1
bw <- 0.3
verfiedUpsp <- "AT1G07660.1-K6-K9-K13;AT1G07820.1-K6-K9-K13;AT2G28740.1-K6-K9-K13;AT3G45930.1-K6-K9-K13;AT3G53730.1-K6-K9-K13;AT3G46320.1-K6-K9-K13;AT5G59690.1-K6-K9-K13;AT5G59970.1-K6-K9-K13"

# read plot
upspTable <- read.xlsx2(upspPath, 1, header = TRUE, colClasses = c(rep("character", 4), rep("numeric", 9), rep("integer", 2), rep("numeric", 4)))

axisLineSize <- 1
curveLineSize <- 1

# plot pretable4
logRatioT <- 0.5 * sd(upspTable$Adjusted_log_ratio_mean, na.rm = TRUE)
cat(sprintf("Log-ratio threshold: +-%f\n", logRatioT))
temp <- fdrT - upspTable$BH_FDR
pvalueT <- max(upspTable[fdrT - upspTable$BH_FDR >= 0,]$p_value)

upspTable$threshold <- 0
upspTable[(upspTable$BH_FDR <= fdrT) & (upspTable$Adjusted_log_ratio_mean >= logRatioT),]$threshold <- 1
upspTable[(upspTable$BH_FDR <= fdrT) & (upspTable$Adjusted_log_ratio_mean <= -1 * logRatioT),]$threshold <- 2
upspTable[grepl(verfiedUpsp, upspTable$UPSP),]$threshold <- 3
upspTable$threshold <- as.factor(upspTable$threshold)

xmax <- ceiling(max(abs(upspTable$Adjusted_log_ratio_mean)) * 1.01)
xmin <- -1 * xmax
pretable4G1 <- ggplot() +
  geom_point(data = subset(upspTable, threshold == 0), aes(x = Adjusted_log_ratio_mean, y = -10 * log10(p_value)), colour = "gray", size = 3, na.rm = TRUE) +
  geom_point(data = subset(upspTable, threshold == 1 | threshold == 2), aes(x = Adjusted_log_ratio_mean, y = -10 * log10(p_value), fill = Adjusted_log_ratio_mean, size = -10 * log10(p_value)), shape = 21, na.rm = TRUE) +
  scale_fill_gradient2(low = "green", mid = "yellow", high = "red", midpoint = -1 * logRatioT) +
  scale_size(range = c(3, 7)) +
  geom_point(data = subset(upspTable, threshold == 3), aes(x = Adjusted_log_ratio_mean, y = -10 * log10(p_value), size = -10 * log10(p_value)), fill = "skyblue", shape = 21, na.rm = TRUE) +
  theme(legend.position = "none", axis.title.x = element_blank(), axis.line.x = element_line(color = "black", size = axisLineSize), axis.line.y = element_line(color = "black", size = axisLineSize), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +
  ylab(expression(paste(-10, Log[10], "(", italic(p), "-value)"))) +
  xlim(c(xmin, xmax)) +
  # scale_color_manual(values = c("gray", "red", "blue", "green")) +
  geom_hline(yintercept = -10 * log10(pvalueT), linetype = "dashed", color = "black") +
  geom_vline(xintercept = logRatioT, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -1 * logRatioT, linetype = "dashed", color = "black")

pretable4G2 <- ggplot(data = upspTable, aes(x = Adjusted_log_ratio_mean)) +
  geom_histogram(binwidth = bw, fill = "green", colour = "blue", alpha = 0.5, na.rm = TRUE) +
  theme(legend.position = "none", axis.line.x = element_line(color = "black", size = axisLineSize), axis.line.y = element_line(color = "black", size = axisLineSize), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) +
  xlim(c(xmin, xmax)) +
  xlab(expression(paste(Log[2], "(treatment / control)"))) +
  ylab("# UPSP") +
  stat_function(fun = function(x, mean, sd, n, bw) { n * bw * dnorm(x = x, mean = mean, sd = sd) }, args = list(mean = mean(upspTable$Adjusted_log_ratio_mean), sd = sd(upspTable$Adjusted_log_ratio_mean), n = nrow(upspTable), bw = bw), colour = "red", size = curveLineSize) +
  geom_vline(xintercept = logRatioT, linetype = "dashed", color = "black") +
  geom_vline(xintercept = -1 * logRatioT, linetype = "dashed", color = "black")

g <- plot_grid(pretable4G1, pretable4G2, align = "v", nrow = 2)
ggsave(outputPlotPath, g, width = 8, height = 5)

cat("Done!\n")