# Merge F and R MS1 pairs for each UPS

library(jpeg)

rm(list = ls())

# parameters
plotsDir <- "D:\\Dropbox\\Results\\Quantification\\Acetylation\\Pretable3.xlsx.plot\\selected_plots\\"

mergedPlotsDir <- paste(plotsDir, "merged\\", sep = "")
if (!dir.exists(mergedPlotsDir)) {
    dir.create(mergedPlotsDir)
}

plotsNames <- sort(list.files(plotsDir, pattern = ".+\\.jpg$", full.names = FALSE, no.. = TRUE))

for (i in seq(1, floor(length(plotsNames) / 2))) {
    temp <- strsplit(plotsNames[i * 2 - 1], "-", fixed = TRUE)[[1]]
    ups1 <- paste(temp[1:(length(temp) - 3)], collapse = "-")
    temp <- strsplit(plotsNames[i * 2], "-", fixed = TRUE)[[1]]
    ups2 <- paste(temp[1:(length(temp) - 3)], collapse = "-")
    if (grepl(ups1, ups2)) {
        fig1 <- readJPEG(paste(plotsDir, plotsNames[i * 2 - 1], sep = ""), native = TRUE)
        fig2 <- readJPEG(paste(plotsDir, plotsNames[i * 2], sep = ""), native = TRUE)
        if (grepl("R", plotsNames[i * 2 - 1], ignore.case = FALSE)) {
            temp <- fig1
            fig1 <- fig2
            fig2 <- temp
        }

        jpeg(paste(mergedPlotsDir, ups1, ".jpg", sep = ""), width = 25, height = 15, units = "cm", res = 300)
        par(mfrow = c(1, 2), yaxs = "i", xaxs = "i", cex.main = 3, mai = c(0, 0, 0, 0), oma = c(0, 0, 2, 0))
        plot(NA, xlim = c(0, 5), ylim = c(0, 6), xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
        rasterImage(fig1, 0, 0, 5, 6)
        plot(NA, xlim = c(0, 5), ylim = c(0, 6), xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
        rasterImage(fig2, 0, 0, 5, 6)
        title(ups1, outer = TRUE)
        dev.off()
    } else {
        cat(sprintf("%s and %s have different UPS.\n", plotsNames[i * 2 - 1], plotsNames[i * 2]))
        stop()
    }
}

cat("Done!")