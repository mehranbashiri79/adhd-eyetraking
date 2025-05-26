# ---------------------------------------------------------------------
#  xplotmcmc_all.R          – refined April 2025
#  Draws density panels for every participant in vpList and
#  saves each figure as Fig_VP####.pdf in mcmc_plots.
# ---------------------------------------------------------------------
rm(list = ls())

library(coda)
library(ggplot2)
library(grid)

# ---------- participant IDs (vp = 1 removed) ----------
vpList <- c(
  1142, 1171, 1212, 1379, 1433, 1475, 1515, 1577, 1748, 1842, 1843, 1869,
  2120, 2211, 2226, 2255, 2259, 2460, 2478, 2575, 2614, 2705, 2944, 2964,
  3038, 3061, 3067, 3212, 3239, 3257, 3323, 3387, 3466, 3491, 3721, 3784,
  3898, 4145, 4291, 4316, 4428, 4482, 4550, 4566, 4620, 4731, 4775, 4802,
  4816, 4896, 4897, 4962, 4975, 5075, 5261, 5420, 5443, 5496, 5656, 6001,
   6880, 7280, 7413, 7705, 7797,
  7799, 7836, 8042, 8088, 8183, 8271, 8287, 8316, 8351, 8444, 8503, 8567,
  8610, 8680, 8831, 8878, 9082, 9113, 9148, 9509, 9526, 9569, 9687, 9691,
  9707, 9810, 9881, 9909
)

# ---------- where to save the PDFs ----------
plotDir <- "D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data/mcmc/mcmc_plots"
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

# ---------- constants ----------
nchains <- 3
nvar    <- 5
varlist <- c("nu", "r", "mt", "iota", "beta")

# ---------- loop over every participant ----------
for (vp in vpList) {
  
  inname <- sprintf("mcmc/mcmc_VP%d.rda", vp)
  if (!file.exists(inname)) {
    warning(sprintf("Chain file for VP %d not found – skipped.", vp))
    next
  }
  load(inname)                     # loads object 'mcmc' (list of 5 chains)
  
  outname <- file.path(plotDir, sprintf("Fig_VP%d.pdf", vp))
  pdf(outname, paper = "a4r", width = 12, height = 6)
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3)))
  
  # ---- draw one density panel per parameter ----
  for (var in seq_len(nvar)) {
    
    varname <- varlist[var]
    alldf   <- data.frame()
    
    for (i in seq_len(nchains)) {
      h      <- unlist(mcmc[[i]])
      chain  <- matrix(h, ncol = nvar + 3, byrow = FALSE)
      xdf    <- data.frame(x = chain[, var])
      alldf  <- rbind(alldf, xdf)
      
      if (i == 1)
        p <- ggplot(xdf) + geom_density(aes(x), colour = "grey")
      else
        p <- p + geom_density(data = xdf, aes(x), colour = "grey")
    }
    
    p <- p +
      geom_density(data = alldf, aes(x), colour = "black") +
      labs(x = varname, y = "density")
    
    row <- round(var / 3 - 0.5) + 1
    col <- var - 3 * (row - 1)
    if (row == 3) { row <- 2; col <- 3 }
    
    print(p, vp = viewport(layout.pos.row = row, layout.pos.col = col))
  }
  
  dev.off()
  message(sprintf("Saved: %s", outname))
}

message("Finished plotting all participants.")

