# --------------------------------------------------------------------
#  xppcheck_all.R     – Posterior-predictive checks for all readers
#  Robust version  •  April 2025
# --------------------------------------------------------------------
rm(list = ls())

# ---------------- 1. Libraries & helper functions -------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(dplyr)
  library(tidyr)
})

source("functions/probskip.R")
source("functions/probrefix.R")
source("functions/probregress.R")
source("functions/singlefixdur.R")
source("functions/gazedur.R")
source("functions/totalreadtime.R")

# ---------------- 2. Participant list (vp ≠ 1) ----------------------
vpList <- c(
  1142, 1171, 1212, 1379, 1433, 1475, 1515, 1577, 1748, 1842, 1843,
  1869, 2120, 2211, 2226, 2255, 2259, 2460, 2478, 2575, 2614, 2705,
  2944, 2964, 3038, 3061, 3067, 3212, 3239, 3257, 3323, 3387, 3466,
  3491, 3721, 3784, 3898, 4145, 4291, 4316, 4428, 4482, 4550, 4566,
  4620, 4731, 4775, 4802, 4816, 4896, 4897, 4962, 4975, 5075, 5261,
  5420, 5443, 5496, 5656, 6001, 6154, 6314, 6423, 6692, 6723, 6741,
  6743, 6880, 7280, 7413, 7705, 7797, 7799, 7836, 8042, 8088, 8183,
  8271, 8287, 8316, 8351, 8444, 8503, 8567, 8610, 8680, 8831, 8878,
  9082, 9113, 9148, 9509, 9526, 9569, 9687, 9691, 9707, 9810, 9881,
  9909
)

# How many of the *last* sentences to use:
last_n_sentences <- 19

# ---------------- 3. File paths -------------------------------------
baseDir  <- "D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data"
corpusFP <- file.path(baseDir, "expdata/corpus.dat")
expFP    <- file.path(baseDir, "expdata/%d.dat")
simFP    <- file.path(baseDir, "ppsim/ppsim_VP%d.dat")

plotDir  <- file.path(baseDir, "ppsim/sim_plots")
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

# ---------------- 4. Load & normalize corpus ------------------------
corpus <- read.table(corpusFP, header = TRUE, stringsAsFactors = FALSE)
names(corpus) <- c("sentID", "nw", "wordID", "length", "freq", "unknown")

corpus$length <- as.numeric(corpus$length)
corpus$freq  <- 1
corpus$lfreq <- 1
# lookup: sentID → number of words
mIDNW <- setNames(corpus$nw, corpus$sentID)

# ---------------- 5. Storage matrix ---------------------------------
metrics <- c(
  "fixdur_exp","fixdur_sim","pskip_exp","pskip_sim",
  "prefix_exp","prefix_sim","pregress_exp","pregress_sim",
  "gazedur_exp","gazedur_sim","totaltime_exp","totaltime_sim"
)
res <- matrix(
  NA_real_, 
  nrow = length(vpList), 
  ncol = length(metrics),
  dimnames = list(vpList, metrics)
)

# ---------------- helper to clean data frames -----------------------
clean_df <- function(df) {
  df %>%
    mutate(across(everything(), as.numeric)) %>%
    drop_na(sentID, word, fixdur)
}

# ---------------- 6. Robust per-participant loop --------------------
message("Running PPC for ", length(vpList), " participants ...")
good_vp <- c()     # participants processed successfully

for (vp in vpList) {
  
  expFile <- sprintf(expFP, vp)
  simFile <- sprintf(simFP, vp)
  if (!file.exists(expFile) || !file.exists(simFile)) {
    warning(" >> Missing data for VP ", vp)
    next
  }
  
  # read & clean
  exp <- clean_df(read.table(expFile, header = TRUE, stringsAsFactors = FALSE))
  sim <- clean_df(
    read.table(simFile, header = TRUE, stringsAsFactors = FALSE) %>%
      select(sentID, word, fixdur)
  )
  
  # pick only their *last* N sentences (by order of appearance in exp)
  last_sents <- tail(unique(exp$sentID), last_n_sentences)
  exp <- filter(exp, sentID %in% last_sents)
  sim <- filter(sim, sentID %in% last_sents)
  
  # compute all metrics on those subsets
  ok <- try({
    res[as.character(vp), ] <- c(
      mean(singlefixdur(exp$sentID, exp$word, exp$fixdur, mIDNW), na.rm = TRUE),
      mean(singlefixdur(sim$sentID, sim$word, sim$fixdur, mIDNW), na.rm = TRUE),
      probskip(exp$sentID, exp$word,      mIDNW),
      probskip(sim$sentID, sim$word,      mIDNW),
      probrefix(exp$sentID, exp$word,     mIDNW),
      probrefix(sim$sentID, sim$word,     mIDNW),
      probregress(exp$sentID, exp$word,   mIDNW),
      probregress(sim$sentID, sim$word,   mIDNW),
      mean(gazedur(exp$sentID, exp$word, exp$fixdur, mIDNW), na.rm = TRUE),
      mean(gazedur(sim$sentID, sim$word, sim$fixdur, mIDNW), na.rm = TRUE),
      mean(totalreadtime(exp$sentID, exp$word, exp$fixdur, mIDNW), na.rm = TRUE),
      mean(totalreadtime(sim$sentID, sim$word, sim$fixdur, mIDNW), na.rm = TRUE)
    )
  }, silent = TRUE)
  
  if (inherits(ok, "try-error")) {
    warning(" >> Metric failure for VP ", vp, "; skipping.")
    next
  }
  good_vp <- c(good_vp, vp)
  
  # individual scatter plots
  indivPDF <- file.path(plotDir, sprintf("PPC_VP%d.pdf", vp))
  pdf(indivPDF, width = 6, height = 7, paper = "a4r")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3,2)))
  specs <- list(
    list("fixdur",   200,450, "Mean Fixation Duration"),
    list("pskip",      0,  1, "Skipping Probability"),
    list("gazedur",  200,450, "Mean Gaze Duration"),
    list("prefix",     0,.25,"Refixation Probability"),
    list("totaltime",200,450, "Total Reading Time"),
    list("pregress",   0,.25,"Regression Probability")
  )
  gpair <- function(met) {
    data.frame(
      exp = res[as.character(vp), paste0(met,"_exp")],
      sim = res[as.character(vp), paste0(met,"_sim")]
    )
  }
  for (k in seq_along(specs)) {
    ps <- specs[[k]]
    df <- gpair(ps[[1]])
    p  <- ggplot(df, aes(exp, sim)) +
      geom_point(size = 3) +
      geom_abline(linetype = "dashed", colour = "grey50") +
      coord_fixed(xlim = c(ps[[2]], ps[[3]]),
                  ylim = c(ps[[2]], ps[[3]])) +
      labs(title = ps[[4]], x="experimental", y="simulated") +
      theme_bw(base_size = 9)
    print(p, vp = viewport(
      layout.pos.row = ceiling(k/2),
      layout.pos.col = ifelse(k%%2==1,1,2)
    ))
  }
  dev.off()
  message("   ✔ VP ", vp)
}

# ---------------- 7. Combined plot ----------------------------
if (length(good_vp) > 1) {
  dfAll <- as.data.frame(res[as.character(good_vp), ], row.names = NULL) %>%
    mutate(vp = good_vp) %>%
    pivot_longer(-vp, names_to = c("metric","source"), names_sep="_") %>%
    pivot_wider(names_from = source, values_from = value) %>%
    drop_na()
  
  combPDF <- file.path(plotDir, "PPC_ALL_participants.pdf")
  pdf(combPDF, width = 6, height = 7, paper = "a4r")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3,2)))
  specs <- list(
    list("fixdur",   200,450, "Mean Fixation Duration"),
    list("pskip",      0,  1, "Skipping Probability"),
    list("gazedur",  200,450, "Mean Gaze Duration"),
    list("prefix",     0,.25,"Refixation Probability"),
    list("totaltime",200,450, "Total Reading Time"),
    list("pregress",   0,.25,"Regression Probability")
  )
  for (k in seq_along(specs)) {
    ps  <- specs[[k]]
    sub <- dfAll %>% filter(metric == ps[[1]])
    p   <- ggplot(sub, aes(exp, sim)) +
      geom_point(size = 2) +
      geom_abline(linetype = "dashed", colour = "grey50") +
      coord_fixed(xlim = c(ps[[2]], ps[[3]]),
                  ylim = c(ps[[2]], ps[[3]])) +
      labs(title = ps[[4]], x="experimental", y="simulated") +
      theme_bw(base_size = 9)
    print(p, vp = viewport(
      layout.pos.row = ceiling(k/2),
      layout.pos.col = ifelse(k%%2==1,1,2)
    ))
  }
  dev.off()
  message("\nCombined plot saved to ", combPDF)
} else {
  warning("No successful participants for combined plot.")
}

# ---------------- 8. Save numeric summary ---------------------------
write.csv(res, file.path(plotDir, "..", "ppcheck_results.csv"), row.names = FALSE)
saveRDS(res, file.path(plotDir, "..", "ppcheck_results.rds"))
message("Numeric table written to ppcheck_results.*\nDone.")
