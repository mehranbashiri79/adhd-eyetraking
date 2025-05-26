#---------------------------------------------------------------------
#  xtoyDREAM (Multi-Participant Version for File Paths, with Updated Frequency & Fixation Processing)
#  (Version 2.2, April 2025) - Modified for iota tuning: adaptation increased from 0.2 to 0.4,
#                           and initial values revised.
#---------------------------------------------------------------------
library(BayesianTools)
library(Rcpp)
library(parallel)
library(dplyr)
library(tidyr)
library(MASS)
library(Rfast)
library(rstudioapi)

# Clear environment
rm(list = ls())

# Compile the C++ likelihood function
# (make sure toyLL.cpp is in your working directory)
sourceCpp("toyLL.cpp")  # Warning about Rcpp attributes can usually be ignored if the code runs correctly

# Optionally, set a seed for reproducibility
# set.seed(1234)

# Define a vector of participant IDs (modify according to your data)
vpList <- c(1142, 1171, 1212, 1379, 1433, 1475)
# Example: replace with your actual participant IDs

#---------------------------------------------------------------------
# Process Corpus Data (common for all participants)
#---------------------------------------------------------------------
corpus_path <- "D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data/expdata/corpus.dat"
corpus <- read.table(corpus_path, header = TRUE, stringsAsFactors = FALSE)
names(corpus) <- c("sentID", "nw", "wordID", "length", "freq", "unknown")
corpus$sentID <- as.numeric(corpus$sentID)
corpus$length <- as.numeric(corpus$length)
# 1. Calculate frequency using length^(-1.6) with scaling
corpus$freq <- 1 + (corpus$length^(-1.6)) / min(corpus$length^(-1.6))
# 2. Find maximum frequency (REQUIRED for normalization)
fmax <- max(corpus$freq)
# 3. Log-normalize to [0,1] range
corpus$lfreq <- log10(corpus$freq) / log10(fmax)


cat("----- Summary of Corpus Data -----\n")
print(summary(corpus))
cat("\nUnique Frequency Values after Transformation:\n")
print(unique(corpus$freq))
cat("\n----- Summary of Normalized Frequencies (lfreq) -----\n")
print(summary(corpus$lfreq))

#---------------------------------------------------------------------
# Define the Log-Likelihood Function (same for all participants)
#---------------------------------------------------------------------
LL <- function(corpus, data, nu, r, mt, iota, beta, gamma) {
  Nruns <- 1  # Experimental data: Nsubj = 1
  common_sentIDs <- intersect(unique(corpus$sentID), unique(data$sentID))
  if (length(common_sentIDs) < 38) {
    warning("Only ", length(common_sentIDs), " common sentences available. Using all available sentences.")
    Nsent <- length(common_sentIDs)
  } else {
    Nsent <- 38
  }
  SIDvec <- common_sentIDs[1:Nsent]
  loglik <- 0
  for (run in 1:Nruns) {
    for (sent in SIDvec) {
      idx <- which(corpus$sentID == sent)
      lfreq <- corpus$lfreq[idx]
      idx <- which(data$run == run & data$sentID == sent)
      fseq <- cbind(data$word[idx], data$fixdur[idx])
      fixpos <- fseq[, 1]
      fixdur <- fseq[, 2]
      loglik <- loglik + toyLL_cpp(lfreq, fixpos, fixdur, nu, r, mt, iota, beta)
    }
  }
  LL <- loglik / (Nruns * Nsent)
  return(LL)
}

#---------------------------------------------------------------------
# Combine Priors Function and Setup (same for all participants)
#---------------------------------------------------------------------
combinePriors <- function(...) {
  priors <- list(...)
  npar <- vapply(priors, function(x) length(x$lower), integer(1))
  lastpar <- cumsum(npar)
  firstpar <- lag(lastpar, default = 0) + 1
  createPrior(
    density = function(x) do.call(sum, lapply(seq_along(priors), function(i) 
      priors[[i]]$density(x[firstpar[i]:lastpar[i]]))),
    sampler = function() {
      result <- double(sum(npar))
      for (i in seq_along(priors)) {
        result[firstpar[i]:lastpar[i]] <- priors[[i]]$sampler()
      }
      result
    },
    lower = do.call(c, lapply(priors, function(x) x$lower)),
    upper = do.call(c, lapply(priors, function(x) x$upper)),
    best  = do.call(c, lapply(priors, function(x) x$best))
  )
}
prior <- combinePriors(
  createUniformPrior(0, 1),      # nu
  createUniformPrior(0, 12),     # r
  createUniformPrior(100, 400),  # mt
  createUniformPrior(-0.5, 2),   # iota
  createUniformPrior(0, 1)       # beta
)

# Set the initial "best" values so that the sampler uses these values as starting points.
# Here, we propose: nu = 0.45, r = 9.7, mt = 230, iota = 0.1, beta = 0.9.
bayesianSetup <- createBayesianSetup(
  likelihood = function(pars) do.call(LL, c(list(corpus), list(data), as.list(pars))),
  prior = prior,
  names = c("nu", "r", "mt", "iota", "beta"),
  #best = c(nu = 0.45, r = 9.7, mt = 230, iota = 0.1, beta = 0.9),
  parallel = FALSE,
  catchDuplicates = FALSE
)

#---------------------------------------------------------------------
# MCMC Settings and Sampling (same for all participants)
#---------------------------------------------------------------------
Nchains <- 3
Niter <- 20000
setlist <- list(iterations = Niter * Nchains,
                nCR = Nchains,
                pCRupdate = FALSE,
                pSnooker = 0.1, 
                adaptation = 0.2,
                consoleUpdates = 100,
                startValue = Nchains,  # This will use the bayesianSetup$best values
                parallel = FALSE)

# Load the custom DREAMzs sampler and assign it
source("D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data/dreamzs2.R")
assignInNamespace("DREAMzs", DREAMzs2, ns = "BayesianTools")

#---------------------------------------------------------------------
# Loop over Participants
#---------------------------------------------------------------------
for (vp in vpList) {
  cat("\nProcessing participant ID:", vp, "\n")
  
  # Load fixation data for the current participant
  fixation_path <- sprintf("D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data/expdata/%d.dat", vp)
  data <- read.table(fixation_path, header = TRUE, stringsAsFactors = FALSE)
  names(data) <- c("sentID", "word", "fixpos", "fixdur")
  data$run <- 1
  data$sentID <- as.numeric(data$sentID)
  data$word   <- as.numeric(data$word)
  data$fixdur <- as.numeric(data$fixdur)
  data$fixpos <- as.numeric(data$fixpos)
  
  # Update the bayesianSetup for the current participant's fixation data.
  # (The corpus remains the same.)
  bayesianSetup <- createBayesianSetup(
    likelihood = function(pars) do.call(LL, c(list(corpus), list(data), as.list(pars))),
    prior = prior,
    names = c("nu", "r", "mt", "iota", "beta"),
    #best = c(nu = 0.45, r = 9.7, mt = 230, iota = 0.1, beta = 0.9),
    parallel = FALSE,
    catchDuplicates = FALSE
  )
  
  # Run MCMC sampler for current participant
  mcmc1 <- runMCMC(bayesianSetup, "DREAMzs", settings = setlist)
  cat("Summary for participant", vp, ":\n")
  print(summary(mcmc1))
  
  # Save the MCMC chain to file
  outname <- sprintf("D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data/mcmc/mcmc_VP%d.rda", vp)
  mcmc <- mcmc1$chain
  save(mcmc, file = outname)
}

