rm(list = ls())
library(coda)
library(dplyr)
library(rstudioapi)
# setwd(dirname(getActiveDocumentContext()$path))
source("toygen.R")

# Configuration
nchains <- 3
nvar    <- 5
nstat   <- 5000
MODE    <- 2
run_id  <- 1
nHPD    <- 0.68
last_n_sentences <- 19

# Participants
vpList <- c(1142, 1171, 1212, 1379, 1433, 1475, 1515, 1577, 1748, 1842, 1843,
            1869, 2120, 2211, 2226, 2255, 2259, 2460, 2478, 2575, 2614, 2705,
            2944, 2964, 3038, 3061, 3067, 3212, 3239, 3257, 3323, 3387, 3466,
            3491, 3721, 3784, 3898, 4145, 4291, 4316, 4428, 4482, 4550, 4566,
            4620, 4731, 4775, 4802, 4816, 4896, 4897, 4962, 4975, 5075, 5261,
            5420, 5443, 5496, 5656, 6001, 6154, 6314, 6423, 6692, 6723, 6741,
            6743, 6880, 7280, 7413, 7705, 7797, 7799, 7836, 8042, 8088, 8183,
            8271, 8287, 8316, 8351, 8444, 8503, 8567, 8610, 8680, 8831, 8878,
            9082, 9113, 9148, 9509, 9526, 9569, 9687, 9691, 9707, 9810, 9881,
            9909)

# Corpus normalization with new frequency calculation
corpusPath <- "D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data/expdata/corpus.dat"
corpus <- read.table(corpusPath, header = TRUE, stringsAsFactors = FALSE)
names(corpus) <- c("sentID", "nw", "wordID", "length", "freq", "unknown")
corpus$length <- as.numeric(corpus$length)

corpus$freq  <- 1
corpus$lfreq <- 1

# Storage for parameters
parma <- data.frame(
  vp    = integer(),
  sentID= integer(),
  nu    = double(),
  r     = double(),
  mt    = double(),
  iota  = double(),
  beta  = double(),
  stringsAsFactors = FALSE
)

# Simulation loop
for (vp in vpList) {
  fixPath <- sprintf(
    "D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data/expdata/%d.dat",
    vp
  )
  if (!file.exists(fixPath)) {
    warning("Fixation data not found for participant ", vp); next
  }
  
  fixData <- read.table(fixPath, header = TRUE, stringsAsFactors = FALSE)
  names(fixData) <- c("sentID", "word", "fixpos", "fixdur")
  
  # determine which sentences they actually read, in order of appearance:
  assigned <- unique(fixData$sentID)
  # keep only the last N:
  if (length(assigned) > last_n_sentences) {
    assigned <- tail(assigned, last_n_sentences)
  }
  
  mcmcPath <- sprintf(
    "D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data/mcmc/mcmc_VP%d.rda",
    vp
  )
  if (!file.exists(mcmcPath)) {
    warning("MCMC file not found for participant ", vp); next
  }
  load(mcmcPath)  # loads list `mcmc`
  if (!is.list(mcmc) || length(mcmc) < nchains) {
    warning("Unexpected MCMC structure for participant ", vp); next
  }
  
  simAll <- NULL
  for (sent in assigned) {
    # pick one chain, pull last nstat samples
    chain_id <- sample.int(nchains, 1)
    mat <- tail(mcmc[[chain_id]], nstat)
    if (!is.matrix(mat) || ncol(mat) < nvar) {
      warning("Bad chain for VP ", vp, ", sentence ", sent); next
    }
    colnames(mat)[1:nvar] <- c("nu","r","mt","iota","beta")
    
    # HPD filter
    hpd  <- matrix(unlist(HPDinterval(as.mcmc(mat[,1:nvar]), prob = nHPD)), ncol=2)
    keep <- apply(mat[,1:nvar], 1, function(x) all(x >= hpd[,1] & x <= hpd[,2]))
    idxs <- if (any(keep)) which(keep) else seq_len(nrow(mat))
    
    # sample one parameter set
    par   <- mat[sample(idxs, 1), 1:nvar]
    nu    <- par["nu"];   r    <- par["r"]
    mt    <- par["mt"];   iota <- par["iota"]
    beta  <- par["beta"]
    gamma <- 1.5;  eta <- -2.0;  kappa <- 0.0
    
    parma <- rbind(
      parma,
      data.frame(vp=vp, sentID=sent, nu, r, mt, iota, beta)
    )
    
    # simulate that sentence
    lfreq <- corpus$lfreq[corpus$sentID == sent]
    if (length(lfreq)==0) next
    traj <- toygen(nu, r, mt, iota, eta, beta, kappa, gamma, lfreq, MODE)[,1:3]
    simAll <- rbind(simAll, cbind(run_id, sent, traj))
  }
  
  # write out sim data
  if (!is.null(simAll)) {
    simDF <- data.frame(
      run    = simAll[,1],
      sentID = simAll[,2],
      word   = simAll[,4],
      fixdur = simAll[,5]
    )
    outFile <- sprintf(
      "D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data/ppsim/ppsim_VP%d.dat",
      vp
    )
    write.table(simDF, file=outFile, row.names=FALSE, sep="\t")
    message(" ✔ Saved simulation for VP ", vp)
  } else {
    warning("No simulation data for VP ", vp)
  }
}

# Save all parameter draws
if (nrow(parma)>0) {
  outPar <- "D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data/ppsim/ppsim.par"
  write.table(parma, file=outPar, row.names=FALSE, sep="\t")
  message("✔ Saved parameter sets to ", outPar)
} else {
  warning("No parameters saved; check warnings above.")
}

