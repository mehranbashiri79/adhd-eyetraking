# -------------------------------
# 1. Load the MCMC Output Object
# -------------------------------
file_path <- "mcmc/mcmc_VP1142.rda"
load(file_path)
if (exists("mcmc")) {
  cat("The object 'mcmc' has been loaded successfully.\n")
  print(summary(mcmc))
  
  if (is.matrix(mcmc) || is.data.frame(mcmc)) {
    cat("Preview of the first 5 rows of mcmc:\n")
    print(head(mcmc, 5))
  } else {
    cat("Structure of the mcmc object:\n")
    str(mcmc)
  }
} else {
  cat("The object 'mcmc' was not found in the file.\n")
}

# -------------------------------
# 2. Set Up the coda Package
# -------------------------------
if (!require(coda)) {
  install.packages("coda")
  library(coda)
} else {
  library(coda)
}

if (!inherits(mcmc, "mcmc.list")) {
  mcmc <- as.mcmc.list(mcmc)
}

# -----------------------------------------------------
# 3. Subset MCMC Chains to Include Only the Parameters
#    of Interest: nu, r, mt, iota, beta
# -----------------------------------------------------
param_names <- c("nu", "r", "mt", "iota", "beta")
mcmc_params <- lapply(mcmc, function(chain) {
  available_params <- intersect(colnames(chain), param_names)
  as.mcmc(chain[, available_params, drop = FALSE])
})
mcmc_params <- mcmc.list(mcmc_params)

# -------------------------------
# 4. Convergence Diagnostics
# -------------------------------

# 4.1 Gelman-Rubin Diagnostic
cat("\nGelman-Rubin Diagnostic:\n")
gelman_results <- gelman.diag(mcmc_params, autoburnin = FALSE)
print(gelman_results)

# 4.2 Trace Plots for Visual Inspection
cat("\nGenerating trace plots...\n")
plot(mcmc_params, main = "Trace Plots for Parameters")

# 4.3 Autocorrelation Plots
cat("\nGenerating autocorrelation plots...\n")
autocorr.plot(mcmc_params, main = "Autocorrelation Plots for Parameters")

# 4.4 Effective Sample Size
cat("\nEffective Sample Sizes:\n")
ess <- effectiveSize(mcmc_params)
print(ess)

# 4.5 Summary Statistics for Each Parameter
cat("\nSummary of the MCMC Chains (Parameters Only):\n")
print(summary(mcmc_params))

# 4.6 Heidelberger and Welch Diagnostics
cat("\nHeidelberger and Welch Diagnostics:\n")
heidel_results <- lapply(mcmc_params, heidel.diag)
print(heidel_results)

library(coda)
# assuming `mcmc_params` is your mcmc.list of nu, r, mt, iota, beta:
ess <- effectiveSize(mcmc_params)
print(ess)

# how quickly does autocorrelation die out?
autocorr.diag(mcmc_params)

