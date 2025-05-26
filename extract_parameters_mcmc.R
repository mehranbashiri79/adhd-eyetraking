library(coda)

# Paths
mcmc_dir <- "D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data/mcmc"
data_dir <- "D:/ACIT Master's Thesis/New data/swift/R-Code-parameter-estimation-from-experimental-data/expdata"
output_path <- file.path(mcmc_dir, "mcmc_features_with_sd.csv")

# List all MCMC files
mcmc_files <- list.files(mcmc_dir, pattern = "^mcmc_VP\\d+\\.rda$", full.names = TRUE)

# Initialize results storage
results <- data.frame()

for (file in mcmc_files) {
  load(file)
  
  if (!exists("mcmc")) next
  
  vp_id <- as.numeric(gsub("\\D", "", basename(file)))
  if (!inherits(mcmc, "mcmc.list")) mcmc <- as.mcmc.list(mcmc)
  
  param_names <- c("nu", "r", "mt", "iota", "beta", "LL")
  mcmc_df <- do.call(rbind, mcmc)
  
  if (!all(param_names %in% colnames(mcmc_df))) next
  
  # Extract means and SDs
  means <- colMeans(mcmc_df[, param_names])
  sds   <- apply(mcmc_df[, param_names], 2, sd)
  
  # Normalize LL by number of fixations
  exp_file <- file.path(data_dir, paste0(vp_id, ".dat"))
  if (!file.exists(exp_file)) next
  
  exp_data <- read.table(exp_file, header = TRUE)
  n_fix <- nrow(exp_data)
  norm_ll <- means["LL"] / n_fix
  norm_ll_sd <- sds["LL"] / n_fix
  
  # Save result row
  result_row <- data.frame(
    participant = vp_id,
    nu_mean = means["nu"],
    nu_sd   = sds["nu"],
    r_mean  = means["r"],
    r_sd    = sds["r"],
    mt_mean = means["mt"],
    mt_sd   = sds["mt"],
    iota_mean = means["iota"],
    iota_sd   = sds["iota"],
    beta_mean = means["beta"],
    beta_sd   = sds["beta"],
    norm_ll_mean = norm_ll,
    norm_ll_sd   = norm_ll_sd
  )
  
  results <- rbind(results, result_row)
}

# Save to CSV
write.csv(results, output_path, row.names = FALSE)
cat("✅ Features saved to:", output_path, "\n")

