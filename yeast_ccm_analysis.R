# Yeast Causality Analysis: Reproducing Pao et al. (2026)
# This script implements Empirical Dynamic Modeling (EDM) and Convergent Cross Mapping (CCM)
# for detecting causation without correlation in transcriptional networks

# Install required packages if needed
required_packages <- c("tidyverse", "purrr", "tidyr", "ggplot2")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ============================================================================
# SECTION 1: UTILITY FUNCTIONS FOR ATTRACTOR RECONSTRUCTION
# ============================================================================

# Takens embedding: Create shadow attractor from time-lagged coordinates
create_shadow_attractor <- function(time_series, embedding_dim = 3, time_delay = 1) {
  """
  Reconstruct a shadow attractor using time-delayed embedding (Takens theorem)
  
  Args:
    time_series: numeric vector of gene expression values
    embedding_dim: number of embedding dimensions (E)
    time_delay: lag in time steps
  
  Returns:
    matrix with embedding_dim columns representing the attractor
  """
  n <- length(time_series)
  max_idx <- n - (embedding_dim - 1) * time_delay
  
  attractor <- matrix(NA, nrow = max_idx, ncol = embedding_dim)
  
  for (i in 1:embedding_dim) {
    idx_start <- 1 + (i - 1) * time_delay
    idx_end <- max_idx + (i - 1) * time_delay
    attractor[, i] <- time_series[idx_start:idx_end]
  }
  
  return(attractor)
}

# Simplex projection: Nearest-neighbor prediction
simplex_projection <- function(attractor, time_series, tp = 1) {
  """
  Use simplex projection to assess predictability
  
  Args:
    attractor: reconstructed attractor matrix
    time_series: original time series
    tp: prediction horizon (time steps ahead)
  
  Returns:
    list with predictions and skill metrics
  """
  n <- nrow(attractor)
  E <- ncol(attractor)
  
  predictions <- numeric(n - tp)
  observations <- numeric(n - tp)
  
  for (i in 1:(n - tp)) {
    current_point <- attractor[i, ]
    
    # Calculate distances to all other points
    distances <- sqrt(rowSums((attractor - matrix(current_point, 
                                                    nrow = n, 
                                                    ncol = E, 
                                                    byrow = TRUE))^2))
    
    # Find E+1 nearest neighbors
    nn_indices <- order(distances)[1:(E + 1)]
    nn_distances <- distances[nn_indices]
    
    # Avoid division by zero
    nn_distances[nn_distances == 0] <- 1e-10
    
    # Weight inversely by distance
    weights <- exp(-nn_distances / mean(nn_distances))
    weights <- weights / sum(weights)
    
    # Predict based on where neighbors go
    target_indices <- nn_indices + tp
    target_indices <- target_indices[target_indices <= length(time_series)]
    
    if (length(target_indices) > 0) {
      predictions[i] <- weighted.mean(
        time_series[target_indices], 
        weights[1:length(target_indices)]
      )
      observations[i] <- time_series[i + tp]
    }
  }
  
  # Remove NA values
  valid_idx <- !is.na(predictions) & !is.na(observations)
  predictions <- predictions[valid_idx]
  observations <- observations[valid_idx]
  
  # Calculate skill as Pearson correlation
  skill <- cor(predictions, observations, use = "complete.obs")
  
  return(list(
    predictions = predictions,
    observations = observations,
    skill = skill,
    rmse = sqrt(mean((predictions - observations)^2, na.rm = TRUE))
  ))
}

# S-map: Nonlinear forecasting with state-dependent coefficients
smap <- function(attractor, time_series, theta_values = seq(0, 2, 0.1), tp = 1) {
  """
  S-map test for nonlinearity using locally-weighted regression
  
  Args:
    attractor: reconstructed attractor
    time_series: original time series
    theta_values: range of nonlinearity parameters to test
    tp: prediction horizon
  
  Returns:
    list with optimal theta and nonlinearity test results
  """
  n <- nrow(attractor)
  E <- ncol(attractor)
  
  skills <- numeric(length(theta_values))
  
  for (theta_idx in seq_along(theta_values)) {
    theta <- theta_values[theta_idx]
    predictions <- numeric(n - tp)
    
    for (i in 1:(n - tp)) {
      current_point <- attractor[i, ]
      
      # Calculate distances
      distances <- sqrt(rowSums((attractor - matrix(current_point, 
                                                      nrow = n, 
                                                      ncol = E, 
                                                      byrow = TRUE))^2))
      
      # Exponential weighting based on theta (nonlinearity parameter)
      avg_distance <- mean(distances)
      weights <- exp(-theta * distances / avg_distance)
      weights <- weights / sum(weights)
      
      # Weighted prediction
      target_indices <- (1:n) + tp
      valid_targets <- target_indices[target_indices <= length(time_series)]
      valid_weights <- weights[1:length(valid_targets)]
      valid_weights <- valid_weights / sum(valid_weights)
      
      if (length(valid_targets) > 0) {
        predictions[i] <- weighted.mean(
          time_series[valid_targets],
          valid_weights
        )
      }
    }
    
    # Calculate skill
    valid_idx <- !is.na(predictions)
    if (sum(valid_idx) > 0) {
      skills[theta_idx] <- cor(predictions[valid_idx], 
                               time_series[1:(n-tp)][valid_idx],
                               use = "complete.obs")
    }
  }
  
  optimal_idx <- which.max(skills)
  optimal_theta <- theta_values[optimal_idx]
  
  # Test for nonlinearity: S-map should improve over linear (theta=0)
  linear_skill <- skills[which.min(abs(theta_values))]
  nonlinear_skill <- max(skills, na.rm = TRUE)
  
  is_nonlinear <- nonlinear_skill > linear_skill
  
  return(list(
    optimal_theta = optimal_theta,
    optimal_skill = nonlinear_skill,
    linear_skill = linear_skill,
    skills = skills,
    theta_values = theta_values,
    is_nonlinear = is_nonlinear,
    nonlinearity_strength = nonlinear_skill - linear_skill
  ))
}

# ============================================================================
# SECTION 2: CONVERGENT CROSS MAPPING (CCM) FOR CAUSAL INFERENCE
# ============================================================================

convergent_cross_mapping <- function(cause_series, effect_series, 
                                     embedding_dim = 3, time_delay = 1,
                                     library_sizes = NULL) {
  """
  Convergent Cross Mapping: Detect causal relationships in time series
  
  The key insight: if Y causes X, then the attractor reconstructed from X's
  time series contains information about Y. This allows us to predict Y from
  the state space of X.
  
  Args:
    cause_series: time series of hypothesized cause (Y)
    effect_series: time series of hypothesized effect (X)
    embedding_dim: embedding dimension
    time_delay: time delay for reconstruction
    library_sizes: library sizes to test for convergence
  
  Returns:
    list with CCM results showing directional causality
  """
  
  if (is.null(library_sizes)) {
    max_lib <- min(length(cause_series), length(effect_series)) - 
               (embedding_dim - 1) * time_delay - 1
    library_sizes <- round(seq(10, max_lib, length.out = 10))
  }
  
  # Reconstruct shadow attractor from effect series
  attractor_effect <- create_shadow_attractor(effect_series, 
                                             embedding_dim, 
                                             time_delay)
  
  ccm_skills <- numeric(length(library_sizes))
  
  for (lib_idx in seq_along(library_sizes)) {
    lib_size <- library_sizes[lib_idx]
    
    if (lib_size > nrow(attractor_effect)) {
      ccm_skills[lib_idx] <- NA
      next
    }
    
    # Sample library subset
    lib_indices <- sample(1:nrow(attractor_effect), lib_size)
    lib_attractor <- attractor_effect[lib_indices, ]
    lib_cause <- cause_series[lib_indices]
    
    # Cross-predict cause from effect's state space
    predictions <- numeric(nrow(attractor_effect))
    
    for (i in 1:nrow(attractor_effect)) {
      point <- attractor_effect[i, ]
      
      # Find nearest neighbors in library
      distances <- sqrt(rowSums((lib_attractor - 
                                matrix(point, nrow = lib_size, 
                                      ncol = ncol(lib_attractor), 
                                      byrow = TRUE))^2))
      
      nn_idx <- which.min(distances)
      predictions[i] <- lib_cause[nn_idx]
    }
    
    # Calculate cross-map skill
    valid_idx <- !is.na(predictions)
    if (sum(valid_idx) > 2) {
      ccm_skills[lib_idx] <- cor(predictions[valid_idx], 
                                 cause_series[valid_idx],
                                 use = "complete.obs")
    }
  }
  
  return(list(
    library_sizes = library_sizes,
    ccm_skills = ccm_skills,
    convergence = mean(diff(ccm_skills), na.rm = TRUE) > 0,
    skill_at_max = max(ccm_skills, na.rm = TRUE)
  ))
}

# ============================================================================
# SECTION 3: TEST DATA GENERATION (Synthetic Yeast-like Data)
# ============================================================================

generate_yeast_oscillations <- function(n_timepoints = 57, n_genes = 100) {
  """
  Generate synthetic yeast cell cycle gene expression data
  with realistic nonlinear dynamics
  
  Args:
    n_timepoints: number of time points
    n_genes: number of genes to simulate
  
  Returns:
    data frame with time series for each gene
  """
  
  set.seed(42)
  
  t <- seq(0, 2 * pi, length.out = n_timepoints)
  
  # Create base oscillations for key cell cycle genes
  CLN3 <- 0.6 + 0.4 * sin(t)
  SWI4 <- 0.5 + 0.5 * sin(t + 0.3)
  CLB2 <- 0.4 + 0.6 * cos(t + 1.0)
  WHI5 <- 0.7 - 0.5 * cos(t + 0.5)
  
  # Add nonlinearity: state-dependent dynamics
  # WHI5 represses SWI4 in a state-dependent way
  SWI4_nonlinear <- SWI4 * (1 - 0.7 * WHI5)
  
  # CLN3 and SWI4 interact nonlinearly at G1/S checkpoint
  CLN3_nonlinear <- CLN3 * (0.5 + 0.5 * SWI4_nonlinear)
  
  # Create additional genes with partial correlation
  data <- data.frame(
    time = 1:n_timepoints,
    CLN3 = CLN3_nonlinear + rnorm(n_timepoints, 0, 0.05),
    SWI4 = SWI4_nonlinear + rnorm(n_timepoints, 0, 0.05),
    CLB2 = CLB2 + rnorm(n_timepoints, 0, 0.05),
    WHI5 = WHI5 + rnorm(n_timepoints, 0, 0.05)
  )
  
  # Add more genes with varying correlation patterns
  for (i in 1:n_genes) {
    # Some genes follow cell cycle
    phase_shift <- runif(1, 0, 2*pi)
    amplitude <- runif(1, 0.3, 0.8)
    data[[paste0("Gene_", i)]] <- amplitude * sin(t + phase_shift) + 
                                   0.5 + rnorm(n_timepoints, 0, 0.05)
  }
  
  return(data)
}

# ============================================================================
# SECTION 4: MAIN ANALYSIS PIPELINE
# ============================================================================

run_yeast_edm_analysis <- function(time_series_data, gene_name, 
                                   embedding_dims = 3:5) {
  """
  Complete EDM analysis for a single gene:
  1. Test predictability (Simplex)
  2. Test for nonlinearity (S-map)
  3. Assess convergence
  """
  
  ts <- time_series_data[[gene_name]]
  
  # Remove NAs
  ts <- na.omit(as.numeric(ts))
  
  results <- list(
    gene = gene_name,
    n_timepoints = length(ts)
  )
  
  # Test different embedding dimensions
  best_skill <- -1
  best_E <- NA
  
  for (E in embedding_dims) {
    if (length(ts) <= (E - 1) + 1) next  # Need enough points
    
    attractor <- create_shadow_attractor(ts, E, time_delay = 1)
    simplex_result <- simplex_projection(attractor, ts, tp = 1)
    
    if (!is.na(simplex_result$skill) && simplex_result$skill > best_skill) {
      best_skill <- simplex_result$skill
      best_E <- E
      results$best_simplex <- simplex_result
    }
  }
  
  results$best_embedding_dim <- best_E
  results$predictability_skill <- best_skill
  results$is_predictable <- best_skill > 0.1  # Significance threshold
  
  # If predictable, test for nonlinearity
  if (results$is_predictable) {
    attractor <- create_shadow_attractor(ts, best_E, time_delay = 1)
    smap_result <- smap(attractor, ts, tp = 1)
    results$smap <- smap_result
    results$is_nonlinear <- smap_result$is_nonlinear
  }
  
  return(results)
}

# ============================================================================
# SECTION 5: EXECUTE ANALYSIS
# ============================================================================

cat("========================================\n")
cat("Yeast CCM Analysis: Causation without Correlation\n")
cat("========================================\n\n")

# Generate synthetic data
cat("Generating synthetic yeast cell cycle data...\n")
yeast_data <- generate_yeast_oscillations(n_timepoints = 57, n_genes = 50)

cat("Data dimensions:", dim(yeast_data), "\n")
cat("Time points:", length(yeast_data$time), "\n\n")

# Analyze key genes
key_genes <- c("CLN3", "SWI4", "CLB2", "WHI5")

cat("Testing for predictability and nonlinearity...\n\n")
edm_results <- list()

for (gene in key_genes) {
  if (gene %in% colnames(yeast_data)) {
    cat(sprintf("Analyzing %s...\n", gene))
    result <- run_yeast_edm_analysis(yeast_data, gene)
    edm_results[[gene]] <- result
    
    cat(sprintf("  Predictability skill: %.3f\n", result$predictability_skill))
    cat(sprintf("  Is predictable: %s\n", result$is_predictable))
    if (result$is_predictable) {
      cat(sprintf("  Is nonlinear: %s\n", result$is_nonlinear))
      cat(sprintf("  Nonlinearity strength: %.3f\n", 
                  result$smap$nonlinearity_strength))
    }
    cat("\n")
  }
}

# Test CCM for causal relationship: WHI5 -> SWI4
cat("Testing causal relationship: WHI5 -> SWI4\n")
cat("(If WHI5 causes SWI4, we should predict SWI4 from WHI5's state space)\n\n")

ccm_result <- convergent_cross_mapping(
  cause_series = yeast_data$WHI5,
  effect_series = yeast_data$SWI4,
  embedding_dim = 3
)

cat("CCM Skills by library size:\n")
for (i in seq_along(ccm_result$library_sizes)) {
  cat(sprintf("  Library size %d: skill = %.3f\n", 
              ccm_result$library_sizes[i],
              ccm_result$ccm_skills[i]))
}

cat(sprintf("\nConvergence detected: %s\n", ccm_result$convergence))
cat(sprintf("Maximum skill: %.3f\n\n", ccm_result$skill_at_max))

# Calculate linear correlation for comparison
correlation <- cor(yeast_data$WHI5, yeast_data$SWI4)
cat(sprintf("For comparison, linear correlation (WHI5, SWI4): %.3f\n", 
            correlation))

cat("\n========================================\n")
cat("Analysis complete!\n")
cat("========================================\n")
