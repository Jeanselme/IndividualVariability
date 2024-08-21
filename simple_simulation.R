library(brms)
library(MASS)
library(matrixStats)
library(rstan)
library(rlist)
library(survival)
library(dplyr)

# Generative functions
## Sampling function
sampling_func <- function(n_individuals, n_points = 5, sampling = 'regular') { 
  if (sampling == 'regular') {
    time <- rep(1:n_points, times = n_individuals)
    time <- data.frame(id = rep(1:n_individuals, each = n_points), time = time)
  }
  return(time)
}

## Generate covariates function
covariates_func <- function(sampling, mean, time_dependent, sigma) { 
  covariates <- mvrnorm(n = nrow(sampling), mu = mean, Sigma = sigma)
  # Ensure time-dependent covariates are temporally consistent for each patient
  if (any(time_dependent)) {
    # Loop through each patient
    for (patient_id in unique(sampling$id)) {
      # Get the indices for the current patient
      patient_indices <- which(sampling$id == patient_id)
      patient_data <- covariates[patient_indices, ]
      
      # Ensure the entire row order reflects the sorted time-dependent covariates
      for (i in which(time_dependent)) {
        covariates[patient_indices, ] <- patient_data[order(patient_data[, i]), ]
      }
    }
  }
  return(covariates)
}

## Random effects function (normally distributed) 
random_effects_func <- function(n_individuals, sd1 = 1, sd2 = 1, corr_effects = 0) { 
  # Create covariance for the random effects
  cov_matrix <- matrix(c(sd1^2, corr_effects*sd1*sd2, corr_effects*sd1*sd2, sd2^2), nrow = 2)
  random_effects <- mvrnorm(n = n_individuals, mu = c(0, 0), Sigma = cov_matrix)
  return(random_effects)
}

## Random effects function (normally distributed) 
random_effects_func_non_normal <- function(n_individuals, sd1 = 1, sd2 = 1, corr_effects = 0) { 
  # Create covariance for the random effects
  cov_matrix <- matrix(c(sd1^2, corr_effects*sd1*sd2, corr_effects*sd1*sd2, sd2^2), nrow = 2)
  random_effects <- mvrnorm(n = n_individuals, mu = c(0, 0), Sigma = cov_matrix)
  return(random_effects)
}


## Outcomes function  
outcomes_func <- function(covariates, random_effects, sample_times, beta, tau, time_slope = FALSE) { 
  # Random noise with the random effect
  if (time_slope) {
    omega <- exp(covariates %*% tau + sample_times$time * random_effects[sample_times$id, 2])
  } else {
    omega <- exp(covariates %*% tau + random_effects[sample_times$id, 2])
  }
  noise <- rnorm(nrow(covariates), 0, omega)
  
  # Arbitrary choice of outcomes
  outcomes <- covariates %*% beta + random_effects[sample_times$id, 1] + noise

  return(outcomes)
}



# Compute error for each model
## Function rmse
rerr <- function(truth, estimate) {
  if (truth != 0) {
    return(mean(abs((truth - estimate) / truth)))
  } else {
    return(mean(abs(truth - estimate)))
  }
}

error_coeff <- function(fit, columns, beta, prefix = 'b_') {
  abs_diffs <- list()
  relative <- list()
  in_bounds <- list()
  for (i in seq_along(beta)) {
    col <- columns[i]
    tryCatch({
      summary_stats <- posterior_summary(fit, variable = paste0(prefix, col))
      estimate <- summary_stats[, 'Estimate']
      lower_bound <- summary_stats[, 'Q2.5']
      upper_bound <- summary_stats[, 'Q97.5']
      
      abs_diffs[[col]] <- abs(beta[i] - estimate)
      relative[[col]] <- rerr(beta[i], estimate)
      in_bounds[[col]] <- as.numeric(beta[i] >= lower_bound & beta[i] <= upper_bound)

    }, error = function(e) {abs_diffs[[col]] <- NA; in_bounds[[col]] <- NA})
  }
  return(list(error=unlist(abs_diffs), relative=unlist(relative), coverage=unlist(in_bounds)))
}

## Function evaluate (silent if error as some methods do not compute all)
evaluate <- function(fit, columns, beta, tau, corr, random_effects, outcomes, indices) {
  errors <- list(time = sum(rstan::get_elapsed_time(fit$fit)))

  # RMSE beta
  beta_estimate = error_coeff(fit, columns, beta)
  errors$beta <- beta_estimate$error
  errors$beta_relative <- beta_estimate$relative
  errors$beta_coverage <- beta_estimate$coverage

  # RMSE tau
  tau_estimate = error_coeff(fit, columns, tau, 'b_sigma_')
  errors$tau <- tau_estimate$error
  errors$tau_relative <- tau_estimate$relative
  errors$tau_coverage <- tau_estimate$coverage

  # RMSE corr
  try(errors$corr <- rerr(corr, posterior_summary(fit, variable='cor_id__Intercept__sigma_Intercept')[, 'Estimate']), silent = TRUE)

  # RMSE random effects
  try(errors$re_intercept <- rerr(random_effects[, 1],  posterior_summary(fit, variable='r_id')[, 'Estimate']), silent = TRUE)
  try(errors$re_sigma <- rerr(random_effects[, 2], posterior_summary(fit, variable='r_id__sigma')[, 'Estimate']), silent = TRUE)
  try(errors$sd_re_intercept <- rerr(1, posterior_summary(fit, variable='sd_id__Intercept')[, 'Estimate']), silent = TRUE)
  try(errors$sd_re_sigma <- rerr(1, posterior_summary(fit, variable='sd_id__sigma_Intercept')[, 'Estimate']), silent = TRUE)

  # RMSE outcomes
  errors$outcomes <- rerr(predict(fit)[, 'Estimate'][indices], outcomes[indices])

  return(errors)
}

## Append function to average across simulations
append <- function(previous, evaluation) {
  errors <- list()
  for (quantity in names(evaluation)) {
    if (!quantity %in% names(previous)) {
      previous[[quantity]] <- c()
    }
    errors[[quantity]] <- rbind(previous[[quantity]], evaluation[[quantity]])
  }
  return(errors)
}

summarise <- function(evaluation) {
  errors <- list()
  for (model in names(evaluation)) {
    errors[[model]] <- list()
    for (quantity in names(evaluation[[model]]))
      errors[[model]][[quantity]] <- list(
        value=colMeans(evaluation[[model]][[quantity]]),
        std=colSds(evaluation[[model]][[quantity]], na.rm = TRUE)
      )
  }
  return(errors)
}

# Simulation studies
## Full loop train and evaluation
simulation_mm <- function(n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov) {
  evaluation = list(melsm = list(), melsm_incorr_sigma = list(), 
                    melsm_incorr_out = list(), melsm_incorr_sigmainter = list(), 
                    melsm_all = list(), mm = list(), mmall = list())
  for (i in 1:n_sim) {

    # Fix seed
    set.seed(i)

    # Generate data
    sample_times <- sampling_func(n_individuals, n_points)
    covariates <- covariates_func(sample_times, covariate_mean, time_dependent, covariate_cov)
    random_effects <- random_effects_func(n_individuals, corr_effects = corr)
    outcomes <- outcomes_func(covariates, random_effects, sample_times, beta, tau)

    last_time_indices <- sample_times %>%
      group_by(id) %>%
      summarize(index_last = (time == max(time)))
    last_time_indices <- last_time_indices$index_last

    data <- data.frame(sample_times, covariates, outcomes = outcomes)

    # Fit MELSM - Correctly specified
    formula <- bf(
      outcomes ~ age + albumin + (1|C|id),
      sigma ~ trig + platelet + (1|C|id), 
      family = gaussian()
    )

    fit <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm <- append(evaluation$melsm, evaluate(fit, columns, beta, tau, corr, random_effects, outcomes, last_time_indices))

    # Fit MM
    formula <- bf(
      outcomes ~ age + albumin + (1|id),
      family = gaussian()
    )

    fit <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$mm <- append(evaluation$mm, evaluate(fit, columns, beta, tau, corr, random_effects, outcomes, last_time_indices))
  }
  return(summarise(evaluation))
}

simulation_melsm <- function(n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov) {
  evaluation = list(melsm = list(), melsm_incorr_sigma = list(), 
                    melsm_incorr_out = list(), melsm_incorr_sigmainter = list(), 
                    melsm_all = list(), mm = list(), mmall = list())
  for (i in 1:n_sim) {

    # Fix seed
    set.seed(i)

    # Generate data
    sample_times <- sampling_func(n_individuals, n_points)
    covariates <- covariates_func(sample_times, covariate_mean, time_dependent, covariate_cov)
    random_effects <- random_effects_func(n_individuals, corr_effects = corr)
    outcomes <- outcomes_func(covariates, random_effects, sample_times, beta, tau)

    last_time_indices <- sample_times %>%
      group_by(id) %>%
      summarize(index_last = (time == max(time)))
    last_time_indices <- last_time_indices$index_last

    data <- data.frame(sample_times, covariates, outcomes = outcomes)

    # Fit MELSM - Correctly specified
    formula <- bf(
      outcomes ~ age + albumin + (1|C|id),
      sigma ~ trig + platelet + (1|C|id), 
      family = gaussian()
    )

    fit <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm <- append(evaluation$melsm, evaluate(fit, columns, beta, tau, corr, random_effects, outcomes, last_time_indices))

    # Fit MELSM - Incorrectly specified sigma
    formula <- bf(
      outcomes ~ age + albumin + (1|C|id),
      sigma ~ age + albumin + (1|C|id), 
      family = gaussian()
    )

    fit <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm_incorr_sigma <- append(evaluation$melsm_incorr_sigma, evaluate(fit, columns, beta, tau, corr, random_effects, outcomes, last_time_indices))

    # Fit MELSM - Incorrectly specified outcome
    formula <- bf(
      outcomes ~ trig + platelet + (1|C|id),
      sigma ~ trig + platelet + (1|C|id), 
      family = gaussian()
    )

    fit <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm_incorr_out <- append(evaluation$melsm_incorr_out, evaluate(fit, columns, beta, tau, corr, random_effects, outcomes, last_time_indices))

    # Fit MELSM - ALL
    formula <- bf(
      outcomes ~ age + albumin + trig + platelet + (1|C|id),
      sigma ~ age + albumin + trig + platelet + (1|C|id), 
      family = gaussian()
    )

    fit <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm_all <- append(evaluation$melsm_all, evaluate(fit, columns, beta, tau, corr, random_effects, outcomes, last_time_indices))


    # Fit MELSM - No sigma intercept
    formula <- bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ trig + platelet, 
      family = gaussian()
    )

    fit <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm_incorr_sigmainter <- append(evaluation$melsm_incorr_sigmainter, evaluate(fit, columns, beta, tau, corr, random_effects, outcomes, last_time_indices))

    fit <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$mm <- append(evaluation$mm, evaluate(fit, columns, beta, tau, corr, random_effects, outcomes, last_time_indices))
  }
  return(summarise(evaluation))
}

simulation_time <- function(n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov) {
  evaluation = list(melsm = list(), melsm_incorr_sigma = list(), 
                    melsm_incorr_out = list(), melsm_incorr_sigmainter = list(), 
                    melsm_all = list(), mm = list(), mmall = list())
  for (i in 1:n_sim) {

    # Fix seed
    set.seed(i)

    # Generate data
    sample_times <- sampling_func(n_individuals, n_points)
    covariates <- covariates_func(sample_times, covariate_mean, time_dependent, covariate_cov)
    random_effects <- random_effects_func(n_individuals, corr_effects = corr)
    outcomes <- outcomes_func(covariates, random_effects, sample_times, beta, tau, TRUE)

    last_time_indices <- sample_times %>%
      group_by(id) %>%
      summarize(index_last = (time == max(time)))
    last_time_indices <- last_time_indices$index_last

    data <- data.frame(sample_times, covariates, outcomes = outcomes)

    # Fit MELSM - Correctly specified
    formula <- bf(
      outcomes ~ age + albumin + (time|C|id),
      sigma ~ trig + platelet + (1|C|id), 
      family = gaussian()
    )

    fit <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm <- append(evaluation$melsm, evaluate(fit, columns, beta, tau, corr, random_effects, outcomes, last_time_indices))

    # Fit MELSM - Incorrectly specified 
    formula <- bf(
      outcomes ~ age + albumin + (1|C|id),
      sigma ~ trig + platelet + (1|C|id), 
      family = gaussian()
    )

    fit <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm_incorr_sigma <- append(evaluation$melsm_incorr_sigma, evaluate(fit, columns, beta, tau, corr, random_effects, outcomes, last_time_indices))
  }
  return(summarise(evaluation))
}


# Default parameters
## Data
columns <- c("age", "albumin", "trig", "platelet")
time_dependent <- c(TRUE, FALSE, FALSE, FALSE) # To unsure time consistency over multiple draw
pbc_scale <- scale(pbc[columns])
covariate_mean <- colMeans(pbc_scale, na.rm = TRUE)
covariate_cov <- cov(pbc_scale, use="complete.obs")

## True model (CANNOT CHANGE 0 without changing experiment)
beta <- c(1, 0.5, 0., 0.)
tau <- c(0., 0, 0.8, 0.5)
corr <- 0

## Population
n_individuals <- 200
n_points <- 15

## Simulation
n_sim <- 10

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  print('Run all')
  run <- -1
} else {
  run <- as.numeric(args[1])
}

# Simulation study: Impact of number of individuals
if ((run == -1)|(run == 1)) {
  n_individuals_list <- c(100, 200, 300)
  for (n_individuals_exp in n_individuals_list) {
    print(paste("Simulating for", n_individuals_exp, "individuals"))
    sim = simulation_mm(n_sim, n_individuals_exp, n_points, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov)
    list.save(sim, file = paste0("results/individuals", n_individuals_exp, '.json'))
  }
}

# Simulation study: Impact of number of points
if ((run == -1)|(run == 2)) {
  n_points_list <- c(5, 10, 20)
  for (n_points_exp in n_points_list) {
    print(paste("Simulating for", n_points_exp, "points"))
    sim = simulation_mm(n_sim, n_individuals, n_points_exp, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov)
    list.save(sim, file = paste0("results/points", n_points_exp, '.json'))
  }
}

# Simulation study: Impact of number of corr
if ((run == -1)|(run == 3)) {
  rho_list <- c(-0.5, -0.25, 0.25, 0.5)
  for (corr_exp in rho_list) {
    print(paste("Simulating for", corr_exp, "correlation"))
    sim = simulation_mm(n_sim, n_individuals, n_points, corr_exp, columns, beta, tau, covariate_mean, time_dependent, covariate_cov)
    list.save(sim, file = paste0("results/corr", corr_exp, '.json'))
  }
}

if ((run == -1)|(run == 4)) {
  print(paste("Simulating for MELSM misspecification"))
  sim = simulation_melsm(n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov)
  list.save(sim, file = paste0("results/misspecification", n_points_exp, '.json'))
}

if ((run == -1)|(run == 5)) {
  print(paste("Simulating for time dependent random effects"))
  sim = simulation_time(n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov)
  list.save(sim, file = paste0("results/time", n_points_exp, '.json'))
}