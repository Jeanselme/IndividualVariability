library(brms)
library(MASS)
library(rstan)
library(rlist)
library(survival)

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
covariates_func <- function(sampling, mean, sigma) { 
  covariates <- mvrnorm(n = nrow(sampling), mu = mean, Sigma = sigma)
  return(covariates)
}

## Random effects function (normally distributed) 
random_effects_func <- function(n_individuals, sd1 = 1, sd2 = 1, corr_effects = 0) { 
  # Create covariance for the random effects
  cov_matrix <- matrix(c(sd1^2, corr_effects*sd1*sd2, corr_effects*sd1*sd2, sd2^2), nrow = 2)
  random_effects <- mvrnorm(n = n_individuals, mu = c(0, 0), Sigma = cov_matrix)
  return(random_effects)
}

## Outcomes function  
outcomes_func <- function(covariates, random_effects, sample_times, beta, tau) { 
  # Random noise with the random effect
  omega <- exp(covariates %*% tau + random_effects[sample_times$id, 2])
  noise <- rnorm(nrow(covariates), 0, omega)
  
  # Arbitrary choice of outcomes
  outcomes <- covariates %*% beta + random_effects[sample_times$id, 1] + noise
  return(outcomes)
}



# Compute error for each model
## Function rmse
rmse <- function(x) {
  sqrt(mean(x^2))
}

error_coeff <- function(fit, columns, beta, prefix = 'b_') {
  abs_diffs <- list()
  in_bounds <- list()
  for (i in seq_along(beta)) {
    col <- columns[i]
    tryCatch({
      summary_stats <- posterior_summary(fit, variable = paste0(prefix, col))
      estimate <- summary_stats[, 'Estimate']
      lower_bound <- summary_stats[, 'Q2.5']
      upper_bound <- summary_stats[, 'Q97.5']
      
      abs_diffs[[col]] <- abs(beta[i] - estimate)
      in_bounds[[col]] <- beta[i] >= lower_bound & beta[i] <= upper_bound

    }, error = function(e) {})
  }
  return(list(error=unlist(abs_diffs), coverage=unlist(in_bounds)))
}

## Function evaluate (silent if error as some methods do not compute all)
evaluate <- function(fit, columns, beta, tau, corr, random_effects, outcomes) {
  errors <- list(time = sum(rstan::get_elapsed_time(fit$fit)))

  # RMSE beta
  beta_estimate = error_coeff(fit, columns, beta)
  errors$beta <- beta_estimate$error
  errors$beta_coverage <- beta_estimate$coverage

  # RMSE tau
  tau_estimate = error_coeff(fit, columns, tau, 'b_sigma_')
  errors$tau <- tau_estimate$error
  errors$tau_coverage <- tau_estimate$coverage

  # RMSE corr
  try(errors$corr <- abs(corr - posterior_summary(fit, variable='cor_id__Intercept__sigma_Intercept')[, 'Estimate']), silent = TRUE)

  # RMSE random effects
  try(errors$re_intercept <- rmse(random_effects[, 1] - posterior_summary(fit, variable='r_id')[, 'Estimate']), silent = TRUE)
  try(errors$re_sigma <- rmse(random_effects[, 2] - posterior_summary(fit, variable='r_id__sigma')[, 'Estimate']), silent = TRUE)

  # RMSE outcomes
  errors$outcomes <- rmse(fitted(fit)[, 'Estimate'] - outcomes)

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
        value=mean(evaluation[[model]][[quantity]]),
        std=sd(evaluation[[model]][[quantity]], na.rm = TRUE)
      )
  }
  return(errors)
}


# Simulation studies
## Full loop train and evaluation
simulation <- function(n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, covariate_cov) {
  evaluation = list(melsm = list(), melsm_incorr_sigma = list(), 
                    melsm_incorr_out = list(), melsm_incorr_sigmainter = list(), 
                    melsm_all = list(), mm = list(), mmall = list())
  for (i in 1:n_sim) {

    # Fix seed
    set.seed(i)

    # Generate data
    sample_times <- sampling_func(n_individuals, n_points)
    covariates <- covariates_func(sample_times, covariate_mean, covariate_cov)
    random_effects <- random_effects_func(n_individuals, corr_effects = corr)
    outcomes <- outcomes_func(covariates, random_effects, sample_times, beta, tau)

    data <- data.frame(sample_times, covariates, outcomes = outcomes)

    # Fit MELSM - Correctly specified
    formula <- bf(
      outcomes ~ age + albumin + (1|C|id),
      sigma ~ trig + platelet + (1|C|id), 
      family = gaussian()
    )

    fit <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm <- append(evaluation$melsm, evaluate(fit, columns, beta, tau, corr, random_effects, outcomes))

    # Fit MELSM - Incorrectly specified sigma
    formula <- bf(
      outcomes ~ age + albumin + (1|C|id),
      sigma ~ age + albumin + (1|C|id), 
      family = gaussian()
    )

    fit_incorr <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm_incorr_sigma <- append(evaluation$melsm_incorr_sigma, evaluate(fit_incorr, columns, beta, tau, corr, random_effects, outcomes))

    # Fit MELSM - Incorrectly specified outcome
    formula <- bf(
      outcomes ~ trig + platelet + (1|C|id),
      sigma ~ trig + platelet + (1|C|id), 
      family = gaussian()
    )

    fit_incorr <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm_incorr_out <- append(evaluation$melsm_incorr_out, evaluate(fit_incorr, columns, beta, tau, corr, random_effects, outcomes))

    # Fit MELSM - ALL
    formula <- bf(
      outcomes ~ . + (1|C|id),
      sigma ~ . + (1|C|id), 
      family = gaussian()
    )

    fit_all <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm_all <- append(evaluation$melsm_all, evaluate(fit_all, columns, beta, tau, corr, random_effects, outcomes))


    # Fit MELSM - No sigma intercept
    formula <- bf(
      outcomes ~ age + albumin + (1|C|id),
      sigma ~ trig + platelet, 
      family = gaussian()
    )

    fit_nointer <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$melsm_incorr_sigmainter <- append(evaluation$melsm_incorr_sigmainter, evaluate(fit_nointer, columns, beta, tau, corr, random_effects, outcomes))


    # Fit MM
    formula <- bf(
      outcomes ~ age + albumin + (1|id),
      family = gaussian()
    )

    fit_mm <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$mm <- append(evaluation$mm, evaluate(fit_mm, columns, beta, tau, corr, random_effects, outcomes))

    # Fit MM all
    formula <- bf(
      outcomes ~ . + (1|id),
      family = gaussian()
    )

    fit_mmall <- brm(formula, data, seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
    evaluation$mmall <- append(evaluation$mmall, evaluate(fit_mmall, columns, beta, tau, corr, random_effects, outcomes))
  }
  return(summarise(evaluation))
}

# Default parameters
## Data
columns <- c("age", "albumin", "trig", "platelet")
covariate_mean <- colMeans(pbc[columns], na.rm = TRUE)
covariate_cov <- cor(pbc[columns], use="complete.obs")

## True model (CANNOT CHANGE 0 without changing experiment)
beta <- c(1, 0.5, 0., 0.)
tau <- c(0., 0, 0.8, 0.5)
corr <- 0

## Population
n_individuals <- 200
n_points <- 15

## Simulation
n_sim <- 25


# Simulation study: Impact of number of individuals
n_individuals_list <- c(100, 200, 300)
for (n_individuals_exp in n_individuals_list) {
  print(paste("Simulating for", n_individuals_exp, "individuals"))
  sim = simulation(n_sim, n_individuals_exp, n_points, corr, beta, tau, covariate_mean, covariate_cov)
  list.save(sim, file = paste0("results/individuals", n_individuals_exp, '.json'))
}

# Simulation study: Impact of number of points
n_points_list <- c(5, 15, 25)
for (n_points_exp in n_points_list) {
  print(paste("Simulating for", n_points_exp, "points"))
  sim = simulation(n_sim, n_individuals, n_points_exp, corr, beta, tau, covariate_mean, covariate_cov)
  list.save(sim, file = paste0("results/points", n_points_exp, '.json'))
}

# Simulation study: Impact of number of corr
rho_list <- c(-0.5, -0.25, 0.25, 0.5)
for (corr_exp in rho_list) {
  print(paste("Simulating for", corr_exp, "correlation"))
  sim = simulation(n_sim, n_individuals, n_points, corr_exp, columns, beta, tau, covariate_mean, covariate_cov)
  list.save(sim, file = paste0("results/corr", corr_exp, '.json'))
}