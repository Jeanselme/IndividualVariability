library(brms)
library(MASS)

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

## Function evaluate (silent if error as some methods do not compute all)
evaluate <- function(fit, beta, tau, corr, random_effects, outcomes) {
  errors <- list()

  # RMSE beta
  errors$beta <- abs(beta - posterior_summary(fit, pars='b_z')[, 'Estimate'])

  # RMSE tau
  try(errors$tau  <- abs(tau - posterior_summary(fit, pars='b_sigma_z')[, 'Estimate']), silent = TRUE)

  # RMSE corr
  try(errors$corr <- abs(corr - posterior_summary(fit, variable='cor_id__Intercept__sigma_Intercept')[, 'Estimate']), silent = TRUE)

  # RMSE random effects
  try(errors$re_intercept <- rmse(random_effects[, 1] - posterior_summary(fit, variable='r_id')[, 'Estimate']), silent = TRUE)
  try(errors$re_sigma <- rmse(random_effects[, 2] - posterior_summary(fit, variable='r_id__sigma')[, 'Estimate']), silent = TRUE)

  # RMSE outcomes
  errors$outcomes <- rmse(fitted(fit)[, 'Estimate'] - outcomes)

  return(errors)
}

## Summarise function to average across simulations
summarise <- function(previous, evaluation, n_sim) {
  errors <- list()
  for (quantity in names(evaluation)) {
    if (!quantity %in% names(previous)) {
      previous[[quantity]] <- 0
    }
    errors[[quantity]] <- previous[[quantity]] + (evaluation[[quantity]] / n_sim)
  }
  return(errors)
}


# Simulation studies
## Full loop train and evaluation
simulation <- function(n_sim, n_individuals, n_points, corr, beta, tau, covariate_mean, covariate_cov) {
  evaluation = list(melsm = list(), melsm_incorr = list(), mm = list())
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
      outcomes ~ zAge + zBase + (1|corr|id),
      sigma ~ zAge + zBase + (1|corr|id), 
      family = gaussian()
    )

    fit <- brm(
      formula,
      data,
      seed = i,
      warmup = 500,
      iter = 1000,
      chains = 4, 
      cores = 4
    )

    evaluation$melsm <- summarise(evaluation$melsm, evaluate(fit, beta, tau, corr, random_effects, outcomes), n_sim)

    # Fit MELSM - Incorrectly specified
    formula <- bf(
      outcomes ~ zAge + zBase + (1|corr|id),
      sigma ~ zAge + zBase, 
      family = gaussian()
    )

    fit_incorr <- brm(
      formula,
      data,
      seed = i,
      warmup = 500,
      iter = 1000,
      chains = 4, 
      cores = 4
    )

    evaluation$melsm_incorr <- summarise(evaluation$melsm_incorr, evaluate(fit_incorr, beta, tau, corr, random_effects, outcomes), n_sim)

    # Fit MM - Correctly specified
    formula <- bf(
      outcomes ~  zAge + zBase + (1|corr|id),
      family = gaussian()
    )

    fit_mm <- brm(
      formula,
      data,
      seed = i,
      warmup = 500,
      iter = 1000,
      chains = 4, 
      cores = 4
    )

    evaluation$mm <- summarise(evaluation$mm, evaluate(fit_mm, beta, tau, corr, random_effects, outcomes), n_sim)
  }
  return(evaluation)
}

# Shared parameters
covariate_mean <- colMeans(epilepsy[c("zAge", "zBase")])
covariate_cov <- cor(epilepsy[c("zAge", "zBase")])
beta <- c(1, 0.1)
tau <- c(0, -0.5)
corr <- 0
n_points <- 5
n_sim <- 10

# Simulation study: Impact of number of individuals
n_individuals_list <- c(100, 200, 300)
for (n_individuals in n_individuals_list) {
  print(paste("Simulating for", n_individuals, "individuals"))
  write.table(simulation(n_sim, n_individuals, n_points, corr, beta, tau, covariate_mean, covariate_cov), paste("ind", n_individuals, ".txt"), sep = "\t", quote = FALSE)
}

# Simulation study: Impact of number of points
n_individuals <- 200
n_points_list <- c(5, 15, 25)
for (n_points in n_points_list) {
  print(paste("Simulating for", n_points, "points"))
  write.table(simulation(n_sim, n_individuals, n_points, corr, beta, tau, covariate_mean, covariate_cov), paste("points", n_points, ".txt"), sep = "\t", quote = FALSE)
}