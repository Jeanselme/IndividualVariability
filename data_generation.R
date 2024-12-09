# Generative functions
## Sampling function: Sample observation per patient
sampling_func <- function(observations, sampling = 'regular') { 
  # observations: list of #points for each patient

  if (sampling == 'regular') {
    times <- unlist(lapply(observations, function(n) seq_len(n))) # Create increasing time for each patient
    patients <- unlist(lapply(seq_along(observations), function(i) rep(i, observations[i]))) # Repeat id for each patient
    time <- data.frame(id = patients, time = times)
  }
  return(time)
}

## Generate covariates function at baseline
covariates_func <- function(sampling, n_individuals, observations, mean, sigma, columns) { 
  # sampling is the result of sampling_func
  # n_individuals: number of individuals
  # observations: list of #points for each patient
  # mean: mean of the covariates
  # sigma: covariance matrix of the covariates

  covariates <- mvrnorm(n = n_individuals, mu = mean, Sigma = sigma)

  # Propagate the covariates over time and make age change
  covariates <- covariates[unlist(lapply(seq_along(observations), function(i) rep(i, observations[i]))), ]
  covariates <- as.data.frame(covariates, col.names = columns)
  covariates$age <- scale(covariates$age + sampling$time)

  return(covariates)
}

## Random effects function (normally distributed) 
random_effects_func <- function(n_individuals, sds, corr_effects = 0, student = FALSE) { 
  # sds - Standard deviations for the different random effects
  # corr_effects - Correlation between the random effects
  # student - If TRUE, the random effects are drawn from a t-distribution

  # Create covariance for the random effects
  cov_matrix_vector <- c(
    sds[1]^2, corr_effects * sds[1] * sds[2], corr_effects * sds[1] * sds[3], corr_effects * sds[1] * sds[4],
    corr_effects * sds[1] * sds[2], sds[2]^2, corr_effects * sds[2] * sds[3], corr_effects * sds[2] * sds[4],
    corr_effects * sds[1] * sds[3], corr_effects * sds[2] * sds[3], sds[3]^2, corr_effects * sds[3] * sds[4],
    corr_effects * sds[1] * sds[4], corr_effects * sds[2] * sds[4], corr_effects * sds[3] * sds[4], sds[4]^2
  )

  cov_matrix <- matrix(cov_matrix_vector, nrow = 4, byrow = TRUE)
  if (student) {
    random_effects <- rmvt(mu = c(0, 0, 0, 0), n = n_individuals, sigma = cov_matrix, df = 3)
  } else {
    random_effects <- mvrnorm(n = n_individuals, mu = c(0, 0, 0, 0), Sigma = cov_matrix)
  }
  return(random_effects)
}



## Outcomes function  
outcomes_func <- function(covariates, random_effects, sample_times, beta, tau, time_slope = FALSE, sinus = FALSE) {    
  # covariates: data frame from covariates_func
  # random_effects: matrix from random_effects_func
  # sample_times: data frame from sampling_func
  # beta: fixed effects
  # tau: fixed effects for the variance
  # time_slope: if TRUE, the variance is time-dependent
  # sinus: if TRUE, the outcome is sinusoidal on age (not variance)

  # Random noise with the random effects
  covariates_out <- covariates
  if (sinus) {
    covariates_out$age <- sin(covariates$age)
  }

  if (time_slope) {
    omega <- exp(as.matrix(covariates) %*% tau + random_effects[sample_times$id, 2] + covariates$age * random_effects[sample_times$id, 3])
  } else {
    omega <- exp(as.matrix(covariates) %*% tau + random_effects[sample_times$id, 2])
  }
  noise <- rnorm(nrow(covariates), mean = 0, sd = omega)

  if (time_slope) {
    outcomes <- as.matrix(covariates_out) %*% beta + random_effects[sample_times$id, 1] + covariates_out$age * random_effects[sample_times$id, 4] + noise
  } else {
    outcomes <- as.matrix(covariates_out) %*% beta + random_effects[sample_times$id, 1] + noise
  }

  return(outcomes)
}

