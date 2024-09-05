library(brms)
library(MASS)
library(matrixStats)
library(rstan)
library(rlist)
library(survival)
library(dplyr)
library(ggplot2)
library(mvnfast)

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

age_diff_func <- function(sample_times, covariates) {
  sample_times <- as.data.frame(sample_times)
  covariates <- as.data.frame(covariates)

  combined_data <- bind_cols(sample_times, covariates) %>%
    group_by(id) %>%
    mutate(time = age - first(age)) %>%
    select(id, time)  # Select only the id and time columns
  
  return(combined_data)
}

## Random effects function (normally distributed) 
random_effects_func <- function(n_individuals, sds = c(2, 1, 0.5), corr_effects = 0, student = FALSE) { 
  # Create covariance for the random effects
  sd1 <- sds[1]
  sd2 <- sds[2]
  sd3 <- sds[3]
  cov_matrix <- matrix(c(sd1^2, corr_effects*sd1*sd2, corr_effects*sd1*sd3, corr_effects*sd1*sd2, sd2^2, corr_effects*sd2*sd3, corr_effects*sd1*sd3, corr_effects*sd2*sd3, sd3^2), nrow = 3)
  if (student) {
    random_effects <- rmvt(mu = c(0, 0, 0), n = n_individuals, sigma = cov_matrix, df = 15)
  } else {
    random_effects <- mvrnorm(n = n_individuals, mu = c(0, 0, 0), Sigma = cov_matrix)
  }
  return(random_effects)
}

## Outcomes function  
outcomes_func <- function(covariates, random_effects, sample_times, beta, tau, time_slope = FALSE) { 
  # Random noise with the random effects
  if (time_slope) {
    omega <- exp(covariates %*% tau + sample_times$time + random_effects[sample_times$id, 2] + sample_times$time * random_effects[sample_times$id, 3])
  } else {
    omega <- exp(covariates %*% tau + random_effects[sample_times$id, 2])
  }
  noise <- rnorm(nrow(covariates), mean = 0, sd = omega)
  
  # Arbitrary choice of outcomes
  outcomes <- covariates %*% beta + random_effects[sample_times$id, 1] + noise

  return(outcomes)
}



# Compute error for each model
## Function rmse
rerr <- function(truth, estimate) {
  return(mean(abs(truth - estimate)))
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
      relative[[col]] <- abs_diffs[[col]] / abs(beta[i])
      in_bounds[[col]] <- as.numeric(beta[i] >= lower_bound & beta[i] <= upper_bound)

    }, error = function(e) {
      abs_diffs[[col]] <<- 0
      relative[[col]] <<- 0
      in_bounds[[col]] <<- 0
    })
  }
  return(list(error=unlist(abs_diffs), relative=unlist(relative), coverage=unlist(in_bounds)))
}

## Function evaluate (silent if error as some methods do not compute all)
evaluate <- function(fit, newdata, columns, beta, tau, sds, corr, random_effects, indices) {
  errors <- list(time = sum(rstan::get_elapsed_time(fit$fit)))

  # RMSE beta
  beta_estimate = error_coeff(fit, columns, beta)
  errors$beta_out <- beta_estimate$error
  errors$beta_out_coverage <- beta_estimate$coverage

  # RMSE tau
  tau_estimate = error_coeff(fit, columns, tau, 'b_sigma_')
  errors$beta_omega <- tau_estimate$error
  errors$beta_omega_coverage <- tau_estimate$coverage

  # RMSE corr
  try(errors$corr <- rerr(corr, posterior_summary(fit, variable='cor_id__Intercept__sigma_Intercept')[, 'Estimate']), silent = TRUE)

  # RMSE random effects
  try(errors$re_out <- rerr(random_effects[, 1],  posterior_summary(fit, variable='r_id')[, 'Estimate']), silent = TRUE)
  try(errors$re_omega_0 <- rerr(random_effects[, 2], posterior_summary(fit, pars='r_id__sigma.*Intercept.*')[, 'Estimate']), silent = TRUE)
  try(errors$re_omega_1 <- rerr(random_effects[, 3], posterior_summary(fit, pars='r_id__sigma.*time.*')[, 'Estimate']), silent = TRUE)

  sd1 <- sds[1]
  sd2 <- sds[2]
  sd3 <- sds[3]
  try(errors$sd_re_out <- rerr(sd1, posterior_summary(fit, variable='sd_id__Intercept')[, 'Estimate']), silent = TRUE)
  try(errors$sd_re_omega_0 <- rerr(sd2, posterior_summary(fit, variable='sd_id__sigma_Intercept')[, 'Estimate']), silent = TRUE)
  try(errors$sd_re_omega_1 <- rerr(sd3, posterior_summary(fit, variable='sd_id__sigma_time')[, 'Estimate']), silent = TRUE)

  # RMSE outcomes
  errors$last <- rerr(fitted(fit, newdata = newdata[indices,])[, 'Estimate'], newdata$outcomes[indices])
  errors$average <- rerr(fitted(fit)[, 'Estimate'], newdata$outcomes[!indices])

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

summarise_perf <- function(evaluation) {
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

spaghetti_plot <- function(data, num_ids = 5) {
  # Randomly select a few unique ids
  selected_ids <- sample(unique(data$id), num_ids)
  
  # Filter data to include only the selected ids
  filtered_data <- data %>%
    filter(id %in% selected_ids)
  
  # Create the spaghetti plot
  p <- ggplot(filtered_data, aes(x = time, y = outcomes, group = id, color = as.factor(id))) +
    geom_line() +
    geom_point() +
    labs(title = "Spaghetti Plot of Outcomes Over Time",
         x = "Time",
         y = "Outcome",
         color = "ID") +
    theme_minimal()

  ggsave('test.png', plot = p, width = 8, height = 6)
}

average_plot <- function(data, name) {
  p <- ggplot(data, aes(x = time, y = outcomes)) +
    geom_smooth(method = "loess", color = "blue", se = FALSE) +
    labs(title = "Smooth Approximation of Outcomes Over Time",
         x = "Time",
         y = "Outcome") +
    theme_minimal()

  ggsave(name, plot = p, width = 8, height = 6)
}

# Create simulation study
simulation <- function(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov, time_slope = FALSE, sds = c(2, 1, 0.5), student = FALSE) {
  json = paste0(path, '.json')
  if (file.exists(json)) {
    # Load the file
    evaluation = list.load(json)
    print("File loaded successfully.")
    print(paste("Already", evaluation$start, "simulations done."))
  } else {
    # Initialise evaluation
    evaluation = list(start = 1)
    for (name in names(formulas)) {
      evaluation[[name]] = list()

    }
    print("File does not exist.")
  }  
  if (evaluation$start > n_sim) {
    return()
  }
  for (i in evaluation$start:n_sim) {
    # Fix seed
    set.seed(i)

    # Generate data
    sample_times <- sampling_func(n_individuals, n_points + 1) # Add one point for evaluation
    covariates <- covariates_func(sample_times, covariate_mean, time_dependent, covariate_cov)
    sample_times <- age_diff_func(sample_times, covariates) # Update times to respect age difference
    random_effects <- random_effects_func(n_individuals, sds, corr_effects = corr, student = student)
    outcomes <- outcomes_func(covariates, random_effects, sample_times, beta, tau, time_slope)

    last_time_indices <- sample_times %>%
      group_by(id) %>%
      summarize(index_last = (time == max(time)))
    last_time_indices <- last_time_indices$index_last

    data <- data.frame(sample_times, covariates, outcomes = outcomes)
    # spaghetti_plot(data, 10)

    # Fit models
    for (method in names(formulas)) {
      fit <- brm(formulas[[method]], data[!last_time_indices,], seed = i, warmup = 1000, iter = 2000, chains = 4, cores = 4)
      evaluation[[method]] <- append(evaluation[[method]], evaluate(fit, data, columns, beta, tau, sds, corr, random_effects, last_time_indices))
    }
    
    # Save
    evaluation$start = evaluation$start + 1
    list.save(evaluation, file = json)
    list.save(summarise_perf(evaluation), file = paste0(path, '_summary.json'))
  }
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


# STUDY 1 
### Formula for MM comparison
formulas = list(
    correct = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ trig + platelet + (1|id), 
      family = gaussian()
    ),
    mm = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ 1,
      family = gaussian()
    )
  )

# Simulation study: Impact of number of individuals
if ((run == -1)|(run == 1)) {
  n_individuals_list <- c(100, 200, 300)
  for (n_individuals_exp in n_individuals_list) {
    print(paste("Simulating for", n_individuals_exp, "individuals"))
    path = paste0("results/individuals", n_individuals_exp)
    simulation(path, formulas, n_sim, n_individuals_exp, n_points, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov)
  }
}

# Simulation study: Impact of number of points
if ((run == -1)|(run == 2)) {
  n_points_list <- c(5, 10, 20)
  for (n_points_exp in n_points_list) {
    print(paste("Simulating for", n_points_exp, "points"))
    path = paste0("results/points", n_points_exp)
    simulation(path, formulas, n_sim, n_individuals, n_points_exp, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov)
  }
}

if ((run == -1)|(run == 3)) {
  print("Simulating when no correlation")
  path = "results/nocorr"
  identity = diag(4)
  simulation(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, time_dependent, identity)
}

# # Simulation study: Impact of number of corr
# if ((run == -1)|(run == 3)) {
#   rho_list <- c(-0.5, -0.25, 0.25, 0.5)
#   for (corr_exp in rho_list) {
#     print(paste("Simulating for", corr_exp, "correlation"))
#     path = paste0("results/corr", corr_exp)
#     simulation(path, formulas, n_sim, n_individuals, n_points, corr_exp, columns, beta, tau, covariate_mean, time_dependent, covariate_cov)
#   }
# }


# STUDY 2
### Formula for MELSM comparison
formulas = list(
    correct = bf(
      outcomes ~ age + albumin + 1 + (1|id),
      sigma ~ trig + platelet + (1|id), 
      family = gaussian()
    ),
    sigma = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ age + albumin + (1|id), 
      family = gaussian()
    ),
    outcomes = bf(
      outcomes ~ trig + platelet + (1|id),
      sigma ~ trig + platelet + (1|id), 
      family = gaussian()
    ),
    all = bf(
      outcomes ~ age + albumin + trig + platelet + (1|id),
      sigma ~ age + albumin + trig + platelet + (1|id), 
      family = gaussian()
    ),
    nore = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ trig + platelet, 
      family = gaussian()
    ),
    sigmanore = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ age + albumin, 
      family = gaussian()
    ),
    mm = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ 1,
      family = gaussian()
    )
  )

if ((run == -1)|(run == 4)) {
  print("Simulating for MELSM misspecification")
  path = "results/misspecification"
  simulation(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov)
}

# STUDY 3
### Formula for time comparison
formulas = list(
    correct = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ trig + platelet + time + (1 + time|id), 
      family = gaussian()
    ),
    melsm_notime = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ trig + platelet + (1|id), 
      family = gaussian()
    )
  )

if ((run == -1)|(run == 5)) {
  print("Simulating for time dependent random effects")
  path = "results/time"
  simulation(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov, time_slope = TRUE)
}

# STUDY 4
### Random effect family
formulas = list(
    correct = bf(
      outcomes ~ age + albumin + (1|gr(id, dist='student')),
      sigma ~ trig + platelet + (1|gr(id, dist='student')), 
      family = gaussian()
    ),
    gaussian = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ trig + platelet + (1|id), 
      family = gaussian()
    )
  )

if ((run == -1)|(run == 6)) {
  print("Simulating for time dependent random effects")
  path = "results/random_effects"
  simulation(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, time_dependent, covariate_cov, student = TRUE)
}

# STUDY 5
### Correlation
formulas = list(
    correct = bf(
      outcomes ~ age + albumin + (1|I|id),
      sigma ~ trig + platelet + (1|I|id), 
      family = gaussian()
    ),
    incorrect = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ trig + platelet + (1|id), 
      family = gaussian()
    )
  )

if ((run == -1)|(run == 7)) {
  print("Simulating for correlated random effects")
  path = "results/corr_random_effects"
  simulation(path, formulas, n_sim, n_individuals, n_points, 0.5, columns, beta, tau, covariate_mean, time_dependent, covariate_cov, student = TRUE)
}