library(brms)
library(MASS)
library(matrixStats)
library(rstan)
library(rlist)
library(survival)
library(dplyr)
library(ggplot2)
library(mvnfast)

source('data_generation.R')

# Create simulation study
simulation <- function(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, covariate_cov, time_slope = FALSE, sds = c(1, 0.5, 0.25, 0.25), student = FALSE, sinus = FALSE) {
  for (i in 1:n_sim) {
    folder_i <- paste0(path, i, '/')
    if (!dir.exists(folder_i)) {
      # Create folder for simulation
      dir.create(folder_i, recursive = TRUE)
    }

    # Fix seed
    set.seed(i)

    # Generate data
    observations <- sample(1:(2*n_points - 1), n_individuals, replace=TRUE) + 1 # Add one point for evaluation
    sample_times <- sampling_func(observations) 
    covariates <- covariates_func(sample_times, n_individuals, observations, covariate_mean, covariate_cov, columns)
    random_effects <- random_effects_func(n_individuals, sds, corr_effects = corr, student = student)
    outcomes <- outcomes_func(covariates, random_effects, sample_times, beta, tau, time_slope, sinus)   

    last_time_indices <- sample_times %>%
      group_by(id) %>%
      summarize(index_last = (time == max(time)))
    last_time_indices <- last_time_indices$index_last

    data <- data.frame(sample_times, covariates, outcomes = outcomes)

    # Fit models
    for (method in names(formulas)) {
      # Run if file does not exist
      if (!file.exists(paste0(folder_i, method, '.rds'))) {
        fit <- brm(formulas[[method]], data[!last_time_indices,], seed = i, warmup = 10, iter = 20, chains = 1, cores = 4)
        # Save model
        saveRDS(fit, file = paste0(folder_i, method, '.rds'))

        # Extract all columns of interest and save
        write.csv(posterior_summary(fit), paste0(folder_i, method, '.csv'), row.names = TRUE) 
      } else {
        print(paste0('Skipping ', folder_i, method))
      }
    }
  }
}


# Default parameters
## Data
columns <- c("age", "albumin", "trig", "platelet")
pbc_scale <- scale(pbc[columns])
covariate_mean <- colMeans(pbc_scale, na.rm = TRUE)
covariate_cov <- cov(pbc_scale, use="complete.obs")

## True model (CANNOT CHANGE 0. without changing experiment)
beta <- c(0.5, 0.5, 0., 0.)
tau <- c(0.8, 0, 0.8, 0.)
corr <- 0

## Population
n_individuals <- 200
n_points <- 15

## Simulation
n_sim <- 2

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
      sigma ~ age + trig + (1|id), 
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
  n_individuals_list <- c(100, 300, 500, 1000)
  for (n_individuals_exp in n_individuals_list) {
    print(paste("Simulating for", n_individuals_exp, "individuals"))
    path = paste0("results/Individuals/", n_individuals_exp, '/')
    simulation(path, formulas, n_sim, n_individuals_exp, n_points, corr, columns, beta, tau, covariate_mean, covariate_cov)
  }
}

# Simulation study: Impact of number of points
if ((run == -1)|(run == 2)) {
  n_points_list <- c(5, 10, 20)
  for (n_points_exp in n_points_list) {
    print(paste("Simulating for", n_points_exp, "points"))
    path = paste0("results/Points/", n_points_exp, '/')
    simulation(path, formulas, n_sim, n_individuals, n_points_exp, corr, columns, beta, tau, covariate_mean, covariate_cov)
  }
}

# STUDY 2
### Formula for MELSM comparison
formulas = list(
    correct = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ age + trig + (1|id), 
      family = gaussian()
    ),
    sigma = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ trig + (1|id), 
      family = gaussian()
    ),
    outcomes = bf(
      outcomes ~ albumin + (1|id),
      sigma ~ age + trig + (1|id), 
      family = gaussian()
    ),
    all = bf(
      outcomes ~ age + albumin + trig + platelet + (1|id),
      sigma ~ age + albumin + trig + platelet + (1|id), 
      family = gaussian()
    ),
    nore = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ age + trig, 
      family = gaussian()
    ),
    sigmanore = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ trig, 
      family = gaussian()
    ),
    mm = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ 1,
      family = gaussian()
    )
  )

if ((run == -1)|(run == 3)) {
  print("Simulating when no correlation")
  path = "results/Uncorrelated Data/"
  identity = diag(4)
  simulation(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, identity)
}

if ((run == -1)|(run == 4)) {
  print("Simulating for MELSM misspecification")
  path = "results/Misspecification/"
  simulation(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, covariate_cov)
}

# STUDY 3
### Formula for time comparison
formulas = list(
    correct = bf(
      outcomes ~ age + albumin + (1 + age|id),
      sigma ~ age + trig + (1 + age|id), 
      family = gaussian()
    ),
    melsm_notimeomega = bf(
      outcomes ~ age + albumin + (1 + age|id),
      sigma ~ age + trig + (1|id), 
      family = gaussian()
    ),
    melsm_notime = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ age + trig + (1|id), 
      family = gaussian()
    )
  )

if ((run == -1)|(run == 5)) {
  print("Simulating for time dependent random effects")
  path = "results/Time/"
  simulation(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, covariate_cov, time_slope = TRUE)
}

# STUDY 4
### Random effect family
formulas = list(
    correct = bf(
      outcomes ~ age + albumin + (1|gr(id, dist='student')),
      sigma ~ age + trig + (1|gr(id, dist='student')), 
      family = gaussian()
    ),
    gaussian = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ age + trig + (1|id), 
      family = gaussian()
    )
  )

if ((run == -1)|(run == 6)) {
  print("Simulating for time dependent random effects")
  path = "results/Random effect/"
  simulation(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, covariate_cov, student = TRUE)
}

# STUDY 5
### Correlation
formulas = list(
    correct = bf(
      outcomes ~ age + albumin + (1|I|id),
      sigma ~ age + trig + (1|I|id), 
      family = gaussian()
    ),
    incorrect = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ age + trig + (1|id), 
      family = gaussian()
    )
  )

if ((run == -1)|(run == 7)) {
  n_corr_list <- c(0.5, 0.9)
  for (n_corr_exp in n_corr_list) {
    print("Simulating for correlated random effects")
    path = paste0("results/Correlated Random Effect/", n_corr_exp, '/')
    simulation(path, formulas, n_sim, n_individuals, n_points, n_corr_exp, columns, beta, tau, covariate_mean, covariate_cov)
  }
}

# STUDY 6
### Sinus
formulas = list(
    correct = bf(
      outcomes ~ sin(age) + albumin + (1|id),
      sigma ~ age + trig + (1|id), 
      family = gaussian()
    ),
    nonsinus = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ age + trig + (1|id), 
      family = gaussian()
    )
  )

if ((run == -1)|(run == 8)) {
  print("Simulating for sinus age")
  path = "results/Sinus/"
  simulation(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, covariate_cov, sinus = TRUE)
}