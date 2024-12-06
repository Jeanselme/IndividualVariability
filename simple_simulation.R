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
source('evaluation.R')

# Create simulation study
simulation <- function(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, covariate_cov, time_slope = FALSE, sds = c(2, 1, 0.5), student = FALSE, sinus = FALSE) {
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
    observations <- sample(1:(2*n_points - 1), n_individuals, replace=TRUE) + 1 # Add one point for evaluation
    sample_times <- sampling_func(observations) 
    covariates <- covariates_func(sample_times, n_individuals, observations, covariate_mean, covariate_cov, columns)
    random_effects <- random_effects_func(n_individuals, sds, corr_effects = corr, student = student)
    if (sinus) {
      covariates <- as.data.frame(covariates)
      covariates$sinage <- sin(covariates$age)
      outcomes <- outcomes_func(as.matrix(covariates[, c("sinage",  "albumin", "trig", "platelet")]), random_effects, sample_times, beta, tau, time_slope)
    } else {
      outcomes <- outcomes_func(covariates, random_effects, sample_times, beta, tau, time_slope)
    }    

    last_time_indices <- sample_times %>%
      group_by(id) %>%
      summarize(index_last = (time == max(time)))
    last_time_indices <- last_time_indices$index_last

    data <- data.frame(sample_times, covariates, outcomes = outcomes)

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
n_sim <- 100

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
    path = paste0("results/individuals", n_individuals_exp)
    simulation(path, formulas, n_sim, n_individuals_exp, n_points, corr, columns, beta, tau, covariate_mean, covariate_cov)
  }
}

# Simulation study: Impact of number of points
if ((run == -1)|(run == 2)) {
  n_points_list <- c(5, 10, 20)
  for (n_points_exp in n_points_list) {
    print(paste("Simulating for", n_points_exp, "points"))
    path = paste0("results/points", n_points_exp)
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
      sigma ~ age + albumin + (1|id), 
      family = gaussian()
    ),
    outcomes = bf(
      outcomes ~ age + trig + (1|id),
      sigma ~ age + trig + (1|id), 
      family = gaussian()
    ),
    all = bf(
      outcomes ~ age + albumin + trig + platelet + (1|id),
      sigma ~ age + albumin + trig + platelet + (1|id), 
      family = gaussian()
    ),
    nore = bf(
      outcomes ~ age + albumin,
      sigma ~ age + trig, 
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

if ((run == -1)|(run == 3)) {
  print("Simulating when no correlation")
  path = "results/nocorr"
  identity = diag(4)
  simulation(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, identity)
}

if ((run == -1)|(run == 4)) {
  print("Simulating for MELSM misspecification")
  path = "results/misspecification"
  simulation(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, covariate_cov)
}

# STUDY 3
### Formula for time comparison
formulas = list(
    correct = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ age + trig + (1 + age|id), 
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
  path = "results/time"
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
  path = "results/random_effects"
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
  print("Simulating for correlated random effects")
  path = "results/corr_random_effects"
  simulation(path, formulas, n_sim, n_individuals, n_points, 0.5, columns, beta, tau, covariate_mean, covariate_cov)
}

# STUDY 6
### Sinus
formulas = list(
    correct = bf(
      outcomes ~ sin(age) + albumin + (1|id),
      sigma ~ sin(age) + trig + (1|id), 
      family = gaussian()
    ),
    incorrect = bf(
      outcomes ~ age + albumin + (1|id),
      sigma ~ age + trig + (1|id), 
      family = gaussian()
    )
  )

if ((run == -1)|(run == 8)) {
  print("Simulating for sinus age")
  path = "results/sinus"
  simulation(path, formulas, n_sim, n_individuals, n_points, corr, columns, beta, tau, covariate_mean, covariate_cov, sinus = TRUE)
}