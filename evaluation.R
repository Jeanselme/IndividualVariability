# Compute error for each model
## Function rmse
rerr <- function(truth, estimate) {
  return(mean(abs(truth - estimate)))
}

coverage <- function(truth, summary_stats) {
  lower_bound <- summary_stats[, 'Q2.5']
  upper_bound <- summary_stats[, 'Q97.5']
  return(as.numeric(truth >= lower_bound & truth <= upper_bound))
}

error_coeff <- function(fit, columns, beta, prefix = 'b_') {
  abs_diffs <- list()
  relative <- list()
  in_bounds <- list()
  mean <- list()
  for (i in seq_along(beta)) {
    col <- columns[i]
    tryCatch({
      summary_stats <- tryCatch({
        posterior_summary(fit, variable = paste0(prefix, col))
      }, error = function(e) {return(posterior_summary(fit, variable = paste0(prefix, 'sin', col)))})
      
      estimate <- summary_stats[, 'Estimate']

      mean[[col]] <- estimate
      abs_diffs[[col]] <- abs(beta[i] - estimate)
      relative[[col]] <- abs_diffs[[col]] / abs(beta[i])
      in_bounds[[col]] <- coverage(beta[i], summary_stats)
    }, error = function(e) {
      mean[[col]] <<- NA
      abs_diffs[[col]] <<- NA
      relative[[col]] <<- NA
      in_bounds[[col]] <<- NA
    })
  }
  return(list(mean=unlist(mean), error=unlist(abs_diffs), relative=unlist(relative), coverage=unlist(in_bounds)))
}

## Function evaluate (silent if error as some methods do not compute all)
evaluate <- function(fit, newdata, columns, beta, tau, sds, corr, random_effects, indices) {
  errors <- list(time = sum(rstan::get_elapsed_time(fit$fit)))

  # RMSE beta
  beta_estimate = error_coeff(fit, columns, beta)
  errors$beta <- beta_estimate$mean
  errors$beta_out_coverage <- beta_estimate$coverage

  # RMSE tau
  tau_estimate = error_coeff(fit, columns, tau, 'b_sigma_')
  errors$beta_omega <- tau_estimate$mean
  errors$beta_omega_coverage <- tau_estimate$coverage

  # RMSE corr
  try(errors$corr <- rerr(corr, posterior_summary(fit, variable='cor_id__Intercept__sigma_Intercept')[, 'Estimate']), silent = TRUE)

  # RMSE random effects
  try(errors$sd_re_out_0_error <- rerr(sds[1], posterior_summary(fit, variable='sd_id__Intercept')[, 'Estimate']), silent = TRUE)
  try(errors$sd_re_out_0_coverage <- coverage(sds[1], posterior_summary(fit, variable='sd_id__Intercept')), silent = TRUE)

  try(errors$sd_re_out_1_error <- rerr(sds[4], posterior_summary(fit, variable='sd_id__age')[, 'Estimate']), silent = TRUE)
  try(errors$sd_re_out_1_coverage <- coverage(sds[4], posterior_summary(fit, variable='sd_id__age')), silent = TRUE)

  try(errors$sd_re_omega_0_error <- rerr(sds[2], posterior_summary(fit, variable='sd_id__sigma_Intercept')[, 'Estimate']), silent = TRUE)
  try(errors$sd_re_omega_0_coverage <- coverage(sds[2], posterior_summary(fit, variable='sd_id__sigma_Intercept')), silent = TRUE)

  try(errors$sd_re_omega_1_error <- rerr(sds[3], posterior_summary(fit, variable='sd_id__sigma_age')[, 'Estimate']), silent = TRUE)
  try(errors$sd_re_omega_1_coverage <- coverage(sds[3], posterior_summary(fit, variable='sd_id__sigma_age')), silent = TRUE)

  try(errors$sd_re_out_0 <- mean(posterior_summary(fit, variable='sd_id__Intercept')[, 'Estimate']), silent = TRUE)
  try(errors$sd_re_out_1 <- mean(posterior_summary(fit, variable='sd_id__age')[, 'Estimate']), silent = TRUE)
  try(errors$sd_re_omega_0 <- mean(posterior_summary(fit, variable='sd_id__sigma_Intercept')[, 'Estimate']), silent = TRUE)
  try(errors$sd_re_omega_1 <- mean(posterior_summary(fit, variable='sd_id__sigma_age')[, 'Estimate']), silent = TRUE)

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
        value=colMeans(evaluation[[model]][[quantity]], na.rm = TRUE),
        std=tryCatch({colSds(evaluation[[model]][[quantity]], na.rm = TRUE)}, error = function(e) {return(NA)})
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
