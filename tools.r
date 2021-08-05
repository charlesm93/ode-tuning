
library(wesanderson)

# TODO: Create a saved_last_iter function which works for every model.
# A function to extract the final sample of warmup (first sample of
# "sampling" phases) and save it.
# Designed for the Michaelis-Menten PK model for one individual.
save_last_iter <- function(samples, nChains, sample_name, parm_index) {
  saved_file_name <- paste0("init/", sample_name, 1:nChains, ".json")
  for (i in 1:nChains) {
    init_saved <- list(ka = samples[1, i, parm_index[1]],
                       V = samples[1, i, parm_index[2]],
                       Vm = samples[1, i, parm_index[3]],
                       Km = samples[1, i, parm_index[4]],
                       sigma = samples[1, i, parm_index[5]])

    write_stan_json(init_saved, saved_file_name[i])
  }

  saved_file_name
}

# Same as above, but for population model with a centered parameterization.
# The first sampling iteration is treated as the last warmup iteration.
save_last_iter_pop <- function(samples, nChains, sample_name, parm_index,
                              nIIV = 4, n_patients = 10) {
  saved_file_name <- paste0("init/", sample_name, 1:nChains, ".json")
  for (i in 1:nChains) {
    init_saved <- list(ka_pop = samples[1, i, parm_index[1]],
                       V_pop = samples[1, i, parm_index[2]],
                       Vm_pop = samples[1, i, parm_index[3]],
                       Km_pop = samples[1, i, parm_index[4]],
                       sigma = samples[1, i, parm_index[5]])

    omega <- rep(NA, nIIV)
    for (j in 1:nIIV) {
      omega[j] <- samples[1, i, parm_index[5 + j]]
    }
    init_saved$omega <- omega

    theta <- matrix(NA, nrow = nIIV, ncol = n_patients)
    for (j in 1:n_patients) {
      theta[, j] <- as.vector(samples[1, i, parm_index[5 + nIIV +
                                            ((j - 1) * nIIV + 1):(j * nIIV)]])
    }
    init_saved$theta <- theta
    write_stan_json(init_saved, saved_file_name[i])
  }

  saved_file_name
}


# A function to run the warmup and save the fit for each phase.
# save_function should either be save_last_iter or save_last_iter_pop.
fit_warmup_phases <- function(mod, model_name, data, save_function,
                              parm_index = 2:6,
                              init = NULL,
                              iter_warmup = 500, 
                              init_buffer = 75,
                              term_buffer = 50, nChains = 4, seed = 123,
                              metric = "diag_e",
                              adapt_delta = 0.8,
                              ...) {
  fit_w1 <- mod$sample(data = data, init = init,
                       chains = nChains, parallel_chains = parallel_chains,
                       iter_warmup = init_buffer,
                       iter_sampling = 1,
                       init_buffer = init_buffer,
                       term_buffer = 0, window = 0,
                       metric = metric,
                       seed = seed, refresh = 10, adapt_delta = adapt_delta)
  
  samples_w1 <- fit_w1$draws()
  init_files <- save_function(samples_w1, nChains,
                              paste0(model_name, "_w1."),
                              parm_index, ...)

  if (metric == "diag_e") {
    metric_is_matrix <- FALSE
  } else {
    metric_is_matrix <- TRUE
  }

  window_size <- iter_warmup - (init_buffer + term_buffer)
  fit_w2 <- mod$sample(data = data, init = init_files,
                       chains = nChains, parallel_chains = parallel_chains,
                       iter_sampling = 1,
                       iter_warmup = window_size,
                       init_buffer = 0, term_buffer = 0,
                       metric = metric,
                       seed = seed,
                       step_size = fit_w1$metadata()$step_size_adaptation,
                       inv_metric = fit_w1$inv_metric(matrix = metric_is_matrix),
                       refresh = 10, adapt_delta = adapt_delta)

  samples_w2 <- fit_w2$draws()
  init_files2 <- save_function(samples_w2, nChains, paste0(model_name, "_w2."),
                               parm_index, ...)
  
  fit_w3 <- mod$sample(data = data, init = init_files2,
                       chains = nChains, parallel_chains = parallel_chains,
                       iter_sampling = 1,
                       iter_warmup = term_buffer, window = 0,
                       init_buffer = 0, term_buffer = term_buffer,
                       metric = metric,
                       seed = seed,
                       step_size = fit_w2$metadata()$step_size_adaptation,
                       inv_metric = fit_w2$inv_metric(matrix = metric_is_matrix),
                       refresh = 10, adapt_delta = adapt_delta)
  
  # Return
  list(fit_w1 = fit_w1, fit_w2 = fit_w2, fit_w3 = fit_w3)
}


ess_summary <- function(fit, parms, nChains, time) {
  n_parm <- length(parms)
  draws <- fit$draws(variables = parms)
  chain_ess <- matrix(NA, n_parm, nChains)
  chain_eff <- matrix(NA, n_parm, nChains)
  for (i in 1:nChains) {
    chain_draw <- as_draws_df(draws[, i, ])
    chain_ess[, i] <- dplyr::pull(summarize_draws(chain_draw, 
                                                  measure = c("ess_bulk")), 2) 
    chain_eff[, i] <- chain_ess[, i] / time[i]
  }
  
  mean_ess <- rep(NA, n_parm)
  sd_ess <- rep(NA, n_parm)
  for (i in 1:n_parm) {
    mean_ess[i] <- mean(chain_eff[i, ])
    sd_ess[i] <- sd(chain_eff[i, ])
  }
  
  list(chain_ess = chain_ess, chain_eff = chain_eff,
       mean_ess = mean_ess, sd_ess = sd_ess)
}


# Use this function to plot warmup and sampling time, if these are
# all included in the fit object.
plot_run_time <- function(fit_object, time_data) {
  time_data <- fit_object$time()$chains
  run_times <- c(time_data$warmup, time_data$sampling)
  chain_id <- rep(time_data$chain_id, 2)
  phase <- rep(c("warmup", "sampling"), each = nChains)
  
  plot_data <- data.frame(run_times = run_times,
                          chain_id = chain_id,
                          phase = phase)
  
  plot <- ggplot(data = plot_data, 
                 aes(x = chain_id, y = run_times, fill = phase)) +
    geom_bar(stat = "identity", position = "stack", width = 0.3) +
    coord_flip() + theme_bw() +
    ylab("Run time (s)") + xlab("Chain ID") +
    theme(text = element_text(size = 16))
  plot
}

plot_phase_time <- function(run_times, chain_id, phase,
                            integrator = NULL) {
  plot_data <- data.frame(run_times = run_times,
                          chain_id = chain_id,
                          phase = phase)
  
  if (!is.null(integrator)) plot_data$integrator <- integrator

  plot <- ggplot(data = plot_data, 
                 aes(x = chain_id, y = run_times, fill = phase)) +
    geom_bar(stat = "identity", position = "stack", width = 0.3) +
    coord_flip() + theme_bw() +
    ylab("Run time (s)") + xlab("Chain ID") +
    theme(text = element_text(size = 16)) +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "coral"))
  
  if (!is.null(integrator)) {
    plot <- plot + facet_wrap(~integrator)
  }

  plot
}

