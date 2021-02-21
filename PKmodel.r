
rm(list = ls())
gc()
set.seed(1954)

setwd("~/Code/ode-tuning")
.libPaths("~/Rlib/")

library(ggplot2)
library(cmdstanr)
set_cmdstan_path("~/Rlib/cmdstan/")
library(rjson)
library(posterior)
library(deSolve)
library(bayesplot)

model_name <- "Michaelis_MentenPK"

#####################################################################
## Simulate data

simulate_data <- TRUE
if (simulate_data) {
  # simulate data
  # ka <- exp(rnorm(1, log(1), 1))
  # V <- exp(rnorm(1, log(35), 0.25))
  # Vm <- exp(rnorm(1, log(10), 0.25))
  # Km <- exp(rnorm(1, log(2.5), 1))
  # sigma <- abs(rnorm(1, 0, 0.1))
  # ka <- exp(rnorm(1, log(2.5), 1))
  # V <- exp(rnorm(1, log(35), 0.5))
  # Vm <- exp(rnorm(1, log(10), 0.5))
  # Km <- exp(rnorm(1, log(2.5), 1.5))
  # sigma <- abs(rnorm(1, 0, 0.1))
  ka <- 3
  V <- 27
  Vm <- 10
  Km <- 14
  sigma <- 0.1
  y0 <- 100

  parm <- c(ka, V, Vm, Km)
  state <- c(Y1 = y0, Y2 = 0)
  dydt <- function(t, state, parm) {
    with(as.list(c(state, parm)), {
      C = Y2 / V
      dY1 <- - ka * Y1
      dY2 <- ka * Y1 - Vm * C / (Km + C)

      list(c(dY1, dY2))
    })
  }

  times <- c(seq(0.5, 1.5, by = 0.25), seq(10, 100, by = 10))
  n <- length(times)

  mass <- ode(y = state, time = times, func = dydt, parms = parm)
  concentration <- mass[, 3] / V
  y <- rnorm(n, concentration, sigma)

  stan_data <- list(N = n,  y = y, t = times, y0 = c(y0, 0))
  write_stan_json(stan_data, "data/PKModel.data.json")
  
  p <- ggplot(data = data.frame(y = y, times = times),
              aes(y = y, x = times)) +
    geom_point() + theme_bw() + theme(text = element_text(size = 20))
  p
  # plot(x = mass[, 1], y = y)
}

stan_data_rk45 <- fromJSON(file = "data/PKModel.data.json")
stan_data_rk45$stiff_solver <- 0
stan_data_bdf <- stan_data_rk45
stan_data_bdf$stiff_solver <- 1

#####################################################################
## Build inits (we'll use the same inits for model fit)

init <- function() {
  list(
    ka = exp(rnorm(1, log(1.5), 3)),
    V = exp(rnorm(1, log(35), 0.5)),
    Vm = exp(rnorm(1, log(10), 0.5)),
    Km = exp(rnorm(1, log(2.5), 3)),
    sigma = abs(rnorm(1, 0, 1))
  )
}

nChains <- 8
parallel_chains <- 8
for (i in 1:nChains) {
  init0 <- init()
  write_stan_json(init0, paste0("init/init", i, ".json"))
}

#####################################################################
## Tuning parameters for MCMC + compile model

iter_warmup <- 500
iter_sampling <- 500

## Compile model
if (TRUE) {
  file <- paste0("model/", model_name, ".stan")
  mod <- cmdstan_model(file)
}

#####################################################################
## Fit model with rk45 solver.

run_model <- TRUE
saved_fit0_file <- paste0("output/", model_name, ".rk45.fit.RDS")
if (run_model) {
  fit0 <- mod$sample(
      data = stan_data_rk45, chains = nChains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      seed = 123, adapt_delta = 0.8,
      init = paste0("init/init", 1:nChains, ".json"))

  fit0$cmdstan_diagnose()
  
  fit0$save_object(saved_fit0_file)
}

fit0 <- readRDS(saved_fit0_file)
fit0$time()
fit0$metadata()$step_size_adaptation

#####################################################################
## Fit model with bdf solver

saved_fit1_file <- paste0("output/", model_name, ".bdf.fit.RDS")
if (run_model) {
  fit <- mod$sample(
    data = stan_data_bdf, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    seed = 123, adapt_delta = 0.8,
    init = paste0("init/init", 1:nChains, ".json"))

  fit$cmdstan_diagnose()
  
  fit$save_object(saved_fit1_file)
}
  
fit <- readRDS(saved_fit1_file)
fit$time()

#####################################################################
## Run sampling with rk45, using tuning parameters warmed up with bdf.

mass_matrix <- fit$inv_metric()
step_size <- fit$metadata()$step_size_adaptation

# stan_data2 <- stan_data
# stan_data2$stiff_solver <- 0

# Create init files, using warmup samples and save them.
# CHECK -- the draws should be ordered
samples <- fit$draws()
parm_index <- 2:6
for (i in 1:nChains) {
  init_saved <- list(ka = samples[1, i, parm_index[1]],
                     V = samples[1, i, parm_index[2]],
                     Vm = samples[1, i, parm_index[3]],
                     Km = samples[1, i, parm_index[4]],
                     sigma = samples[1, i, parm_index[5]])

  write_stan_json(init_saved, paste0("init/warm_init", i, ".json"))
}

saved_fit2_file <- paste0("output/", model_name, ".rk45_sampling.fit.RDS")
if (run_model) {
  fit2 <- mod$sample(
    data = stan_data_rk45, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = 0, iter_sampling = iter_sampling,
    seed = 123,
    init = paste0("init/warm_init", 1:nChains, ".json"),
    step_size = step_size,
    inv_metric = fit$inv_metric(matrix = F),
    adapt_engaged = FALSE  # just in case...
  )

  fit2$cmdstan_diagnose()

  fit2$save_object(saved_fit2_file)
}


fit2 <- readRDS(saved_fit2_file)
fit2$time()

#####################################################################
## Same as above, except we use a short warmup with bdf, and do more
## warmup with rk45.

partial_warmup <- 200

saved_fit_warmup_file <- paste0("output/", model_name, ".bdf_warmup.fit.RDS")
saved_fit3_file <- paste0("output/", model_name, ".rk45_cool.fit.RDS")
if (run_model) {
  fit_warmup <- mod$sample(
    data = stan_data_bdf, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = partial_warmup, iter_sampling = 1,
    seed = 123,
    init = paste0("init/init", 1:nChains, ".json"))

  fit_warmup$save_object(saved_fit_warmup_file)
  
  samples <- fit_warmup$draws()
  parm_index <- 2:6
  for (i in 1:nChains) {
    init_saved <- list(ka = samples[1, i, parm_index[1]],
                       V = samples[1, i, parm_index[2]],
                       Vm = samples[1, i, parm_index[3]],
                       Km = samples[1, i, parm_index[4]],
                       sigma = samples[1, i, parm_index[5]])
  
    write_stan_json(init_saved, paste0("init/cool_init", i, ".json"))
  }

  mass_matrix <- fit_warmup$inv_metric()
  step_size <- fit_warmup$metadata()$step_size_adaptation

  fit3 <- mod$sample(
    data = stan_data_rk45, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup - partial_warmup,
    iter_sampling = iter_sampling,
    seed = 123,
    init = paste0("init/cool_init", 1:nChains, ".json"),
    step_size = step_size,
    inv_metric = fit_warmup$inv_metric(matrix = F)
    # inv_metric = diag(mass_matrix[[1]])
  )

  fit3$cmdstan_diagnose()  
  fit3$save_object(saved_fit3_file)
}

fit_warmup <- readRDS(saved_fit_warmup_file)
fit_warmup$time()

fit3 <- readRDS(saved_fit3_file)
fit3$time()

#####################################################################
# Make sure the posterior samples are consistent between methods.
parms <- c("ka", "V", "Vm", "Km", "sigma", "lp__")
samples0 <- as_draws_df(fit0$draws(variables = parms))
samples1 <- as_draws_df(fit$draws(variables = parms))
samples2 <- as_draws_df(fit2$draws(variables = parms))
samples3 <- as_draws_df(fit3$draws(variables = parms))

samples_rk45 <- with(samples0, c(ka, V, Vm, Km, sigma, lp__))
samples_bdf <- with(samples1, c(ka, V, Vm, Km, sigma, lp__))
samples_warm <- with(samples2, c(ka, V, Vm, Km, sigma, lp__))
samples_cool <- with(samples2, c(ka, V, Vm, Km, sigma, lp__))

samples_all <- c(samples_rk45, samples_bdf, samples_warm, samples_cool)

# parameters <- c("ka", "V", "Vm", "Km", "sigma")
total_samples <- length(samples_rk45)
parm_name <- rep(rep(parms, each = iter_sampling * nChains), 4)
method <- rep(c("rk45", "bdf", "warm", "cool"), each = total_samples)

plot_data <- data.frame(samples = samples_all,
                        parm = parm_name,
                        method = method)

plot <- ggplot(data = plot_data,
               aes(x = samples, color = method, fill = method)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_bw() + facet_wrap(~ parm, scale = "free")
plot

#####################################################################
# Examine global efficiency (include lp__) -- meaning the run time
# is taken to be the total, i.e. worst run time.
ess_bulk <- fit0$summary()$ess_bulk[1:6]
eff0 <- ess_bulk / fit0$time()$total # max(fit0$time()$chains[, 4])

ess_bulk <- fit$summary()$ess_bulk[1:6]
eff1 <- ess_bulk / fit$time()$total # max(fit$time()$chains[, 4])

ess_bulk2 <- fit2$summary()$ess_bulk[1:6]
eff2 <- ess_bulk2 / (max(fit$time()$chains[, 2]) + fit2$time()$total)

ess_bulk3 <- fit3$summary()$ess_bulk[1:6]
eff3 <- ess_bulk3 / (max(fit_warmup$time()$chains[, 2]) + fit3$time()$total)

eff <- c(eff0, eff1, eff2, eff3)
parm <- rep(parms, 4)
method <- rep(c("rk45", "bdf", "warm_start", "cool_start"), each = length(parms))
plot_data <- data.frame(eff = eff, parm = parm, method = method)

plot <- ggplot(data = plot_data,
               aes(x = parm, y = eff, fill = method)) +
  geom_bar(stat = "identity", width = 0.3, alpha = 0.8,
           position = "dodge") + 
  theme_bw() + theme(text = element_text(size = 10)) + coord_flip() +
  ylab("ESS / s") + xlab(" ")
plot


#####################################################################
# Let's now compute the mean ESS (with sd)

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

time0 <- fit0$time()$chains[, 4]
ess0 <- ess_summary(fit = fit0, parms, nChains, time0)

time1 <- fit$time()$chains[, 4]
ess1 <- ess_summary(fit = fit, parms, nChains, time1)

time2 <- fit2$time()$chains[, 3] + fit$time()$chains[, 2] 
ess2 <- ess_summary(fit = fit2, parms, nChains, time2)

time3 <- fit3$time()$chains[, 4] + fit_warmup$time()$chains[, 2]
ess3 <- ess_summary(fit = fit2, parms, nChains, time3)

# Let's examine only parameter (assuming it tells the same story)
parm_index <- 6
mean_ess <- c(ess0$mean_ess[parm_index],
              ess1$mean_ess[parm_index],
              ess2$mean_ess[parm_index], 
              ess3$mean_ess[parm_index])

sd_ess <- c(ess0$sd_ess[parm_index],
            ess1$sd_ess[parm_index],
            ess2$sd_ess[parm_index], 
            ess3$sd_ess[parm_index])

method <- c("rk45", "bdf", "warm start", "room start")

plot <- ggplot(data = data.frame(method = method, mean_ess = mean_ess),
               aes(x = method, y = mean_ess)) + theme_bw() +
  geom_bar(stat = "identity", width = 0.3, alpha = 0.8,
           position = "dodge") + theme_bw() +
  geom_errorbar(aes(ymin = mean_ess - sd_ess, ymax = mean_ess + sd_ess))
plot

# central metrics are not good summaries. Let's try plotting all the points.

recorded_eff <- c(ess0$chain_eff[parm_index, ],
                  ess1$chain_eff[parm_index, ],
                  ess2$chain_eff[parm_index, ],
                  ess3$chain_eff[parm_index, ])

method_names <- c("rk45", "bdf", "warm start", "cool start")
method <- rep(method_names, each = nChains)
method <- factor(method, levels = method_names)

plot <- ggplot(data = data.frame(method = method, eff = recorded_eff),
               aes(x = method, y = eff)) + theme_bw() +
  geom_point()
plot

# compute the relaxation time
recorded_tau <- c(ess0$chain_ess[parm_index, ] / ess0$chain_eff[parm_index, ],
                  ess1$chain_ess[parm_index, ] / ess1$chain_eff[parm_index, ],
                  ess2$chain_ess[parm_index, ] / ess2$chain_eff[parm_index, ],
                  ess3$chain_ess[parm_index, ] / ess3$chain_eff[parm_index, ])

median_tau <- c(mean(ess0$chain_ess[parm_index, ] / ess0$chain_eff[parm_index, ]),
                mean(ess1$chain_ess[parm_index, ] / ess1$chain_eff[parm_index, ]),
                mean(ess2$chain_ess[parm_index, ] / ess2$chain_eff[parm_index, ]),
                mean(ess3$chain_ess[parm_index, ] / ess3$chain_eff[parm_index, ]))

plot_data <- data.frame(method = method, tau = recorded_tau)

plot <- ggplot(data = plot_data,
               aes(x = method, y = tau)) + theme_bw() +
  geom_point() + ylab("Relaxation time (s)") + coord_flip() +
  geom_bar(stat = "identity", position = "dodge",
           aes(x = method[1], y = median_tau[1]),
           alpha = 0.01, width = 0.3) +
  geom_bar(stat = "identity", position = "dodge",
           aes(x = method[nChains + 1], y = median_tau[2]),
           alpha = 0.01, width = 0.3) +
  geom_bar(stat = "identity", position = "dodge",
           aes(x = method[2 * nChains + 1], y = median_tau[3]),
           alpha = 0.01, width = 0.3) +
  geom_bar(stat = "identity", position = "dodge",
           aes(x = method[3 * nChains + 1], y = median_tau[4]),
           alpha = 0.01, width = 0.3) +
  theme(text = element_text(size = 16))
  
plot

#####################################################################
## Additional plots (for presentation)

# Only plot rk45 and BDF
plot <- ggplot(data = plot_data[method == "rk45" | method == "bdf", ],
               aes(x = method, y = tau)) + theme_bw() +
  geom_point() + ylab("Relaxation time (s)") + coord_flip() +
  geom_bar(stat = "identity", position = "dodge",
           aes(x = method[1], y = median_tau[1]),
           alpha = 0.01, width = 0.3) +
  geom_bar(stat = "identity", position = "dodge",
           aes(x = method[nChains + 1], y = median_tau[2]),
           alpha = 0.01, width = 0.3)
plot


plot_run_time <- function(fit_object) {
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

plot_run_time(fit0) + ylim(0, 125)
plot_run_time(fit) + ylim(0, 125)

bayesplot::mcmc_trace

plot_run_time(fit2)
plot_run_time(fit3)



