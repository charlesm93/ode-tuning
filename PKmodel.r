
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
source("tools.r")


model_name <- "Michaelis_MentenPK"

#####################################################################
## Simulate data

simulate_data <- FALSE
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
## Build inits (we'll use the same inits for all the fits)

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

stan_seed <- 123

metric <- "dense_e"  # "diag_e"

init0_files <- paste0("init/init", 1:nChains, ".json")

#####################################################################
## Fit model with rk45 solver.

run_model <- TRUE  # FALSE
saved_fit0_file <- paste0("output/", model_name, "rk45")

if (run_model) {
  fit0_warmups <- fit_warmup_phases(mod = mod, 
                                    model_name = model_name,
                                    save_function = save_last_iter,
                                    data = stan_data_rk45,
                                    init = init0_files,
                                    iter_warmup = iter_warmup,
                                    init_buffer = 75,
                                    term_buffer = 50,
                                    nChains = nChains,
                                    seed = stan_seed,
                                    metric = metric)

  fit0_warmups$fit_w1$save_object(paste0(saved_fit0_file, ".phase1.fit.RDS"))
  fit0_warmups$fit_w2$save_object(paste0(saved_fit0_file, ".phase2.fit.RDS"))
  # ERROR: below should be fit_w3, not fit_w1 (FIXED 05/05)
  fit0_warmups$fit_w3$save_object(paste0(saved_fit0_file, ".phase3.fit.RDS"))
  
  samples_w1 <- fit0_warmups$fit_w1$draws()
  parm_index <- 2:6
  init_files <- save_last_iter(samples_w1, nChains, paste0(model_name, "_w3."),
                               parm_index)

  fit0 <- mod$sample(
    data = stan_data_rk45, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = 0, iter_sampling = iter_sampling,
    seed = stan_seed, adapt_delta = 0.8,
    init = init_files,
    metric = metric,
    step_size = fit0_warmups$fit_w3$metadata()$step_size_adaptation,
    inv_metric = fit0_warmups$fit_w3$inv_metric(matrix = F),
    adapt_engaged = FALSE)
  
  fit0$save_object(paste0(saved_fit0_file, ".RDS"))
}

fit0_w1 <- readRDS(paste0(saved_fit0_file, ".phase1.fit.RDS"))
fit0_w2 <- readRDS(paste0(saved_fit0_file, ".phase2.fit.RDS"))
fit0_w3 <- readRDS(paste0(saved_fit0_file, ".phase3.fit.RDS"))
fit0 <- readRDS(paste0(saved_fit0_file, ".RDS"))


#####################################################################
## Fit model with bdf solver

saved_fit1_file <- paste0("output/", model_name, "bdf")
if (run_model) {
  fit1_warmups <- fit_warmup_phases(mod = mod, 
                                    model_name = paste0(model_name, "_bdf"),
                                    data = stan_data_bdf,
                                    init = init0_files,
                                    save_function = save_last_iter,
                                    iter_warmup = iter_warmup,
                                    init_buffer = 75,
                                    term_buffer = 50,
                                    nChains = nChains,
                                    seed = stan_seed,
                                    metric = metric)
  
  fit1_warmups$fit_w1$save_object(paste0(saved_fit1_file, ".phase1.fit.RDS"))
  fit1_warmups$fit_w2$save_object(paste0(saved_fit1_file, ".phase2.fit.RDS"))
  # FIX ME -- below should fit_w3, not fit_w1 (fixed 05 / 05)
  fit1_warmups$fit_w3$save_object(paste0(saved_fit1_file, ".phase3.fit.RDS"))
  
  samples_warmup <- fit1_warmups$fit_w1$draws()
  parm_index <- 2:6
  init_files <- save_last_iter(samples_warmup, nChains, paste0(model_name, "_w3."),
                               parm_index)
  
  fit1 <- mod$sample(
    data = stan_data_bdf, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = 0, iter_sampling = iter_sampling,
    seed = stan_seed + 1, 
    adapt_delta = 0.8,
    init = init_files,
    metric = metric,
    step_size = fit1_warmups$fit_w3$metadata()$step_size_adaptation,
    inv_metric = fit1_warmups$fit_w3$inv_metric(matrix = F),
    adapt_engaged = FALSE)

  fit1$save_object(paste0(saved_fit1_file, ".RDS"))
}

fit1_w1 <- readRDS(paste0(saved_fit1_file, ".phase1.fit.RDS"))
fit1_w2 <- readRDS(paste0(saved_fit1_file, ".phase2.fit.RDS"))
fit1_w3 <- readRDS(paste0(saved_fit1_file, ".phase3.fit.RDS"))
fit1 <- readRDS(paste0(saved_fit1_file, ".RDS"))

#####################################################################
## Run sampling with rk45, using tuning parameters warmed up with bdf.

saved_fit2_file <- paste0("output/", model_name, "rk45_sampling")
if (run_model) {
  fit2 <- mod$sample(
    data = stan_data_rk45, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = 0, iter_sampling = iter_sampling,
    seed = 123, adapt_delta = 0.8,
    init = init_files,
    step_size = fit1_warmups$fit_w3$metadata()$step_size_adaptation,
    inv_metric = fit1_warmups$fit_w3$inv_metric(matrix = F),
    metric = metric,
    adapt_engaged = FALSE)

  fit2$save_object(paste0(saved_fit2_file, ".RDS"))
}

fit2 <- readRDS(paste0(saved_fit2_file, ".RDS"))
  
#####################################################################
## Run rk45 for phase 3 and sampling.
saved_fit3_file <- paste0("output/", model_name, "rk45_cool")
if (run_model) {
  init_files_w2_bdf <- paste0("init/", model_name, "_bdf_w2.", 1:nChains, ".json")

  fit_rk45_w3 <- mod$sample(data = stan_data_rk45, 
                            init = init_files_w2_bdf,
                            chains = nChains, parallel_chains = parallel_chains,
                            iter_sampling = 1,
                            iter_warmup = 50, window = 0,
                            init_buffer = 0, term_buffer = 50,
                            seed = 123,
                            step_size = fit1_warmups$fit_w2$metadata()$step_size_adaptation,
                            inv_metric = fit1_warmups$fit_w2$inv_metric(matrix = F),
                            metric = metric)
  
  fit_rk45_w3$save_object(paste0(saved_fit3_file, ".phase3.fit.RDS"))
  
  fit3 <- mod$sample(
    data = stan_data_rk45, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = 0, iter_sampling = iter_sampling,
    seed = 123, adapt_delta = 0.8,
    init = init_files,
    step_size = fit_rk45_w3$metadata()$step_size_adaptation,
    inv_metric = fit_rk45_w3$inv_metric(matrix = F),
    metric = metric,
    adapt_engaged = FALSE)

  fit3$save_object(paste0(saved_fit3_file, ".RDS"))
}

fit3_w3 <- readRDS(paste0(saved_fit3_file, ".phase3.fit.RDS"))
fit_rk45_w3 <- readRDS(paste0(saved_fit3_file, ".phase3.fit.RDS"))
fit3 <- readRDS(paste0(saved_fit3_file, ".RDS"))

#####################################################################
# Make sure the posterior samples are consistent between methods.
parms <- c("ka", "V", "Vm", "Km", "sigma", "lp__")
samples0 <- as_draws_df(fit0$draws(variables = parms))
samples1 <- as_draws_df(fit1$draws(variables = parms))
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
  geom_histogram(alpha = 0.25, position = "identity", bins = 30) +
  theme_bw() + facet_wrap(~ parm, scale = "free")
plot

#####################################################################
# Examine global efficiency (include lp__) -- meaning the run time
# is taken to be the total, i.e. worst run time.

ess_bulk <- fit0$summary()$ess_bulk[1:6]
time0_total <- (fit0$time()$total +
                fit0_w1$time()$total +
                fit0_w2$time()$total +
                fit0_w3$time()$total)
# time0_total <- (fit0$time()$total +
#                 fit0_warmups$fit_w1$time()$total +
#                 fit0_warmups$fit_w2$time()$total +
#                 fit0_warmups$fit_w3$time()$total)
eff0 <- ess_bulk / time0_total


bdf_warmup_time <- fit1_w1$time()$total +
                   fit1_w2$time()$total +
                   fit1_w3$time()$total
# bdf_warmup_time <- fit1_warmups$fit_w1$time()$total +
#                      fit1_warmups$fit_w2$time()$total +
#                      fit1_warmups$fit_w3$time()$total
ess_bulk <- fit1$summary()$ess_bulk[1:6]

time1_total <- fit1$time()$total + bdf_warmup_time
eff1 <- ess_bulk / time1_total

ess_bulk2 <- fit2$summary()$ess_bulk[1:6]
time2_total <- fit2$time()$total + bdf_warmup_time
eff2 <- ess_bulk2 / time2_total

ess_bulk3 <- fit3$summary()$ess_bulk[1:6]
time3_total <- fit3$time()$total +
  fit1_w1$time()$total +
  fit1_w2$time()$total +
  fit_rk45_w3$time()$total
# time3_total <- fit3$time()$total +
#   fit1_warmups$fit_w1$time()$total +
#   fit1_warmups$fit_w2$time()$total +
#   fit_rk45_w3$time()$total
eff3 <- ess_bulk3 / time3_total

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


###############################################################################
# Central and worse-case scenario metrics don't tell the whole story.
# Let's plot all the points (for one "representative parameter")
# Need now to extract the run time per chain.

parm_index <- 1

# time0 <- (fit0$time()$chains[, 4] +
#             fit0_warmups$fit_w1$time()$chains[, 4] +
#             fit0_warmups$fit_w2$time()$chains[, 4] +
#             fit0_warmups$fit_w3$time()$chains[, 4])
time0 <- (fit0$time()$chains[, 4] +
            fit0_w1$time()$chains[, 4] +
            fit0_w2$time()$chains[, 4] +
            fit0_w3$time()$chains[, 4])
ess0 <- ess_summary(fit = fit0, parms, nChains, time0)

bdf_warmup_time <- fit1_w1$time()$chains[, 4] +
  fit1_w2$time()$chains[, 4] +
  fit1_w3$time()$chains[, 4]
# bdf_warmup_time <- fit1_warmups$fit_w1$time()$chains[, 4] +
#   fit1_warmups$fit_w2$time()$chains[, 4] +
#   fit1_warmups$fit_w3$time()$chains[, 4]

time1 <- fit1$time()$chains[, 4] + bdf_warmup_time
ess1 <- ess_summary(fit = fit1, parms, nChains, time1)

time2 <- fit2$time()$chains[, 4] + bdf_warmup_time
ess2 <- ess_summary(fit = fit2, parms, nChains, time2)

# time3 <- fit3$time()$chains[, 4] + fit1_warmups$fit_w1$time()$chains[, 4] +
#   fit1_warmups$fit_w2$time()$chains[, 4] +
#   fit_rk45_w3$time()$chains[, 4]
time3 <- fit3$time()$chains[, 4] + fit1_w1$time()$chains[, 4] +
  fit1_w2$time()$chains[, 4] +
  fit_rk45_w3$time()$chains[, 4]
ess3 <- ess_summary(fit = fit3, parms, nChains, time3)

recorded_eff <- c(ess0$chain_eff[parm_index, ],
                  ess1$chain_eff[parm_index, ],
                  ess2$chain_eff[parm_index, ],
                  ess3$chain_eff[parm_index, ])

method_names <- c("RK45", "BDF", "Late switch", "Early switch")
method <- rep(method_names, each = nChains)
method <- factor(method, levels = method_names)

plot <- ggplot(data = data.frame(method = method, eff = recorded_eff),
               aes(x = method, y = eff)) + theme_bw() +
  geom_point()
plot

# compute the relaxation time
# recorded_tau <- c(ess0$chain_ess[parm_index, ] / ess0$chain_eff[parm_index, ],
#                   ess1$chain_ess[parm_index, ] / ess1$chain_eff[parm_index, ],
#                   ess2$chain_ess[parm_index, ] / ess2$chain_eff[parm_index, ],
#                   ess3$chain_ess[parm_index, ] / ess3$chain_eff[parm_index, ])
# 
# median_tau <- c(median(ess0$chain_ess[parm_index, ] / ess0$chain_eff[parm_index, ]),
#                 median(ess1$chain_ess[parm_index, ] / ess1$chain_eff[parm_index, ]),
#                 median(ess2$chain_ess[parm_index, ] / ess2$chain_eff[parm_index, ]),
#                 median(ess3$chain_ess[parm_index, ] / ess3$chain_eff[parm_index, ]))

recorded_tau <- c(1 / ess0$chain_eff[parm_index, ],
                  1 / ess1$chain_eff[parm_index, ],
                  1 / ess2$chain_eff[parm_index, ],
                  1 / ess3$chain_eff[parm_index, ])

median_tau <- c(median(1 / ess0$chain_eff[parm_index, ]),
                median(1 / ess1$chain_eff[parm_index, ]),
                median(1 / ess2$chain_eff[parm_index, ]),
                median(1 / ess3$chain_eff[parm_index, ]))

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

# To capture outliers, want to use plot on the log scale.
plot <- ggplot(data = plot_data,
               aes(x = method, y = tau)) + theme_bw() +
  geom_point() + scale_y_continuous(trans = 'log10') +
  geom_point(aes(x = method[1], y = median_tau[1]), shape = 13,
             size = 3, color = "orange") +
  geom_point(aes(x = method[nChains + 1], y = median_tau[2]), shape = 13,
             size = 3, color = "orange") +
  geom_point(aes(x = method[2 * nChains + 1], y = median_tau[3]), shape = 13,
             size = 3, color = "orange") +
  geom_point(aes(x = method[3 * nChains + 1], y = median_tau[4]), shape = 13,
             size = 3, color = "orange") +
  ylab("Relaxation time (s)") + xlab("") + coord_flip() +
  theme(text = element_text(size = 16))
plot

# Write data for plot
write.csv(plot_data, paste0("plot_data/", model_name, ".data.csv"))
write.csv(median_tau, paste0("plot_data/", model_name,"_tau_median",
                             ".data.csv"))

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


# time0_data <- c(fit0_warmups$fit_w1$time()$chains[, 4],
#                 fit0_warmups$fit_w2$time()$chains[, 4],
#                 fit0_warmups$fit_w3$time()$chains[, 4],
#                 fit0$time()$chains[, 4])
time0_data <- c(fit0_w1$time()$chains[, 4],
                fit0_w2$time()$chains[, 4],
                fit0_w3$time()$chains[, 4],
                fit0$time()$chains[, 4])

# time1_data <- c(fit1_warmups$fit_w1$time()$chains[, 4],
#                 fit1_warmups$fit_w2$time()$chains[, 4],
#                 fit1_warmups$fit_w3$time()$chains[, 4],
#                 fit1$time()$chains[, 4])
time1_data <- c(fit1_w1$time()$chains[, 4],
                fit1_w2$time()$chains[, 4],
                fit1_w3$time()$chains[, 4],
                fit1$time()$chains[, 4])

chain_id <- rep(1:nChains, 4)
phase <- rep(c("phase 1", "phase 2", "phase 3", "sampling"), each = nChains)
phase <- factor(phase, levels = c("sampling", "phase 3", "phase 2", "phase 1"))
# phase <- factor(phase, levels = c("phase 1", "phase 2", "phase 3", "sampling"))

plot_phase_time(run_times = time0_data, chain_id, phase) 
plot_phase_time(run_times = time1_data, chain_id, phase)

# Put plots together for PAGE poster
time_all_data <- c(time0_data, time1_data)
chain_id_all <- c(chain_id, chain_id)
phase_all <- c(phase, phase)
integrator <- c(rep("RK45", length(time0_data)),
                rep("BDF", length(time0_data)))

plot_phase_time(time_all_data, chain_id, phase, integrator)

# TODO: change color scheme of ggplot2
if (FALSE) {
  plot_run_time(fit0) + ylim(0, 125)
  plot_run_time(fit) + ylim(0, 125)

  plot_run_time(fit2)
  plot_run_time(fit3)
}
