
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

model_name <- "Michaelis_MentenPK_pop_centered"

#####################################################################
## Simulate data

simulate_data <- FALSE
if (simulate_data) {
  ka_pop <- 3.2  # 3
  V_pop <- 27
  Vm_pop <- 10
  km_pop <- 14
  sigma <- 0.1 # 0.1
  y0 <- 100
  
  parm_pop <- c(ka_pop, V_pop, Vm_pop, km_pop)
  state <- c(Y1 = y0, Y2 = 0)
  n_parm <- length(parm_pop)

  dydt <- function(t, state, parm) {
    with(as.list(c(state)), {
      ka = parm[1]
      V = parm[2]
      Vm = parm[3]
      Km = parm[4]

      C = Y2 / V
      dY1 <- - ka * Y1
      dY2 <- ka * Y1 - Vm * C / (Km + C)
      list(c(dY1, dY2))
    })
  }

  n_patients <- 3
  times <- rep(c(seq(0.25, 1.5, by = 0.25), seq(10, 100, by = 10)),
               n_patients)
  # times <- rep(c(seq(0.25, 1.5, by = 0.25), seq(10, 30, by = 10)),
  #              n_patients)
  N <- length(times)
  n_obs <- N / n_patients

  start <- seq(from = 1, to = N, by = n_obs)
  end <- seq(from = n_obs, to = N, by = n_obs)

  omega <- exp(rnorm(1, 0.25, 1))
  theta <- rep(NA, n_parm)
  concentration <- rep(NA, N)
  for (i in 1:n_patients) {
    eta <- rnorm(n_parm, 0, 1)
    theta <- parm_pop + omega * eta
    
    mass <- ode(y = state, time = c(0, times[start[i]:end[i]]),
                func = dydt, parms = theta)
    concentration[start[i]:end[i]] <- mass[-1, 3] / theta[2]
  }
  
  y <- exp(rnorm(N, log(concentration), sigma))
  
  patient_ID <- rep(1:n_patients, each = n_obs)
  
  p <- ggplot(data = data.frame(y = y, times = times, patient_ID = patient_ID),
              aes(y = y, x = times)) +
    geom_point() + theme_bw() + facet_wrap(~patient_ID)
  p
  
  stan_data <- list(N = N, n_obs = n_obs, n_patients = n_patients,
                    y = y, t = times, y0 = c(y0, 0),
                    start = start, end = end)
  write_stan_json(stan_data, "data/PKModel_pop.data.json")
}

stan_data_rk45 <- fromJSON(file = "data/PKModel_pop.data.json")
stan_data_rk45$stiff_solver <- 0
stan_data_bdf <- stan_data_rk45
stan_data_bdf$stiff_solver <- 1

# Build inits (we'll use the same inits for all the fits)
if (FALSE) {
init <- function() {
  list(
    ka_pop = exp(rnorm(1, log(2.5), 3)),
    V_pop = exp(rnorm(1, log(35), 0.5)),
    Vm_pop = exp(rnorm(1, log(10), 0.5)),
    Km_pop = exp(rnorm(1, log(2.5), 3)),
    sigma = abs(rnorm(1, 0, 1)),
    omega = exp(rnorm(5, log(0.25), 0.1)),
    eta = matrix(rnorm(stan_data_rk45$n_patients * 4,
                       0, 1), 4, stan_data_rk45$n_patients),
    theta = exp(matrix(rnorm(stan_data_rk45$n_patients * 4,
                         c(log(1), log(35), log(10), log(2.5)),
                         0.1), 4,
                         stan_data_rk45$n_patients))
  )
}
}

nChains <- 8
parallel_chains <- 8

create_init <- FALSE
if (create_init) {
  for (i in 1:nChains) {
    init0 <- init()
    write_stan_json(init0, paste0("init/init_pop", i, ".json"))
  }
}

#####################################################################
# Tuning parameters for MCMC + compile model

iter_warmup <- 1000
iter_sampling <- 1000

## Compile model
if (TRUE) {
  file <- paste0("model/", model_name, ".stan")
  mod <- cmdstan_model(file)
}

stan_seed <- 123
metric <- "dense_e"  # "diag_e" # CHECK -- which metric should we use?
init0_files <- paste0("init/init_pop", 1:nChains, ".json")

#####################################################################
## Fitmodel with bdf solver

run_model <- FALSE
saved_fit0_file <- paste0("output/", model_name, "rk45")
if (run_model) {
  fit0_warmups <- fit_warmup_phases(mod = mod,
                                    model_name = model_name,
                                    data = stan_data_rk45,
                                    save_function = save_last_iter_pop,
                                    parm_index = 2:50,
                                    init = init0_files,
                                    iter_warmup = iter_warmup,  # 500,
                                    init_buffer = 75,           # 75,
                                    term_buffer = 50 ,          #50,
                                    nChains = nChains,
                                    seed = stan_seed,
                                    metric = metric,
                                    adapt_delta = 0.95)

  fit0_warmups$fit_w1$save_object(paste0(saved_fit0_file, ".phase1.fit.RDS"))
  fit0_warmups$fit_w2$save_object(paste0(saved_fit0_file, ".phase2.fit.RDS"))
  fit0_warmups$fit_w3$save_object(paste0(saved_fit0_file, ".phase3.fit.RDS"))

  samples_w1 <- fit0_warmups$fit_w1$draws()
  parm_index <- 2:50
  init_files <- save_last_iter_pop(samples_w1, nChains, paste0(model_name, "_w3."),
                                   parm_index)
  
  fit0 <- mod$sample(
    data = stan_data_rk45, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = 0, iter_sampling = iter_sampling,
    seed = stan_seed + 1, adapt_delta = 0.95,
    init = init_files,
    metric = metric,
    step_size = fit0_warmups$fit_w3$metadata()$step_size_adaptation,
    inv_metric = fit0_warmups$fit_w3$inv_metric(matrix = F),
    adapt_engaged = FALSE, refresh = 10)

  fit0$save_object(paste0(saved_fit0_file, ".RDS"))
}

fit0_w1 <- readRDS(paste0(saved_fit0_file, ".phase1.fit.RDS"))
fit0_w2 <- readRDS(paste0(saved_fit0_file, ".phase2.fit.RDS"))
fit0_w3 <- readRDS(paste0(saved_fit0_file, ".phase3.fit.RDS"))
fit0 <- readRDS(paste0(saved_fit0_file, ".RDS"))


#####################################################################
## Fit model with bdf solver

saved_fit1_file <- paste0("output/", model_name, "_bdf")
if (run_model) {
  fit1_warmups <- fit_warmup_phases(mod = mod, 
                                    model_name = paste0(model_name, "_bdf"),
                                    data = stan_data_bdf,
                                    init = init0_files,
                                    save_function = save_last_iter_pop,
                                    parm_index = 2:50,
                                    iter_warmup = iter_warmup,
                                    init_buffer = 75,
                                    term_buffer = 50,
                                    nChains = nChains,
                                    seed = stan_seed,
                                    metric = metric,
                                    adapt_delta = 0.99)
  
  fit1_warmups$fit_w1$save_object(paste0(saved_fit1_file, ".phase1.fit.RDS"))
  fit1_warmups$fit_w2$save_object(paste0(saved_fit1_file, ".phase2.fit.RDS"))
  fit1_warmups$fit_w3$save_object(paste0(saved_fit1_file, ".phase3.fit.RDS"))

  samples_warmup <- fit1_warmups$fit_w1$draws()
  parm_index <- 2:50
  init_files <- save_last_iter_pop(samples_warmup, nChains, paste0(model_name, "_w3."),
                               parm_index)

  fit1 <- mod$sample(
    data = stan_data_bdf, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = 0, iter_sampling = iter_sampling,
    seed = stan_seed + 1, 
    adapt_delta = 0.9,
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
    seed = 123, adapt_delta = 0.9,
    init = init_files,
    step_size = fit1_warmups$fit_w3$metadata()$step_size_adaptation,
    inv_metric = fit1_warmups$fit_w3$inv_metric(matrix = F),
    metric = metric,
    adapt_engaged = FALSE,
    refresh = 10)

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
                            adapt_delta = 0.95,
                            step_size = fit1_warmups$fit_w2$metadata()$step_size_adaptation,
                            inv_metric = fit1_warmups$fit_w2$inv_metric(matrix = F),
                            metric = metric)
  
  fit_rk45_w3$save_object(paste0(saved_fit3_file, ".phase3.fit.RDS"))
  
  fit3 <- mod$sample(
    data = stan_data_rk45, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = 0, iter_sampling = iter_sampling,
    seed = 123, adapt_delta = 0.9,
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
parms <- c("ka_pop", "V_pop", "Vm_pop", "Km_pop", "sigma", "lp__")
samples0 <- as_draws_df(fit0$draws(variables = parms))
samples1 <- as_draws_df(fit1$draws(variables = parms))
samples2 <- as_draws_df(fit2$draws(variables = parms))
samples3 <- as_draws_df(fit3$draws(variables = parms))

samples_rk45 <- with(samples0, c(ka_pop, V_pop, Vm_pop, Km_pop, sigma, lp__))
samples_bdf <- with(samples1, c(ka_pop, V_pop, Vm_pop, Km_pop, sigma, lp__))
samples_warm <- with(samples2, c(ka_pop, V_pop, Vm_pop, Km_pop, sigma, lp__))
samples_cool <- with(samples2, c(ka_pop, V_pop, Vm_pop, Km_pop, sigma, lp__))

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

###############################################################################
# Central and worse-case scenario metrics don't tell the whole story.
# Let's plot all the points (for one "representative parameter", i.e. lp__)
# Need now to extract the run time per chain.

parm_index <- 1  # lp__
# time0 <- (fit0$time()$chains[, 4] +
#             fit0_warmups$fit_w1$time()$chains[, 4] +
#             fit0_warmups$fit_w2$time()$chains[, 4] +
#             fit0_warmups$fit_w3$time()$chains[, 4])
time0 <- (fit0$time()$chains[, 4] +
            fit0_w1$time()$chains[, 4] +
            fit0_w2$time()$chains[, 4] +
            fit0_w3$time()$chains[, 4])
ess0 <- ess_summary(fit = fit0, parms, nChains, time0)

# bdf_warmup_time <- fit1_warmups$fit_w1$time()$chains[, 4] +
#   fit1_warmups$fit_w2$time()$chains[, 4] +
#   fit1_warmups$fit_w3$time()$chains[, 4]
bdf_warmup_time <- fit1_w1$time()$chains[, 4] +
  fit1_w2$time()$chains[, 4] +
  fit1_w3$time()$chains[, 4]

time1 <- fit1$time()$chains[, 4] + bdf_warmup_time
ess1 <- ess_summary(fit = fit1, parms, nChains, time1)


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
recorded_tau <- c(1 / ess0$chain_eff[parm_index, ],
                  1 / ess1$chain_eff[parm_index, ],
                  1 / ess2$chain_eff[parm_index, ],
                  1 / ess3$chain_eff[parm_index, ])

# median_tau <- c(median(ess0$chain_ess[parm_index, ] / ess0$chain_eff[parm_index, ]),
#                 median(ess1$chain_ess[parm_index, ] / ess1$chain_eff[parm_index, ]),
#                 median(ess2$chain_ess[parm_index, ] / ess2$chain_eff[parm_index, ]),
#                 median(ess3$chain_ess[parm_index, ] / ess3$chain_eff[parm_index, ]))
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

write.csv(plot_data, paste0("plot_data/", model_name, ".data.csv"))
write.csv(median_tau, paste0("plot_data/", model_name,"_tau_median",
                             ".data.csv"))

#####################################################################
## Plot run time

# Only plot rk45 and BDF
time0_data <- c(fit0_warmups$fit_w1$time()$chains[, 4],
                fit0_warmups$fit_w2$time()$chains[, 4],
                fit0_warmups$fit_w3$time()$chains[, 4],
                fit0$time()$chains[, 4])

time1_data <- c(fit1_warmups$fit_w1$time()$chains[, 4],
                fit1_warmups$fit_w2$time()$chains[, 4],
                fit1_warmups$fit_w3$time()$chains[, 4],
                fit1$time()$chains[, 4])

chain_id <- rep(1:nChains, 4)
phase <- rep(c("Phase 1", "Phase 2", "Phase 3", "sampling"), each = nChains)
phase <- factor(phase, levels = c("sampling", "Phase 3", "Phase 2", "Phase 1"))

plot_phase_time(run_times = time0_data, chain_id, phase)  # rk45
plot_phase_time(run_times = time1_data, chain_id, phase)  # bdf

