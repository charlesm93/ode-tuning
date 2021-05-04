
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
  ka_pop <- 3
  V_pop <- 27
  Vm_pop <- 10
  km_pop <- 14
  sigma <- 0.1
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

  n_patients <- 10
  times <- rep(c(seq(0.5, 1.5, by = 0.25), seq(10, 100, by = 10)),
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
    geom_line() + geom_point() + theme_bw() + facet_wrap(~patient_ID)
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
init <- function() {
  list(
    ka_pop = exp(rnorm(1, log(1), 0.5)),
    V_pop = exp(rnorm(1, log(35), 0.5)),
    Vm_pop = exp(rnorm(1, log(10), 0.5)),
    Km_pop = exp(rnorm(1, log(2.5), 0.25)),
    sigma = abs(rnorm(1, 0, 1)),
    omega = exp(rnorm(5, log(0.25), 0.1)),
    eta = matrix(rnorm(stan_data_rk45$n_patients * 4,
                       0, 1), 4, stan_data_rk45$n_patients),
    # theta = matrix(rep(c(1, 35, 10, 2.5), stan_data_rk45$n_patients),
    #                4, stan_data_rk45$n_patients)
    theta = exp(matrix(rnorm(stan_data_rk45$n_patients * 4,
                         c(log(1), log(35), log(10), log(2.5)),
                         0.1), 4,
                         stan_data_rk45$n_patients))
  )
}

nChains <- 4
parallel_chains <- 4

if (FALSE) {
  for (i in 1:nChains) {
    init0 <- init()
    write_stan_json(init0, paste0("init/init_pop", i, ".json"))
  }
}

#####################################################################
# Tuning parameters for MCMC + compile model

iter_warmup <- 500
iter_sampling <- 500

## Compile model
if (TRUE) {
  file <- paste0("model/", model_name, ".stan")
  mod <- cmdstan_model(file)
}

stan_seed <- 123
metric <- "diag_e"  # dense_e  CHECK -- which metric should we use?
init0_files <- paste0("init/init_pop", 1:nChains, ".json")

#####################################################################
## Fitmodel with bdf solver

run_model <- TRUE
# saved_fit0_file <- paste0("output/", model_name, "bdf")
if (FALSE) {
  fit_attempt <- mod$sample(iter_warmup = 500, iter_sampling = 500,
                            data = stan_data_rk45, init = init0_files,
                            chains = 4, parallel_chains = 4,
                            seed = stan_seed, max_treedepth = 11)

  fit_attempt$save_object("output/fit_attempt_pop_rk45.RDS")

  fit_attempt$time()
}

if (FALSE) {
  fit_attempt <- mod$sample(iter_warmup = 500, iter_sampling = 500,
                            data = stan_data_bdf, init = init0_files,
                            chains = 4, parallel_chains = 4,
                            seed = stan_seed, max_treedepth = 11)

  # Centered
  if (FALSE) fit_attempt$save_object("output/fit_attempt_pop_c_bdf.RDS")

  # Non-centered
  if (FALSE) fit_attempt$save_object("output/fit_attempt_pop_bdf.RDS")

  fit_attempt$time()
}
  

fit_rk45 <- readRDS("output/fit_attempt_pop_rk45.RDS")
fit_bdf <- readRDS("output/fit_attempt_pop_bdf.RDS")
fit_bdf_c <- readRDS("output/fit_attempt_pop_c_bdf.RDS")

fit_rk45$summary()
fit_bdf$summary()
fit_bdf_c$summary()

fit_rk45$time()
fit_bdf$time()

# Let's examine the wamred up tuning parameters
step_size_rk45 <- fit_rk45$metadata()$step_size_adaptation
inv_metric_raw <- fit_rk45$inv_metric()
inv_metric_rk45 <- list()
inv_metric_rk45$`1` <- diag(inv_metric_raw$`1`)
inv_metric_rk45$`2` <- diag(inv_metric_raw$`2`)
inv_metric_rk45$`3` <- diag(inv_metric_raw$`3`)
inv_metric_rk45$`4` <- diag(inv_metric_raw$`4`)

step_size_bdf <- fit_bdf$metadata()$step_size_adaptation
inv_metric_raw <- fit_bdf$inv_metric()
inv_metric_bdf <- list()
inv_metric_bdf$`1` <- diag(inv_metric_raw$`1`)
inv_metric_bdf$`2` <- diag(inv_metric_raw$`2`)
inv_metric_bdf$`3` <- diag(inv_metric_raw$`3`)
inv_metric_bdf$`4` <- diag(inv_metric_raw$`4`)

# Step sizes are similar. Difficult to compare inverse metric...
# Let's try running the rk45 with the bdf warmup

fit_warm <- mod$sample(iter_warmup = 0, iter_sampling = 50,
                       data = stan_data_rk45, init = init0_files,
                       adapt_engaged = FALSE,
                       chains = 4, parallel_chains = 4,
                       seed = stan_seed, max_treedepth = 11,
                       step_size = step_size_bdf,
                       inv_metric = inv_metric_bdf)

