
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
# source("stanTools.r")

model_name <- "Michaelis_MentenPK"

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
  ka <- 2.8
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
  # plot(x = mass[, 1], y = y)
}

stan_data <- fromJSON(file = "data/PKModel.data.json")

init <- function() {
  list(
    ka = exp(rnorm(1, log(1.5), 3)),
    V = exp(rnorm(1, log(35), 0.5)),
    Vm = exp(rnorm(1, log(10), 0.5)),
    Km = exp(rnorm(1, log(2.5), 3)),
    sigma = abs(rnorm(1, 0, 1))
  )
}

# init <- function() { 
#   list(
#     ka = exp(rnorm(1, log(2.5), 1)),
#     V = exp(rnorm(1, log(35), 0.5)), 
#     Vm = exp(rnorm(1, log(10), 0.5)), 
#     Km = exp(rnorm(1, log(2.5), 1.5)), 
#     sigma = abs(rnorm(1, 0, 0.1))) 
# }

nChains <- 8
iter_warmup <- 1000
iter_sampling <- 1000

## Compile model
if (TRUE) {
  file <- paste0("model/", model_name, ".stan")
  mod <- cmdstan_model(file)
}


stan_data$stiff_solver <- 0
fit0 <- mod$sample(
    data = stan_data, chains = nChains,
    parallel_chains = nChains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    seed = 123, adapt_delta = 0.8,
    init = init)

fit0$save_object(paste0("output/", model_name, ".rk45.fit.RDS"))

fit0$cmdstan_diagnose()
fit0$time()

stan_data$stiff_solver <- 1
fit <- mod$sample(
  data = stan_data, chains = nChains,
  parallel_chains = nChains,
  iter_warmup = iter_warmup, iter_sampling = iter_sampling,
  seed = 123, adapt_delta = 0.8, init = init)

fit$save_object(paste0("output/", model_name, ".bdf.fit.RDS"))

fit$cmdstan_diagnose()
fit$time()

mass_matrix <- fit$inv_metric()
step_size <- fit$metadata()$step_size_adaptation

stan_data2 <- stan_data
stan_data2$stiff_solver <- 0

# Create init files, using warmup samples and save them.
# CHECK -- the draws should be ordered
samples <- fit$draws()
parm_index <- 2:6
for (i in 1:nChains) {
  init <- list(ka = samples[1, i, parm_index[1]],
               V = samples[1, i, parm_index[2]],
               Vm = samples[1, i, parm_index[3]],
               Km = samples[1, i, parm_index[4]],
               sigma = samples[1, i, parm_index[5]])

  write_stan_json(init, paste0("init/init", i, ".json"))
}

fit2 <- mod$sample(
  data = stan_data2, chains = nChains,
  parallel_chains = nChains,
  iter_warmup = 0, iter_sampling = iter_sampling,
  seed = 123,
  init = paste0("init/init", 1:nChains, ".json"),
  step_size = step_size,
  inv_metric = diag(mass_matrix[[1]]),
  adapt_engaged = FALSE  # just in case...
)

fit2$save_object(paste0("output/", model_name,
                        ".rk45_sampling.fit.RDS"))

fit2$cmdstan_diagnose()
fit2$time()

#####################################################################
# Make the posterior samples are consistent between methods.
parms <- c("ka", "V", "Vm", "Km", "sigma")
samples <- as_draws_df(fit$draws(variables = parms))
samples2 <- as_draws_df(fit2$draws(variables = parms))

samples_bdf <- with(samples, c(ka, V, Vm, Km, sigma))
samples_mixed <- with(samples2, c(ka, V, Vm, Km, sigma))
samples_all <- c(samples_bdf, samples_mixed)

# parameters <- c("ka", "V", "Vm", "Km", "sigma")
total_samples <- length(samples_bdf)
parm_name <- rep(rep(parms, each = iter_warmup * nChains), 2)
method <- rep(c("bdf", "mixed"), each = total_samples)

plot_data <- data.frame(samples = samples_all,
                        parm = parm_name,
                        method = method)

plot <- ggplot(data = plot_data,
               aes(x = samples, color = method, fill = method)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_bw() + facet_wrap(~ parm, scale = "free")
plot

#####################################################################
# Examine efficiency (we care about parm 2 - 6)
ess_bulk <- fit0$summary()$ess_bulk[2:6]
eff0 <- ess_bulk / fit0$time()$total # max(fit0$time()$chains[, 4])

ess_bulk <- fit$summary()$ess_bulk[2:6]
eff1 <- ess_bulk / fit$time()$total # max(fit$time()$chains[, 4])

ess_bulk2 <- fit2$summary()$ess_bulk[2:6]
eff2 <- ess_bulk2 / max(fit$time()$chains[, 2]) + fit2$time()$total

eff <- c(eff0, eff1, eff2)
parm <- rep(parms, 3)
method <- rep(c("rk45", "bdf", "mixed"), each = length(parms))
plot_data <- data.frame(eff = eff, parm = parm, method = method)

plot <- ggplot(data = plot_data,
               aes(x = parm, y = eff, fill = method)) +
  geom_bar(stat = "identity", width = 0.3, alpha = 0.8,
           position = "dodge") + 
  theme_bw() + theme(text = element_text(size = 10)) + coord_flip() +
  ylab("ESS / s") + xlab(" ")
plot
