
rm(list = ls())
gc()
set.seed(1954)

setwd("~/Desktop/Code/ode-experiment/")
.libPaths("~/Rlib/")

library(cmdstanr)
set_cmdstan_path("~/Desktop/classes/Teaching/Stan-and-Torsten_2/Torsten/cmdstan/")
library(rstan)
library(parallel)
source("stanTools.r")

data <- read_rdump("data/FKModel.data.r")
data$stiff_solver <- 1

# Specify initial conditions: check if needed.
init <- function() {
  list(CL = rlnorm(1, log(10), 0.2),
       Q = rlnorm(1, log(15), 0.2),
       VC = rlnorm(1, log(35), 0.2),
       VP = rlnorm(1, log(105), 0.2),
       ka = rlnorm(1, log(2), 0.2),
       sigma = abs(rcauchy(1, 0, 1)),
       #
       mtt = rlnorm(1, log(125), 0.2),
       circ0 = rlnorm(1, log(5), 0.05),
       gamma = rlnorm(1, log(0.17), 0.1),
       sigmaNeut = abs(rcauchy(1, 0, 1)))
}

# compile model using cmdstanr
if (TRUE) {
  file <- file.path("model/FKModel.stan")
  mod <- cmdstan_model(file)
}
 
nChains <- 3
num_cores <- min(nChains, detectCores())

fit <- mod$sample(
  data = data, num_chains = nChains, num_cores = num_cores,
  num_warmup = 1000, num_samples = 500, seed = 234,
  adapt_delta = 0.8)

stanfit <- read_stan_csv(fit$output_files())
saveRDS(stanfit, file = file.path("output", "fit.RSave"))

# extract warmed up parameters
sampler_parameters <- get_sampler_params(stanfit, inc_warmup = FALSE)
divergence_by_chain <- sapply(sampler_parameters, 
                              function(x) sum(x[, "divergent__"]))
step_size_by_chain <- sapply(sampler_parameters,
                             function(x) mean(x[, "stepsize__"]))
mass_matrix <- get_mass_matrix(stanfit, nChains = 1)

get_elapsed_time(stanfit)
# warmup  sample
# chain:1 680.028 743.755

data2 <- data
data2$stiff_solver <- 0
samples <- extract(stanfit)
init <- list(Cl = samples$CL[1],
       Q = samples$Q[1],
       VC = samples$VC[1],
       VP = samples$VP[1],
       ka = samples$ka[1],
       sigma = samples$sigma[1],
       mtt = samples$mtt[1],
       circ0 = samples$circ0[1],
       sigmaNeut = samples$sigmaNeut[1])

with(init, stan_rdump(ls(init), "output/init.R"))

fit2 <- mod$sample(
  data = data2, num_chains = 2, num_cores = 2,
  num_warmup = 0, num_samples = 500, seed = 123,
  init = file.path("output", "init.R"),
  stepsize = step_size_by_chain,
  inv_metric = mass_matrix[[1]],
  adapt_engaged = FALSE  # just in case...
)

stanfit2 <- read_stan_csv(fit2$output_files())
get_elapsed_time(stanfit2)

# Let's compare the samples
samples2 <- extract(stanfit2)

samples_bdf <- with(samples, c(CL, Q, VC, VP, ka, log(sigma), mtt, circ0, 
                               log(sigmaNeut)))
samples_rk45 <- with(samples2, c(CL, Q, VC, VP, ka, log(sigma), mtt, circ0,
                                 log(sigmaNeut)))
samples_total <- c(samples_bdf, samples_rk45)
parameters <- c("CL", "Q", "VC", "VP", "ka", "log_sigma", "mtt", "circ0",
                "log_sigmaNeut")
parm <- rep(rep(parameters, each = 500), 2)
method <- rep(c("bdf", "rk45"), each = 500 * length(parameters))

plot <- ggplot(data = data.frame(samples = samples_total, parm = parm,
                                 method = method),
               aes(x = samples, color = method, fill = method)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_bw() + facet_wrap(~ parm, scale = "free")
plot

parameters[6] <- "sigma"
parameters[9] <- "sigmaNeut"
summary(stanfit, parameters)[[1]]
summary(stanfit2, parameters)[[1]]

