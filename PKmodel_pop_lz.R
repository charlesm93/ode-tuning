
rm(list = ls())
gc()
set.seed(1954)

library(ggplot2)
library(cmdstanr)
library(rjson)
library(posterior)
library(deSolve)
library(bayesplot)

model_name <- "Michaelis_MentenPK_pop"

#####################################################################
## load data
stan_data_rk45 <- fromJSON(file = "data/PKModel_pop.data.json")
stan_data_rk45$stiff_solver <- 0
stan_data_bdf <- stan_data_rk45
stan_data_bdf$stiff_solver <- 1

nChains <- 4
parallel_chains <- 4

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
## Fit model with bdf solver

run_model <- FALSE
saved_fit0_file <- paste0("output/", model_name, "bdf")
if (run_model) {
  fit_attempt <- mod$sample(iter_warmup = 500, iter_sampling = 500,
                            data = stan_data_bdf, init = init0_files,
                            chains = 4, parallel_chains = 4,
                            seed = stan_seed, max_treedepth = 11)
  
  fit_attempt$save_object("output/fit_attempt_pop_bdf.RDS")
  fit_attempt$time()
}


fit_rk45 <- readRDS("output/fit_attempt_pop_rk45.RDS")
fit_bdf <- readRDS("output/fit_attempt_pop_bdf.RDS")

fit_rk45$summary()
fit_bdf$summary()

fit_rk45$time()
fit_bdf$time()

mcmc_trace(fit_rk45$draws("lp__"), iter1 = 1) 
mcmc_trace(fit_bdf$draws("lp__"), iter1 = 1) 

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


## #######################################################################
### initialize with pathfinder ###
library(rstan)
library(parallel)
library(foreach)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("./utils/sim_pf.R")
N1 = 1000    # maximum iters in optimization
factr_tol = 1e2 # relative tolerance = 1-4 is not enough, should use at least 1e7
N_sam_DIV = 3   # samples for ELBO evaluation
N_sam = 100
lmm = 6 # histogram size
mc.cores = parallel::detectCores() - 2
MC = 10    
init_bound = 2
seed_list = 1:MC

model <- stan_model(file)

fit_pathfinder = FALSE
run_model_pf = FALSE

#######################################################################
## fit pathfinder with bdf ##
if(fit_pathfinder){
  
  t <- proc.time()
  opath_bdf <- opt_path_stan_parallel(seed_list, seed_list, mc.cores, model,
                                  stan_data_bdf, N1, N_sam_DIV, N_sam,
                                  init_bound = 2.0, factr_tol, lmm)
  pf_bdf_time <- (proc.time() - t)
  
  pick_samples_bdf <- Imp_Resam_WOR(opath_bdf, n_inits = 4, seed = 1)
  
  posterior <- to_posterior(model, stan_data_bdf)
  init_pf = apply(pick_samples_bdf, 2, 
                  f <- function(x){constrain_pars(posterior, x)})
  
  for (i in 1:nChains) {
    init0 <- init_pf[[i]]
    write_stan_json(init0, 
                    paste0("./init/initpfbdf_pop", i, ".json"))
  }
  
  save(list = c("pf_bdf_time"),
       file = paste0("./output/pf_bdf_time_pop.RData"))
  
}


#####################################################################
## Fit model with bdf solver and initials from Pathfinder

saved_fit_pf_bdf_file <- 
  paste0("./output/", model_name, ".pf_bdf.fit.RDS")
if (run_model_pf) {
  fit_pf_bdf <- mod$sample(
    data = stan_data_bdf, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    seed = stan_seed, max_treedepth = 11,
    init = paste0("./init/initpfbdf_pop", 1:nChains, ".json"))
  
  fit_pf_bdf$cmdstan_diagnose()
  fit_pf_bdf$save_object(saved_fit_pf_bdf_file)
  
  samples <- fit_pf_bdf$draws()
  mass_matrix <- fit_pf_bdf$inv_metric(matrix = F)
  step_size <- fit_pf_bdf$metadata()$step_size_adaptation
  
}

fit_pf_bdf <- readRDS(saved_fit_pf_bdf_file)
fit_pf_bdf$time()
fit_bdf$time()

table(fit_pf_bdf$sampler_diagnostics()[, , "treedepth__"])
table(fit_bdf$sampler_diagnostics()[, , "treedepth__"])


