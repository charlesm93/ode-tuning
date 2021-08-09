
rm(list = ls())
gc()
set.seed(1954)

#setwd("~/Code/ode-tuning")
#.libPaths("~/Rlib/")

library(ggplot2)
library(cmdstanr)
#set_cmdstan_path("~/Rlib/cmdstan/")
library(rjson)
library(posterior)
library(deSolve)
library(bayesplot)
source("tools.r")

model_name <- "Michaelis_MentenPK_pop_centered"

#####################################################################
## load Simulate data
stan_data_rk45 <- fromJSON(file = "data/PKModel_pop.data.json")
stan_data_rk45$stiff_solver <- 0
stan_data_bdf <- stan_data_rk45
stan_data_bdf$stiff_solver <- 1

nChains <- 8
parallel_chains <- 8

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

#######################################################################
### warmup with pathfinder ###
source("./utils/sim_pf.R")
init_bound = 2.0 # bandwidth of initial distribution 
N1 = 1000       # maximum iters in optimization
factr_tol = 1e2 
N_sam_DIV = 5   # sample size for ELBO evaluation
N_sam = 100
lmm = 6         # histogram size
mc.cores = parallel::detectCores() - 2
MC = 20    
init_bound = 2
seed_list = 1:MC

model <- stan_model(file)

#######################################################################
## fit with pathfinder ##
fit_pathfinder = FALSE
if(fit_pathfinder){
  
  t <- proc.time()
  opath_bdf <- opt_path_stan_parallel(seed_list, seed_list, mc.cores, model,
                                      stan_data_bdf, init_bound = init_bound, 
                                      N1, N_sam_DIV, N_sam, 
                                      factr_tol, lmm)
  pf_bdf_time <- (proc.time() - t)
  
  print(pf_bdf_time)
  
  pick_samples_bdf <- Imp_Resam_WR(opath_bdf, n_sam = nChains, seed = 1)
  
  # inverse metric estimation #
  pick_ind <- select_normal_apx(opath_bdf)  # return the index of the best normal approximation
  K_hat <- check_pareto_k(opath_bdf[[pick_ind]])   # conduct pareto k hat diagnostic
  K_hat # 0.9244372
  pick_Hk <- est_inv_metric(opath_bdf) # estimate of the mass matrix
  #pick_Hk <- opath_bdf[[pick_ind]]$sample_pkg_save$inv_metric_est # retrieve inverse metric estimate from a specific normal approximation
  
  # generate initials based on Pathfinder #
  posterior <- to_posterior(model, stan_data_bdf)
  init_pf = apply(pick_samples_bdf, 2, 
                  f <- function(x){constrain_pars(posterior, x)})
  
  for (i in 1:nChains) {
    init0 <- init_pf[[i]]
    write_stan_json(init0, 
                    paste0("./init/initpfbdf_pop", i, ".json"))
  }
  
  inv_mass <- list()
  for (i in 1:nChains){
    inv_mass[[i]] <- pick_Hk
  }
  
  save(list = c("pf_bdf_time"), 
       file = paste0("./output/pf_warm_time_pop.RData"))
  
  save(list = c("inv_mass"), 
       file = paste0("./output/pf_warm_inv_mass_pop.RData"))
}


## Fit model with rk45 solver and initials suggested by pathfinder
run_model_pf <- FALSE
saved_fit_pf_file <- paste0("./output/", model_name, ".rk45pf.fit.RDS")
load(paste0("./output/pf_warm_time_pop.RData"))

if (run_model_pf) {
  fit_pf <- mod$sample(
    data = stan_data_rk45, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
    seed = stan_seed, adapt_delta = 0.8,
    init = paste0("./init/initpfbdf_pop", 1:nChains, ".json"),
    save_warmup = TRUE)
  
  fit_pf$cmdstan_diagnose()
  
  fit_pf$save_object(saved_fit_pf_file)
}

fit_pf <- readRDS(saved_fit_pf_file)
fit_pf$time()

## Fit model with rk45 solver and initials and inv-metric suggested by pathfinder

saved_fit_pf_lit_file <- paste0("./output/", model_name, ".rk45pflit.fit.RDS")

load(paste0("./output/pf_warm_inv_mass_pop.RData"))


if (run_model_pf) {
  fit_pf_lit <- mod$sample(
    data = stan_data_rk45, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = 50, iter_sampling = iter_sampling,
    seed = stan_seed, adapt_delta = 0.8,
    init = paste0("./init/initpfbdf_pop",  1:nChains, ".json"),
    inv_metric = inv_mass,
    metric = "dense_e",
    init_buffer = 0,
    term_buffer = 50,
    save_warmup = TRUE)
  
  fit_pf_lit$cmdstan_diagnose()
  
  fit_pf_lit$save_object(saved_fit_pf_lit_file)
}

fit_pf_lit <- readRDS(saved_fit_pf_lit_file)
fit_pf_lit$time()

parms <- c("ka_pop", "V_pop", "Vm_pop", "Km_pop", "sigma", "lp__")
summary(fit_pf$draws(parms))
summary(fit_pf_lit$draws(parms))
mcmc_trace(fit_pf$draws(parms, inc_warmup = TRUE))
mcmc_trace(fit_pf_lit$draws(parms, inc_warmup = TRUE))
mcmc_trace(fit_pf_lit$draws(parms))

#####################################################################
## initial from pf, warm-up with BDF, sample with rk45

saved_fit_pf_bdf_warmup_file <- paste0("./output/", model_name, ".pf_bdf_warmup.fit.RDS")
saved_fit6_file <- paste0("./output/", model_name, ".pf_bdf_rk45.fit.RDS")
if (run_model) {
  fit_pf_bdf_warmup <- mod$sample(
    data = stan_data_bdf, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup, iter_sampling = 1,
    seed = stan_seed,
    init = paste0("./init/initpfbdf_pop", 1:nChains, ".json"))
  
  fit_pf_bdf_warmup$save_object(saved_fit_pf_bdf_warmup_file)
  
  samples_warmup <- fit_pf_bdf_warmup$draws()
  parm_index <- 2:50
  init_files <- save_last_iter_pop(samples_warmup, nChains, 
                                   paste0(model_name, "_w3_pf."),
                                   parm_index)

  mass_matrix <- fit_pf_bdf_warmup$inv_metric()
  step_size <- fit_pf_bdf_warmup$metadata()$step_size_adaptation
  
  fit6 <- mod$sample(
    data = stan_data_rk45, chains = nChains,
    parallel_chains = parallel_chains,
    iter_warmup = 0,
    iter_sampling = iter_sampling,
    seed = stan_seed,
    init = init_files,
    step_size = fit_pf_bdf_warmup$metadata()$step_size_adaptation,
    inv_metric = fit_pf_bdf_warmup$inv_metric(matrix = F),
    # inv_metric = diag(mass_matrix[[1]])
    adapt_engaged = FALSE
  )
  
  fit6$cmdstan_diagnose()  
  fit6$save_object(saved_fit6_file)
}

fit_pf_bdf_warmup <- readRDS(saved_fit_pf_bdf_warmup_file)
fit_pf_bdf_warmup$time()

fit6 <- readRDS(saved_fit6_file)
fit6$time()


###############################################################################
# Central and worse-case scenario metrics don't tell the whole story.
# Let's plot all the points (for one "representative parameter")
# Need now to extract the run time per chain.

parm_index <- 6 # check the relax time for log-density

## Fit model with rk45 solver and initials suggested by pathfinder
time_pf <- (fit_pf$time()$chains[, 4]) + pf_bdf_time[3]
ess_pf <- ess_summary_2(fit = fit_pf, parms, nChains, time_pf)

1/ess_pf$chain_eff[parm_index, ]
# 1.680897 1.324314 1.491202 1.508716 2.182576 1.191092 1.535192 1.734745

## Fit model with rk45 solver and initials and inv-metric suggested by pathfinder
time_pf_lit <- (fit_pf_lit$time()$chains[, 4]) + pf_bdf_time[3]
ess_pf_lit <- ess_summary_2(fit = fit_pf_lit, parms, nChains, time_pf_lit)

1/ess_pf_lit$chain_eff[parm_index, ]
# 0.3603418 0.3168348 0.3135876 0.3463095 0.4112261 0.4114141 0.4400490 0.3892146

## initial from pf, warm-up with BDF, sample with rk45
time_pf_bdf_warmup <- pf_bdf_time[3] + fit_pf_bdf_warmup$time()$chains[, 4]
time_rk45_sampling <- fit6$time()$chains[, 4]
time_pf_bdf_rk45 <- time_pf_bdf_warmup + time_rk45_sampling
ess_pf_bdf_rk45 <- ess_summary_2(fit = fit6, parms, nChains, time_pf_bdf_rk45)
1 / ess_pf_bdf_rk45$chain_eff[parm_index,]
# 1.594780 1.942493 1.447381 1.556838 1.814276 1.398494 1.890933 2.324157

