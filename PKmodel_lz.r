
#rm(list = ls())
gc()
library(ggplot2)
library(cmdstanr)
#set_cmdstan_path("~/Rlib/cmdstan/")
library(rjson)
library(posterior)
library(deSolve)
library(bayesplot)

seed_list <- c(1954, 1234, 1111, 4321, 2222, 3333)
#for(k in 1:length(seed_list)){
  set.seed(seed_list[k]) 
  
  
  #setwd("~/Code/ode-tuning")
  #.libPaths("~/Rlib/")
  
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
    write_stan_json(stan_data, paste0("../example/reslowwarmupodeexample/data/", k,
                    "/PKModel.data.json"))
    
    p <- ggplot(data = data.frame(y = y, times = times),
                aes(y = y, x = times)) +
      geom_point() + theme_bw() + theme(text = element_text(size = 20))
    p
    # plot(x = mass[, 1], y = y)
  }
  
  stan_data_rk45 <- fromJSON(file = paste0("../example/reslowwarmupodeexample/data/", 
                             k, "/PKModel.data.json"))
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
    write_stan_json(init0, paste0("../example/reslowwarmupodeexample/init/", k,
                                  "/init", i, ".json"))
  }
  
  #####################################################################
  ## Tuning parameters for MCMC + compile model
  
  iter_warmup <- 500
  iter_sampling <- 500
  
  ## Compile model
  if (TRUE) {
    file <- paste0("../example/reslowwarmupodeexample/model/", model_name, ".stan")
    mod <- cmdstan_model(file)
  }
  
  #####################################################################
  ## Fit model with rk45 solver.
  
  run_model <- FALSE
  saved_fit0_file <- paste0("../example/reslowwarmupodeexample/output/", k, "/", 
                            model_name, ".rk45.fit.RDS")
  if (run_model) {
    fit0 <- mod$sample(
      data = stan_data_rk45, chains = nChains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      seed = 123, #123 321
      adapt_delta = 0.8, 
      init = paste0("../example/reslowwarmupodeexample/init/", k, "/init",
                    1:nChains, ".json"))
    
    fit0$cmdstan_diagnose()
    
    fit0$save_object(saved_fit0_file)
  }
  
  fit0 <- readRDS(saved_fit0_file)
  fit0$time()
  fit0$metadata()$step_size_adaptation
  
  #####################################################################
  ## Fit model with bdf solver
  
  saved_fit1_file <- paste0("../example/reslowwarmupodeexample/output/", k, "/",
                            model_name, ".bdf.fit.RDS")
  if (run_model) {
    fit <- mod$sample(
      data = stan_data_bdf, chains = nChains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      seed = 123, adapt_delta = 0.8,
      init = paste0("../example/reslowwarmupodeexample/init/", k , "/init", 
                    1:nChains, ".json"))
    
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
    
    write_stan_json(init_saved, 
                    paste0("../example/reslowwarmupodeexample/init/", k, 
                           "/warm_init", i, ".json"))
  }
  
  saved_fit2_file <- paste0("../example/reslowwarmupodeexample/output/", k, "/",
                            model_name, ".rk45_sampling.fit.RDS")
  if (run_model) {
    fit2 <- mod$sample(
      data = stan_data_rk45, chains = nChains,
      parallel_chains = parallel_chains,
      iter_warmup = 0, iter_sampling = iter_sampling,
      seed = 123,
      init = paste0("../example/reslowwarmupodeexample/init/", k, "/warm_init", 
                    1:nChains, ".json"),
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
  
  saved_fit_warmup_file <- paste0("../example/reslowwarmupodeexample/output/",
                                  k, "/", model_name, ".bdf_warmup.fit.RDS")
  saved_fit3_file <- paste0("../example/reslowwarmupodeexample/output/", k, "/", 
                            model_name, ".rk45_cool.fit.RDS")
  if (run_model) {
    fit_warmup <- mod$sample(
      data = stan_data_bdf, chains = nChains,
      parallel_chains = parallel_chains,
      iter_warmup = partial_warmup, iter_sampling = 1,
      seed = 123,
      init = paste0("../example/reslowwarmupodeexample/init/", k, "/init", 
                    1:nChains, ".json"))
    
    fit_warmup$save_object(saved_fit_warmup_file)
    
    samples <- fit_warmup$draws()
    parm_index <- 2:6
    for (i in 1:nChains) {
      init_saved <- list(ka = samples[1, i, parm_index[1]],
                         V = samples[1, i, parm_index[2]],
                         Vm = samples[1, i, parm_index[3]],
                         Km = samples[1, i, parm_index[4]],
                         sigma = samples[1, i, parm_index[5]])
      
      write_stan_json(init_saved, 
                      paste0("../example/reslowwarmupodeexample/init/", k, 
                             "/cool_init", i, ".json"))
    }
    
    mass_matrix <- fit_warmup$inv_metric()
    step_size <- fit_warmup$metadata()$step_size_adaptation
    
    fit3 <- mod$sample(
      data = stan_data_rk45, chains = nChains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup - partial_warmup,
      iter_sampling = iter_sampling,
      seed = 123,
      init = paste0("../example/reslowwarmupodeexample/init/", k, "/cool_init", 
                    1:nChains, ".json"),
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
  
  #######################################################################
  ### warmup with pathfinder ###
  
  source("../utils/sim_variational_pf_mean.R")
  N1 = 300    # Maximum iters in optimization
  N_mode_max = 1
  N_sam = 3   # ELBO
  N_mass = 5
  factr_tol = 1e2
  mc.cores = parallel::detectCores() - 2
  MC = 20    # 10 iterations
  init_bound = 2
  seed_list = 1:MC 
  
  model <- stan_model(file)
  
  fit_pathfinder = FALSE
  
  
  #######################################################################
  ## fit with pathfinder ##
  # set up lmm #
  if(fit_pathfinder){
    posterior <- to_posterior(model, stan_data_bdf)
    D <- get_num_upars(posterior)
    lmm = min(max(D, 5), N1)
    cat("No. pars:", D," lmm in L-BFGS: ", lmm, "\n")
    
    t <- proc.time()
    opath1 <- opt_path_stan_parallel(seed_list, mc.cores, model, stan_data_bdf,
                                     N1, N_mode_max, 
                                     N_sam, N_mass,
                                     init_bound, factr_tol, lmm)
    
    pf_bdf_time <- (proc.time() - t) 
    
    print(pf_bdf_time)
    #0.027   0.139   3.433
    
    pick_mode1 <- filter_mode2(opath1, fit_info)
    pick_mode2 <- filter_mode2(opath1, fit_info_DIV)
    pick_mode <- intersect(pick_mode1, pick_mode2)
    pick_samples1 <- filter_samples(opath1[pick_mode])
    
    init_pf1 = apply(pick_samples1, 2, 
                     f <- function(x){constrain_pars(posterior, x)})
    
    inv_mass <- lapply(opath1[pick_mode], f <- function(x){
      x$sample_pkg_save[[1]]$inv_mass_est
    })
    
    if(length(init_pf1) < nChains){
      sample_ind <- sample(1:length(init_pf1), nChains, replace = TRUE)
      init_pf1 <- init_pf1[sample_ind]
      inv_mass <- inv_mass[sample_ind]
    }else{
      inv_mass <- inv_mass[1:nChains]
    }
    
    for (i in 1:nChains) {
      init0 <- init_pf1[[i]]
      write_stan_json(init0, 
                      paste0("../example/reslowwarmupodeexample/init/", k, 
                             "/initpfbdf", i, ".json"))
    }
    
    save(list = c("pf_bdf_time"), 
         file = paste0("../example/reslowwarmupodeexample/output/", k, 
                       "/pf_warm_time.RData"))
    
    save(list = c("inv_mass"), 
         file = paste0("../example/reslowwarmupodeexample/output/", k, 
                       "/pf_warm_inv_mass.RData"))
  }
  
  ## Fit model with rk45 solver and initials suggested by pathfinder
  
  saved_fit_pf_file <- paste0("../example/reslowwarmupodeexample/output/", k, 
                              "/", model_name, ".rk45pf.fit.RDS")
  load(paste0("../example/reslowwarmupodeexample/output/", k, 
              "/pf_warm_time.RData"))
  
  if (run_model) {
    fit_pf <- mod$sample(
      data = stan_data_rk45, chains = nChains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
      seed = 123, adapt_delta = 0.8,
      init = paste0("../example/reslowwarmupodeexample/init/", k, "/initpfbdf",
                    1:nChains, ".json"),
      save_warmup = TRUE)
    
    fit_pf$cmdstan_diagnose()
    
    fit_pf$save_object(saved_fit_pf_file)
  }
  
  fit_pf <- readRDS(saved_fit_pf_file)
  fit_pf$time()
  
  ## Fit model with rk45 solver and initials and inv mass suggested by pathfinder
  
  saved_fit_pf_lit_file <- paste0("../example/reslowwarmupodeexample/output/", 
                                  k,  "/", model_name, ".rk45pflit.fit.RDS")
  
  load(paste0("../example/reslowwarmupodeexample/output/", k, 
              "/pf_warm_inv_mass.RData"))
  
  
  if (run_model) {
    fit_pf_lit <- mod$sample(
      data = stan_data_rk45, chains = nChains,
      parallel_chains = parallel_chains,
      iter_warmup = 50, iter_sampling = iter_sampling,
      seed = 123, adapt_delta = 0.8,
      init = paste0("../example/reslowwarmupodeexample/init/", 
                    k,  "/initpfbdf",  1:nChains, ".json"),
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
  
  # parms <- c("ka", "V", "Vm", "Km", "sigma", "lp__")
  # summary(fit0$draws(parms))
  # summary(fit_pf$draws(parms))
  # summary(fit_pf_lit$draws(parms))
  # mcmc_trace(fit_pf$draws(parms, inc_warmup = TRUE))
  # mcmc_trace(fit_pf_lit$draws(parms, inc_warmup = TRUE))
  # mcmc_trace(fit_pf2$draws(parms))
  # mcmc_trace(fit_pf_lit$draws(parms))
  # mcmc_trace(fit$draws(parms))
  
  
  #####################################################################
  ## initial from pf, warm-up with dbf, sample with rk45
  ## warmup with rk45.
  
  saved_fit_pf_bdf_warmup_file <- paste0("../example/reslowwarmupodeexample/output/", 
                                         k, "/", model_name, ".pf_bdf_warmup.fit.RDS")
  saved_fit6_file <- paste0("../example/reslowwarmupodeexample/output/", k, "/",
                            model_name, ".pf_bdf_rk45.fit.RDS")
  if (run_model) {
    fit_pf_bdf_warmup <- mod$sample(
      data = stan_data_bdf, chains = nChains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup, iter_sampling = 1,
      seed = 123,
      init = paste0("../example/reslowwarmupodeexample/init/", k, "/initpfbdf", 
                    1:nChains, ".json"))
    
    fit_pf_bdf_warmup$save_object(saved_fit_pf_bdf_warmup_file)
    
    samples <- fit_pf_bdf_warmup$draws()
    parm_index <- 2:6
    for (i in 1:nChains) {
      init_saved <- list(ka = samples[1, i, parm_index[1]],
                         V = samples[1, i, parm_index[2]],
                         Vm = samples[1, i, parm_index[3]],
                         Km = samples[1, i, parm_index[4]],
                         sigma = samples[1, i, parm_index[5]])
      
      write_stan_json(init_saved, 
                      paste0("../example/reslowwarmupodeexample/init/", k, 
                             "/pf_bdf_init", i, ".json"))
    }
    
    mass_matrix <- fit_warmup$inv_metric()
    step_size <- fit_warmup$metadata()$step_size_adaptation
    
    fit6 <- mod$sample(
      data = stan_data_rk45, chains = nChains,
      parallel_chains = parallel_chains,
      iter_warmup = 0,
      iter_sampling = iter_sampling,
      seed = 123,
      init = paste0("../example/reslowwarmupodeexample/init/", k, 
                    "/pf_bdf_init", 1:nChains, ".json"),
      step_size = step_size,
      inv_metric = fit_warmup$inv_metric(matrix = F),
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
  
#}





