library(loo)
library(Matrix)
library(matrixStats)
printf <- function(msg, ...) cat(sprintf(msg, ...), "\n")

to_posterior <- function(model, data) {
  sampling(model, data = data, chains = 1, iter = 1, refresh = 0,
           algorithm = "Fixed_param")
}


opt_path <- function(init, fn, gr, 
                     init_bound = 2.0,
                     N1 = 1000,
                     N_sam_DIV = 5, N_sam = 100,
                     factr_tol = 1e2, lmm = 6, seed = 1234) {
  
  #' Run one-path Pathfinder for initialization init 
  #' 
  #' @param init        initial parameter values
  #' @param fn          negative of log-density 
  #' @param gr          gradient function of negative log-density 
  #' @param init_bound  the boundwith of random initials for each dimension
  #' @param N1          maxium number of iterations for L-BFGS (default = 1000)
  #' @param N_sam_DIV   number of samples used in estimating ELBO
  #' @param N_sam       number of samples from the approximating Normal returned
  #'                    will not do extra samples when < N_sam_DIV 
  #' @param factr_tol   the option factr in optim() (default = 1e2)
  #' @param lmm         the option lmm in optim() (default = 6)
  #' @param seed        random seed for running one-path Pathfinder
  #' @return 
  
  
  set.seed(seed)
  D <- length(init)
  # preallocation
  sample_pkg_save <- NA # save the sample info for each mode
  DIV_save <- -Inf # save the divergence related info for each mode
  fn_call = 0     # No. calls to fn
  gr_call = 0     # No. calls to gr
  lgnorms = c()     # log of gradient's norm
  sknorm_ls = c()    # distance of updates
  thetak_ls = c()     # curvature of updates
  
  
  ## retrieve the optimization trajectory from optim() ##
  y <- matrix(NA, nrow = 1, ncol = D + 1) # record the optimization trajectory and log-density
  LBFGS_fail <- TRUE
  N_try <- 0
  while(LBFGS_fail){   # if L-BFGS cannot run with the initial, find another inital
    LBFGS_fail <- FALSE
    if(N_try == 0){
      N_try = N_try + 1
      y[1, 1:D] <- init
    }else{
      print("\n reinitialize \n")
      N_try = N_try + 1
      y[1, 1:D] <- runif(D, -init_bound, init_bound)
    }
    tryCatch(y[1, D + 1] <- -fn(y[1, 1:D]), 
             error = function(e) { LBFGS_fail <<- TRUE})
    if(LBFGS_fail){next}
    lgnorms[1] <- log(sqrt(sum(gr(y[1, 1:D])^2)))
    tryCatch(
      my_data <- capture.output(
        tt <- optim(par = y[1, 1:D],
                    fn = fn,  # negate for maximization
                    gr = gr,
                    method = "L-BFGS-B",
                    control = list(maxit = N1, #reltol = 0.0, #factr_tol,
                                   pgtol = 0.0, #abstol = 0.0, 
                                   factr = factr_tol,
                                   #ndeps = 1e-8 #,
                                   trace = 6, REPORT = 1, lmm = lmm))
      ), 
      error = function(e) { LBFGS_fail <<- TRUE})
  }
  
  # recover the optimization trajectory X and gradient G.
  L = length(my_data); L
  splited_output = unlist(lapply(my_data, f <- function(x){
    strsplit(as.character(x),split = " ")}))
  
  G_ind = which(splited_output == "G") + 2
  Iter = length(G_ind);
  if(tt$convergence == 0){Iter = Iter - 1}
  X = matrix(NA, nrow = Iter + 1, ncol = D)
  G = matrix(NA, nrow = Iter + 1, ncol = D)
  X[1, ] = y[1, 1:D]
  G[1, ] = gr(y[1, 1:D])  # add this since we cannot retrieve the gradient of the initial 
  for(g in 1:Iter){
    X[g + 1, ] = as.numeric(splited_output[(G_ind[g] - D - 2):(G_ind[g] - 3)])
    G[g + 1, ] = as.numeric(splited_output[G_ind[g]:(G_ind[g] + D - 1)])
  }
  
  ### record geometry info 
  Ykt = as.matrix(G[-1, ] - G[-nrow(G), ])
  Skt = as.matrix(X[-1, ] - X[-nrow(X), ])
  Dk = c()
  thetak = c()
  for(i in 1:Iter){
    Dk[i] = sum(Ykt[i, ] * Skt[i, ])
    thetak[i] = sum(Ykt[i, ]^2) / Dk[i]   # curvature checking
  }
  
  sknorm_ls = c(sknorm_ls, sqrt(rowSums(Skt^2)))       
  thetak_ls = c(thetak_ls, thetak)   
  
  # save results
  lgnorm = log(sqrt(rowSums(G^2)))   
  lgnorms = c(lgnorms, lgnorm)     # record log of the norm of gradients
  fn_ls <- apply(X, 1, fn)         # record fn of the optimization trajectory
  y = rbind(y, cbind(X, -fn_ls))   # update y
  fn_call <- tt$counts[1]       # record the calls of fn in L-BFGS
  gr_call <- tt$counts[2]       # record the calls of gr in L-BFGS
  fn_draws <-  c()
  lp_apx_draws <- c()
  sample_draws <- rep(0.0, D)
  
  # estimate DIV for all approximating Gaussians and save results of the one with maximum DIV
  DIV_ls <-  c()
  DIV_fit_pick <- NA
  DIV_max <- -Inf
  sample_pkg_pick <- list()
  E <- rep(1, D)    # initial diag inv Hessian
  t <- proc.time()
  for (iter in 1:Iter){
    #cat(iter, "\t")
    E <- Form_init_Diag(E, Ykt[iter, ], Skt[iter, ]) # initial estimate of diagonal inverse Hessian
    ill_distr = FALSE
    tryCatch(
      sample_pkg <- 
        Form_N_apx(X[iter + 1, ], 
                   Ykt[1:iter, ],  
                   Skt[1:iter, ],
                   Dk[1:iter], 
                   thetak[1:iter], 
                   E,
                   lmm),
      error = function(e) { ill_distr <<- TRUE})
    if(ill_distr){ next }
    if(is.na(sample_pkg[1])){ next }
    DIV_label = "ELBO"; 
    DIV_fit <- est_DIV(sample_pkg, N_sam_DIV, fn, label = DIV_label) #lCDIV #lADIV  #lIKL #ELBO
    fn_call <- fn_call + DIV_fit$fn_calls_DIV
    fn_draws <- c(fn_draws, DIV_fit$fn_draws)                 # record the log-density of all samples
    lp_apx_draws <- c(lp_apx_draws, DIV_fit$lp_approx_draws)  # record the log-density of approximating Normal for all samples
    sample_draws <- cbind(sample_draws, DIV_fit$repeat_draws) # record all samples 
    DIV_ls <- c(DIV_ls, DIV_fit$DIV)
    if(DIV_fit$DIV > DIV_max){
      DIV_max <- DIV_fit$DIV
      DIV_fit_pick <- DIV_fit
      sample_pkg_pick <- sample_pkg
    }
  } 
  proc.time() -t
  #plot(1:length(DIV_ls), DIV_ls, ylim = c(-200, 200))
  sample_pkg_save <- sample_pkg_pick
  DIV_save <- DIV_fit_pick
  
  if(is.null(DIV_ls)){ # if no normal approximation generated, return the set of parameters found by optim
    return(list(sample_pkg_save = sample_pkg_save, 
                DIV_save = DIV_save, y = y, fn_call = fn_call, 
                gr_call = gr_call, lgnorms = lgnorms, 
                sknorm_ls = sknorm_ls, thetak_ls = thetak_ls,
                fn_draws =  fn_ls[length(fn_ls)], 
                sample_draws = y[length(fn_ls), 1:D],
                lp_apx_draws = Inf))
  }
  
  ## Generate upto N_sam samples from the picked approximating Normal ##
  if(N_sam > length(DIV_fit_pick$fn_draws)){
    draws_N_apx <- Sam_N_apx(sample_pkg_save, 
                             N_sam - length(DIV_fit_pick$fn_draws))
    fn_draws_apx <- apply(draws_N_apx$samples, 2, 
                          f <- function(x){
                            ill_fn = FALSE
                            tryCatch(
                              nld <- fn(x),
                              error = function(e) { ill_fn <<- TRUE})
                            ifelse(ill_fn, Inf, nld)
                          })
    ## update the samples in DIV_save ##
    DIV_save$repeat_draws <- cbind(DIV_save$repeat_draws, draws_N_apx$samples)
    DIV_save$lp_approx_draws <- c(DIV_save$lp_approx_draws,
                                  draws_N_apx$lp_apx_draws)
    DIV_save$fn_draws <- c(DIV_save$fn_draws, fn_draws_apx)
    fn_call = fn_call +  N_sam - length(DIV_fit_pick$fn_draws)
  }
  
  return(list(sample_pkg_save = sample_pkg_save, 
              DIV_save = DIV_save, y = y, fn_call = fn_call, 
              gr_call = gr_call, lgnorms = lgnorms, 
              sknorm_ls = sknorm_ls, thetak_ls = thetak_ls,
              fn_draws =  fn_draws, sample_draws = sample_draws,
              lp_apx_draws = lp_apx_draws))
}

opt_path_stan <- function(model, data, N1, N_sam_DIV, N_sam,
                          factr_tol = 1e7, lmm = 6,
                          init_bound = 2, seed = 1234) {
  
  #' Run one-path Pathfinder for specified Stan model and data for
  #' the specified number of iterations using the specified bound on
  #' uniform initialization on the unconstrained scale.  See opt_path()
  #' for a description of output format and algorithm.
  #'
  #' @param model          Stan model (compiled using rstan::stan_model())
  #' @param data           data list for call to rstan::sampling for specified model
  #' @param N1             maxium number of iterations 
  #' @param N_sam_DIV      number of samples used in estimating ELBO
  #' @param N_sam          number of samples from the approximating Normal returned
  #'                       will not do extra samples when < N_sam_DIV 
  #' @param factr_tol      the option factr in optim() (default = 1e2)
  #' @param lmm            the option lmm in optim() (default = 6)
  #' @param init_bound     the boundwith of random initials for each dimension
  #' @param seed           random seed for generating initials
  #' @return 
  
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  set.seed(seed)
  init <- runif(D, -init_bound, init_bound)
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                  gradient = TRUE)[1] 
  gr <- function(theta) -grad_log_prob(posterior, theta, 
                                       adjust_transform = TRUE)
  
  out <- opt_path(init, fn = fn, gr = gr, N1 = N1, N_sam_DIV = N_sam_DIV,
                  N_sam = N_sam, factr_tol = factr_tol,
                  lmm = lmm, seed = seed)
  return(out)
}

opt_path_stan_parallel <- function(seed_init, seed_list, mc.cores, model, data, 
                                   N1, N_sam_DIV, N_sam, init_bound, 
                                   factr_tol, lmm){
  
  #' Run one-path Pathfinder for multiple random initials in parallel
  #'
  #' @param seed_init      array of random seeds for generating random initials
  #' @param seed_list      array of random seeds for running one-path Pathfinders
  #' @param model          Stan model (compiled using rstan::stan_model())
  #' @param data           data list for call to rstan::sampling for specified model
  #' @param N1             maxium number of iterations 
  #' @param N_sam_DIV      number of samples used in estimating ELBO
  #' @param N_sam          number of samples from the approximating Normal returned
  #'                       will not do extra samples when < N_sam_DIV 
  #' @param factr_tol      the option factr in optim() (default = 1e2)
  #' @param lmm            the option lmm in optim() (default = 6)
  #' @param init_bound     the boundwith of random initials for each dimension
  #' @return 
  
  posterior <- to_posterior(model, data)
  D <- get_num_upars(posterior)
  fn <- function(theta) -log_prob(posterior, theta, adjust_transform = TRUE, 
                                  gradient = TRUE)[1]
  gr <- function(theta) -grad_log_prob(posterior, theta, 
                                       adjust_transform = TRUE)
  
  MC = length(seed_init)
  init = c()
  list_ind = c()
  for(i in 1:MC){
    set.seed(seed_init[i])
    init[[i]] <- runif(D, -init_bound, init_bound)
    list_ind[[i]] <- i
  }
  
  out <- mclapply(list_ind, f <- function(x){
    opt_path(init = init[[x]] ,fn = fn, gr = gr, N1 = N1, N_sam_DIV = N_sam_DIV,
             N_sam = N_sam,  factr_tol = factr_tol,
             lmm = lmm, seed = seed_list[x])
  }, mc.cores = mc.cores)
}

Form_init_Diag <- function(E0, Yk, Sk){
  
  #' Form the initial diagonal inverse Hessian in the L-BFGS update 
  #' 
  #' @param E0       initial diagonal inverse Hessian before updated
  #' @param Yk       update in parameters 
  #' @param Skt      update in gradient 
  #' 
  #' @return 
  
  Dk = sum(Yk * Sk)
  if(Dk == 0){
    return(E0)
  }else{
    thetak = sum(Yk^2) / Dk   # curvature checking
    if((Dk <= 0 | abs(thetak) > 1e12)){ #2.2*e^{-16}
      E = E0
    }else{
      a <- (sum(E0 * Yk^2) / Dk)
      E = 1 / (a / E0 + Yk^2 / Dk - a * (Sk / E0)^2 / sum(Sk^2 / E0))
    }
  }
  return(E)
}

Form_N_apx <- function(x_center, Ykt, Skt, Dk, thetak, E, lmm){
  
  #' Returns sampling metrics of the approximating Gaussian given the history
  #' of optimization trajectory 
  #' 
  #' @param x_center The center of the approximating Gaussian
  #' @param Ykt      All history of changes along optimization trajectory
  #' @param Skt      All history of changes of gradients along optimization trajectory
  #' @param E        initial diagonal inverse Hessian
  #' @param lmm      The size of history
  #' 
  #' @return 
  
  D = length(x_center)
  
  Neg_D_ind <- c()
  #curvature condition: abs(thetak) > 1/.Machine$double.eps
  if(any(Dk < 0 | abs(thetak) > 1e12)){   #any(Dk < 0 | abs(thetak) > 1e7)
    print("Negative or unstable Hessian")
    Neg_D_ind = which(Dk < 0 | abs(thetak) > 1e12)
    Ykt = as.matrix(Ykt[-Neg_D_ind, ])
    Skt = as.matrix(Skt[-Neg_D_ind, ])
    Dk = Dk[-Neg_D_ind]
    thetak = thetak[-Neg_D_ind]
  }
  m = length(Dk)
  
  if(m < 2){ # cannot approximate Hessian
    return(NA)
  }else if(m > lmm){ 
    
    # if the size of total history is greater than lmm, only keep upto lmm latest history 
    truc_ind = (1:(m - lmm))
    Ykt = as.matrix(Ykt[-truc_ind, ])
    Skt = as.matrix(Skt[-truc_ind, ])
    Dk = Dk[-truc_ind]
    thetak = thetak[-truc_ind]
  }
  m = min(m, lmm)
  # simulation samples of approximated posterior distribution and estimate E_lp
  
  Rk = matrix(0.0, nrow = m, ncol = m)
  for(s in 1:m){
    for(i in 1:s){
      Rk[i, s] = sum(Skt[i, ] * Ykt[s, ])
    }
  }
  ninvRST = -backsolve(Rk, Skt)
  
  if( 2*m >= D){
    # directly calculate inverse Hessian and the cholesky decomposition
    Hk = diag(E, nrow = D) + 
      crossprod(Ykt %*% diag(E, nrow = D), ninvRST) + 
      crossprod(ninvRST, Ykt %*% diag(E, nrow = D))  + 
      crossprod(ninvRST, 
                (diag(Dk) + 
                   tcrossprod(Ykt %*% diag(sqrt(E), nrow = D))) %*% 
                  ninvRST)
    cholHk = chol(Hk)
    logdetcholHk = determinant(cholHk)$modulus
    
    sample_pkg <- list(label = "full", cholHk = cholHk, 
                       logdetcholHk = logdetcholHk,
                       x_center = x_center)
    
  } else {
    # use equation ?? to sample
    Wkbart = rbind(Ykt %*% diag(sqrt(E)),
                   ninvRST %*% diag(sqrt(1 / E)))
    Mkbar = rbind(cbind(matrix(0.0, nrow = m, ncol = m), diag(m)),
                  cbind(diag(m),
                        (diag(Dk) +
                           tcrossprod(Ykt %*% diag(sqrt(E))))))
    qrW = qr(t(Wkbart))
    Qk = qr.Q(qrW)
    Rkbar = qr.R(qrW)
    Rktilde = chol(Rkbar %*% Mkbar %*% t(Rkbar) + diag(nrow(Rkbar)))
    logdetcholHk = sum(log(diag(Rktilde))) + 0.5 * sum(log(E))
    
    sample_pkg <- list(label = "sparse", theta_D = 1 / E,
                       Qk = Qk, Rktilde = Rktilde,
                       logdetcholHk = logdetcholHk,
                       Mkbar = Mkbar,
                       Wkbart = Wkbart,
                       x_center = x_center)
  }
  return(sample_pkg)
}


est_DIV <- function(sample_pkg, N_sam, fn, label = "ELBO"){
  
  #' estimate divergence based on Monte Carlo samples given the output of 
  #' subroutine Form_N_apx
  #' 
  #' @param sample_pkg  The output of subroutine Form_N_apx
  #' @param N_sam       number of samples used in estimating divergence
  #' @param fn          negative of log density
  #' @param label       type of divergence, (default = "ELBO")
  #' 
  #' @return 
  
  D <- length(sample_pkg$x_center)
  fn_draws <-  rep(Inf, N_sam)
  lp_approx_draws <- rep(0.0, N_sam)
  repeat_draws <- matrix(0, nrow = D, ncol = N_sam)
  draw_ind = 1
  fn_calls_DIV = 0
  f_test_DIV = -Inf
  
  for(l in 1:(2*N_sam)){
    if(sample_pkg$label == "full"){
      u = rnorm(D)
      u2 = crossprod(sample_pkg$cholHk, u) + sample_pkg$x_center
    }else{
      u = rnorm(D)
      u1 = crossprod(sample_pkg$Qk, u)
      u2 = diag(sqrt(1 / sample_pkg$theta_D)) %*%
        (sample_pkg$Qk %*% crossprod(sample_pkg$Rktilde, u1) + 
           (u - sample_pkg$Qk %*% u1)) + sample_pkg$x_center
    }
    # skip bad samples
    skip_flag = FALSE
    tryCatch(f_test_DIV <- fn(u2),
             error = function(e) { skip_flag <<- TRUE})
    if(is.nan(f_test_DIV)){skip_flag <<- TRUE}
    if(skip_flag){
      next
    } else {
      fn_draws[draw_ind] <- f_test_DIV
      lp_approx_draws[draw_ind] <- - sample_pkg$logdetcholHk - 0.5 * sum(u^2)
      repeat_draws[, draw_ind] <- u2
      draw_ind = draw_ind + 1
    }
    fn_calls_DIV = fn_calls_DIV + 1
    if(draw_ind == N_sam + 1){break}
  }
  
  ### Divergence estimation ###
  ELBO <- -mean(fn_draws) - mean(lp_approx_draws)
  if(is.nan(ELBO)){ELBO <- -Inf}
  if(label == "ELBO"){
    DIV <- ELBO
  }else if(label == "lIKL"){
    DIV <- ELBO + log(mean(exp(-fn_draws - lp_approx_draws - ELBO)))    # log Inclusive-KL E[p(x_i)/q(x_k)]
  }else if(label == "lADIV"){
    DIV <- 0.5 * ELBO + 
      log(mean(exp(0.5 * (-fn_draws - lp_approx_draws - ELBO)))) # log alpha-divergence E[(p(x_i)/q(x_k))^alpha], e.g. with 1/2
  }else if(label == "lCDIV"){
    DIV <- 2.0 * ELBO + 
      log(mean(exp(2.0 * (-fn_draws - lp_approx_draws - ELBO))))# log Chi^2-divergence is alpha-divergence with alpha=2
  } else{
    stop("The divergence is misspecified")
  }
  
  if(is.nan(DIV)){DIV <- -Inf}
  if(is.infinite(DIV)){DIV <- -Inf}
  
  if(any(is.nan(fn_draws))){
    fn_draws[is.nan(fn_draws)] <- Inf
  }
  
  return(list(DIV = DIV, 
              repeat_draws = repeat_draws, fn_draws = fn_draws, 
              lp_approx_draws = lp_approx_draws, fn_calls_DIV = fn_calls_DIV))
}

Sam_N_apx <- function(sample_pkg, N_sam){
  
  #' Generate N_sam samples from the approximating Normal given the output of 
  #' subroutine Form_N_apx
  #' 
  #' @param sample_pkg The output of subroutine Form_N_apx
  #' @param N_sam      Number of samples 
  #' 
  #' @return \item{samples}{The samples from the approximating Normal, 
  #'   each column stores a sample} 
  #'   \item{lp_approx_draws}{The log-density of generated samples under 
  #'   approximating Normal} 
  
  
  D <- length(sample_pkg$x_center)
  lp_approx_draws <- rep(0.0, N_sam)
  repeat_draws <- matrix(0, nrow = D, ncol = N_sam)
  
  if(sample_pkg$label == "full"){
    u = matrix(rnorm(D * N_sam), nrow = D)
    u2 = crossprod(sample_pkg$cholHk, u) + sample_pkg$x_center
  }else{
    u = matrix(rnorm(D * N_sam), nrow = D)
    u1 = crossprod(sample_pkg$Qk, u)
    u2 = diag(sqrt(1 / sample_pkg$theta_D)) %*%
      (sample_pkg$Qk %*% crossprod(sample_pkg$Rktilde, u1) + 
         (u - sample_pkg$Qk %*% u1)) + sample_pkg$x_center
  }
  
  lp_apx_draws <- - sample_pkg$logdetcholHk - 0.5 * colSums(u^2)
  
  return(list(samples = u2, lp_apx_draws = lp_apx_draws))
}


params_only <- function(path) {
  
  #' Return optimization path with last column (objective function value)
  #' removed.
  #'
  #' @param path   optimization path with last column for objective
  #'               function value
  #' @return
  
  N <- dim(path)[1]
  D <- dim(path)[2]
  path[1:N, 1:(D - 1)]
}

fit_info_DIV <- function(param_path) {
  # extract fitting information #
  DIV = param_path$DIV_save$DIV
  
  return(DIV)
}

Imp_Resam_WOR <- function(param_path, n_inits, seed = 123){
  
  #' SIR without replacement, index of apporoximating distribution as an auxilary variable
  #' Return n_inits distinct samples
  #'  
  #' @param param_path output of function opt_path_stan_parallel()
  #' @param n_inits    number of distinct samples
  #' @param seed       random seed for importance resampling
  
  # remove the failed Pathfinder
  check <- sapply(param_path, f <- function(x){
    work <- is.finite(x$lp_apx_draws[[1]])
    work
  })
  
  filter_mode <- c(1:length(check))[check]
  
  param_path <- param_path[filter_mode]
  
  ## extract samples and log ratios ##
  samples <- do.call("cbind", lapply(param_path, extract_samples))
  lrms <- c(sapply(param_path, extract_log_ratio))
  
  ## take off samples with infinite log ratios
  finit_ind <- is.finite(lrms)
  samples <- samples[, finit_ind]
  lrms <- lrms[finit_ind]
  
  ## compute the importance weight ##
  sample_weights_psis <- suppressWarnings(weights(psis(lrms, r_eff = 1),
                                                  log = FALSE))
  #sample_weights_IS <- exp(lrms - max(lrms))/sum(exp(lrms - max(lrms)))
  
  ## Importance resampling ##
  set.seed(seed)
  sample_ind <-sample.int(ncol(samples), size = n_inits, replace = FALSE, 
                          prob = sample_weights_psis)
  
  return(samples[, sample_ind])
}

extract_samples <- function(param_path){
  sample <- param_path$DIV_save$repeat_draws
  return(sample)
}

extract_samples_center <- function(param_path){
  sample <- param_path$sample_pkg_save$x_center
  return(sample)
}

extract_lps <- function(param_path){
  lps <- -param_path$DIV_save$fn_draws
  return(lps)
}

constrain_draws <- function(unconstrain_draws, posterior){
  lapply(1:ncol(unconstrain_draws), 
         f <- function(s){
           constrain_pars(posterior, unconstrain_draws[, s])
         })
}

extract_log_ratio <- function(param_path){
  lrm <- - param_path$DIV_save$fn_draws -
    param_path$DIV_save$lp_approx_draws
  return(lrm)
}


#' pick_inv_mass <- function(param_path){
#'   
#'   #' Extract the inverse Hessian approximation of the normal approximation
#'   #' with the highest ELBO 
#'   #' param_path: output of function opt_path_stan_parallel()
#' 
#'   # remove the failed Pathfinder
#'   check <- sapply(param_path, f <- function(x){
#'     work <- is.finite(x$lp_apx_draws[[1]])
#'     work
#'   })
#'   
#'   filter_mode <- c(1:length(check))[check]
#'   param_path <- param_path[filter_mode]
#'   
#'   # find the index of the normal approximation with the higest ELBO 
#'   ELBO <- apply(sapply(param_path, extract_log_ratio), 2, mean)
#'   pick_index <- which.max(ELBO)
#'   
#'   return(param_path[[pick_index]]$sample_pkg_save$inv_mass_est)
#'   
#' }


