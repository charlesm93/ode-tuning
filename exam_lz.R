parms <- c("ka", "V", "Vm", "Km", "sigma", "lp__")
samples0 <- as_draws_df(fit0$draws(variables = parms))
samples1 <- as_draws_df(fit$draws(variables = parms))
samples2 <- as_draws_df(fit2$draws(variables = parms))
samples3 <- as_draws_df(fit3$draws(variables = parms))
samples4 <- as_draws_df(fit_pf$draws(variables = parms))
samples5 <- as_draws_df(fit_pf_lit$draws(variables = parms))
samples6 <- as_draws_df(fit6$draws(variables = parms))

samples_rk45 <- with(samples0, c(ka, V, Vm, Km, sigma, lp__))
samples_bdf <- with(samples1, c(ka, V, Vm, Km, sigma, lp__))
samples_warm <- with(samples2, c(ka, V, Vm, Km, sigma, lp__))
samples_cool <- with(samples3, c(ka, V, Vm, Km, sigma, lp__))
samples_pf <- with(samples4, c(ka, V, Vm, Km, sigma, lp__))
samples_pf_lit <- with(samples5, c(ka, V, Vm, Km, sigma, lp__))
samples_pf_bdf_rk45 <- with(samples6, c(ka, V, Vm, Km, sigma, lp__))

samples_all <- c(samples_rk45, samples_bdf, samples_warm, samples_cool,
                 samples_pf, samples_pf_lit, samples_pf_bdf_rk45)

# parameters <- c("ka", "V", "Vm", "Km", "sigma")
total_samples <- length(samples_rk45)
parm_name <- rep(rep(parms, each = iter_sampling * nChains), 7)
method <- rep(c("rk45", "bdf", "warm", "cool", "pf", "pf_lit", "pf_bdf_rk45"), 
              each = total_samples)

plot_data <- data.frame(samples = samples_all,
                        parm = parm_name,
                        method = method)

plot <- ggplot(data = plot_data,
               aes(x = samples, color = method, fill = method)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  theme_bw() + facet_wrap(~ parm, scale = "free")
plot

jpeg(filename = paste0("../example/reslowwarmupodeexample/pics/hist6.jpeg"),
     width = 600, height = 500, #740, 
     units = "px", pointsize = 16)
plot
dev.off()

#####################################################################
# Examine global efficiency (include lp__) -- meaning the run time
# is taken to be the total, i.e. worst run time.
ess_bulk <- fit0$summary()$ess_bulk[1:6]
eff0 <- ess_bulk / fit0$time()$total # max(fit0$time()$chains[, 4])

ess_bulk1 <- fit$summary()$ess_bulk[1:6]
eff1 <- ess_bulk1 / fit$time()$total # max(fit$time()$chains[, 4])

ess_bulk2 <- fit2$summary()$ess_bulk[1:6]
eff2 <- ess_bulk2 / (max(fit$time()$chains[, 2]) + fit2$time()$total)

ess_bulk3 <- fit3$summary()$ess_bulk[1:6]
eff3 <- ess_bulk3 / (max(fit_warmup$time()$chains[, 2]) + fit3$time()$total)

ess_bulk4 <- fit_pf$summary()$ess_bulk[1:6]
eff4 <- ess_bulk4 / (fit_pf$time()$total + pf_bdf_time[3]) # max(fit0$time()$chains[, 4])

ess_bulk5 <- fit_pf_lit$summary()$ess_bulk[1:6]
eff5 <- ess_bulk5 / (fit_pf_lit$time()$total + pf_bdf_time[3]) # max(fit0$time()$chains[, 4])

ess_bulk6 <- fit6$summary()$ess_bulk[1:6]
eff6 <- ess_bulk6 / (max(fit_pf_bdf_warmup$time()$chains[, 2]) + 
                       pf_bdf_time[3] + fit6$time()$total) # max(fit0$time()$chains[, 4])



eff <- c(eff0, eff1, eff2, eff3, eff4, eff5, eff6)
parm <- rep(parms, 7)
method <- rep(c("rk45", "bdf", "warm_start", "cool_start", "pf", "pf_lit",
                "pf_bdf_rk45"), each = length(parms))
plot_data <- data.frame(eff = eff, parm = parm, method = method)

plot <- ggplot(data = plot_data,
               aes(x = parm, y = eff, fill = method)) +
  geom_bar(stat = "identity", width = 0.3, alpha = 0.8,
           position = "dodge") + 
  theme_bw() + theme(text = element_text(size = 25)) + coord_flip() +
  ylab("ESS / s") + xlab(" ")
plot

jpeg(filename = paste0("../example/reslowwarmupodeexample/pics/ESSperS6.jpeg"),
     width = 600, height = 800, #740, 
     units = "px", pointsize = 12)
plot
dev.off()

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

time4 <- fit_pf$time()$chains[, 4] + pf_bdf_time[3]
ess4 <- ess_summary(fit = fit_pf, parms, nChains, time4)

time5 <- fit_pf_lit$time()$chains[, 4] + pf_bdf_time[3]
ess5 <- ess_summary(fit = fit_pf_lit, parms, nChains, time5)

time6 <- fit_pf_bdf_warmup$time()$chains[, 2] + fit6$time()$chains[, 3] + 
  pf_bdf_time[3]
ess6 <- ess_summary(fit = fit6, parms, nChains, time6)

# Let's examine only parameter (assuming it tells the same story)
parm_index <- 6
mean_ess <- c(ess0$mean_ess[parm_index],
              ess1$mean_ess[parm_index],
              ess2$mean_ess[parm_index], 
              ess3$mean_ess[parm_index],
              ess4$mean_ess[parm_index],
              ess5$mean_ess[parm_index],
              ess6$mean_ess[parm_index])

sd_ess <- c(ess0$sd_ess[parm_index],
            ess1$sd_ess[parm_index],
            ess2$sd_ess[parm_index], 
            ess3$sd_ess[parm_index],
            ess4$sd_ess[parm_index],
            ess5$sd_ess[parm_index],
            ess6$sd_ess[parm_index])

method <- c("rk45", "bdf", "warm start", "room start", "pf", "pf_lit",
            "pf_bdf_rk45")

plot <- ggplot(data = data.frame(method = method, mean_ess = mean_ess),
               aes(x = method, y = mean_ess)) + theme_bw() +
  geom_bar(stat = "identity", width = 0.3, alpha = 0.8,
           position = "dodge") + theme_bw() + 
  theme(text = element_text(size = 25)) +
  geom_errorbar(aes(ymin = mean_ess - sd_ess, ymax = mean_ess + sd_ess))
plot

jpeg(filename = paste0("../example/reslowwarmupodeexample/pics/meanESS6.jpeg"),
     width = 800, height = 600, #740, 
     units = "px", pointsize = 12)
plot
dev.off()

# central metrics are not good summaries. Let's try plotting all the points.

recorded_eff <- c(ess0$chain_eff[parm_index, ],
                  ess1$chain_eff[parm_index, ],
                  ess2$chain_eff[parm_index, ],
                  ess3$chain_eff[parm_index, ],
                  ess4$chain_eff[parm_index, ],
                  ess5$chain_eff[parm_index, ],
                  ess6$chain_eff[parm_index, ])

method_names <- c("rk45", "bdf", "warm start", "cool start", 
                  "pf", "pf_lit", "pf_bdf_rk45")
method <- rep(method_names, each = nChains)
method <- factor(method, levels = method_names)

plot <- ggplot(data = data.frame(method = method, eff = recorded_eff),
               aes(x = method, y = eff)) + theme_bw() +
  theme(text = element_text(size = 25)) +
  geom_point()
plot

jpeg(filename = paste0("../example/reslowwarmupodeexample/pics/eff6.jpeg"),
     width = 800, height = 600, #740, 
     units = "px", pointsize = 12)
plot
dev.off()


# compute the relaxation time
recorded_tau <- c(ess0$chain_ess[parm_index, ] / ess0$chain_eff[parm_index, ],
                  ess1$chain_ess[parm_index, ] / ess1$chain_eff[parm_index, ],
                  ess2$chain_ess[parm_index, ] / ess2$chain_eff[parm_index, ],
                  ess3$chain_ess[parm_index, ] / ess3$chain_eff[parm_index, ],
                  ess4$chain_ess[parm_index, ] / ess4$chain_eff[parm_index, ],
                  ess5$chain_ess[parm_index, ] / ess5$chain_eff[parm_index, ],
                  ess6$chain_ess[parm_index, ] / ess6$chain_eff[parm_index, ])

median_tau <- c(mean(ess0$chain_ess[parm_index, ] / ess0$chain_eff[parm_index, ]),
                mean(ess1$chain_ess[parm_index, ] / ess1$chain_eff[parm_index, ]),
                mean(ess2$chain_ess[parm_index, ] / ess2$chain_eff[parm_index, ]),
                mean(ess3$chain_ess[parm_index, ] / ess3$chain_eff[parm_index, ]),
                mean(ess4$chain_ess[parm_index, ] / ess4$chain_eff[parm_index, ]),
                mean(ess5$chain_ess[parm_index, ] / ess5$chain_eff[parm_index, ]),
                mean(ess6$chain_ess[parm_index, ] / ess6$chain_eff[parm_index, ]))

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
  geom_bar(stat = "identity", position = "dodge",
           aes(x = method[4 * nChains + 1], y = median_tau[5]),
           alpha = 0.01, width = 0.3) +
  geom_bar(stat = "identity", position = "dodge",
           aes(x = method[5 * nChains + 1], y = median_tau[6]),
           alpha = 0.01, width = 0.3) +
  geom_bar(stat = "identity", position = "dodge",
           aes(x = method[6 * nChains + 1], y = median_tau[7]),
           alpha = 0.01, width = 0.3) +
  theme(text = element_text(size = 25))

plot

jpeg(filename = paste0("../example/reslowwarmupodeexample/pics/relax6.jpeg"),
     width = 600, height = 800, #740, 
     units = "px", pointsize = 12)
plot
dev.off()

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
plot_run_time(fit_pf) + ylim(0, 125)
plot_run_time(fit_pf_lit) + ylim(0, 125)
plot_run_time(fit) + ylim(0, 125)
plot_run_time(fit3) + ylim(0, 125)

bayesplot::mcmc_trace

plot_run_time(fit2)
plot_run_time(fit3)




