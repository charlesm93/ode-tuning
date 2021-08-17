
.libPaths("~/Rlib")
setwd("/Users/charlesm/Code/ode-tuning")

library(ggplot2)

model_name <- "Michaelis_MentenPK"
plot_data_single_hmc <- read.csv(paste0("plot_data/", model_name, ".data.csv"))
plot_data_single_pf <- read.csv(paste0("plot_data/", model_name, ".pf.data.csv"))
plot_data_single <- rbind(plot_data_single_hmc, plot_data_single_pf)

median_tau_single <- read.csv(paste0("plot_data/", model_name,"_tau_median",
                              ".data.csv"))$x

median_tau_single <- c(median_tau_single,
                       median(plot_data_single$tau[plot_data_single$method == "Pf, RK45"]),
                       median(plot_data_single$tau[plot_data_single$method == "Pf, approx tuning"]),
                       median(plot_data_single$tau[plot_data_single$method == "Pf, late switch"])
                       )

max_tau_single <- c(max(plot_data_single$tau[plot_data_single$method == "RK45"]),
                    max(plot_data_single$tau[plot_data_single$method == "BDF"]),
                    max(plot_data_single$tau[plot_data_single$method == "Late switch"]),
                    max(plot_data_single$tau[plot_data_single$method == "Early switch"]),
                    max(plot_data_single$tau[plot_data_single$method == "Pf, RK45"]),
                    max(plot_data_single$tau[plot_data_single$method == "Pf, approx tuning"]),
                    max(plot_data_single$tau[plot_data_single$method == "Pf, late switch"])
                    )

model_name <- "Michaelis_MentenPK_pop_centered"
plot_data_pop_hmc <- read.csv(paste0("plot_data/", model_name, ".data.csv"))
plot_data_pop_pf <- read.csv(paste0("plot_data/", model_name, ".pf.data.csv"))
plot_data_pop <- rbind(plot_data_pop_hmc, plot_data_pop_pf)


median_tau_pop <- read.csv(paste0("plot_data/", model_name,"_tau_median",
                                     ".data.csv"))$x

median_tau_pop <- c(median_tau_pop,
                    median(plot_data_pop$tau[plot_data_single$method == "Pf, RK45"]),
                    median(plot_data_pop$tau[plot_data_single$method == "Pf, approx tuning"]),
                    median(plot_data_pop$tau[plot_data_single$method == "Pf, late switch"])
)

max_tau_pop <- c(max(plot_data_pop$tau[plot_data_pop$method == "RK45"]),
                 max(plot_data_pop$tau[plot_data_pop$method == "BDF"]),
                 max(plot_data_pop$tau[plot_data_pop$method == "Late switch"]),
                 max(plot_data_pop$tau[plot_data_pop$method == "Early switch"]),
                 max(plot_data_pop$tau[plot_data_pop$method == "Pf, RK45"]),
                 max(plot_data_pop$tau[plot_data_pop$method == "Pf, approx tuning"]),
                 max(plot_data_pop$tau[plot_data_pop$method == "Pf, late switch"]))

plot_data_all <- rbind(plot_data_single, plot_data_pop)
plot_data_all$model <- rep(factor(c("Single", "Population"), 
                                  levels = c("Single", "Population")),
                           each = dim(plot_data_single)[1])
plot_data_all$method <- factor(plot_data_all$method, 
                levels = c("RK45", "BDF", "Late switch", "Early switch",
                           "Pf, RK45", "Pf, approx tuning", "Pf, late switch"))

plot_data_all$median_tau <- c(rep(median_tau_single, each = 8),
                              rep(median_tau_pop, each = 8))
plot_data_all$max_tau <- c(rep(max_tau_single, each = 8),
                           rep(max_tau_pop, each = 8))

# To capture outliers, want to use plot on the log scale.
plot <- ggplot(data = plot_data_all,
               aes(x = method, y = tau)) + theme_bw() +
  geom_point() + scale_y_continuous(trans = 'log10') +
  facet_wrap(~model, scale = "free", ncol = 1)  + 
  ylab("Relaxation time (s)") + xlab("") + coord_flip() +
  theme(text = element_text(size = 16)) +
  geom_point(aes(x = method, y = median_tau), 
             shape = 13,
             size = 3, color = "orange") +
  geom_point(aes(x = method, y = max_tau), 
             shape = 21,
             size = 3, color = "red") 

plot

#   geom_point(aes(x = method[1], y = median_tau[1]), shape = 13,
#              size = 3, color = "orange") +
#   geom_point(aes(x = method[nChains + 1], y = median_tau[2]), shape = 13,
#              size = 3, color = "orange") +
#   geom_point(aes(x = method[2 * nChains + 1], y = median_tau[3]), shape = 13,
#              size = 3, color = "orange") +
#   geom_point(aes(x = method[3 * nChains + 1], y = median_tau[4]), shape = 13,
#              size = 3, color = "orange") +
#   ylab("Relaxation time (s)") + xlab("") + coord_flip() +
#   theme(text = element_text(size = 16))
# plot






