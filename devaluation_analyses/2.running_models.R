# running directed action


# load the packages
library(rstan)
library(tidyverse)
library(R.matlab)
library(bayesplot)
library(wesanderson)
library(waterfalls)

# load the data
load('modelling/datalists/online_patient_datalist.RData')

# running stan
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

modelFile <- 'modelling/models/DA_model_5.stan'

nIter     <- 2000
nChains   <- 4 
nWarmup   <- floor(nIter/2)
nThin     <- 1

fit_rl <- stan(modelFile, 
               data    = datalist, 
               chains  = nChains,
               iter    = nIter,
               warmup  = nWarmup,
               thin    = nThin,
               init    = "random",
               seed    = 3210
)

# checking diagnostics with plots
plot_dens_tau <- stan_plot(fit_rl, pars=c('mu_tau','sigma_tau','tau'), show_density=T, fill_color = 'skyblue')
plot_dens_lr <- stan_plot(fit_rl, pars=c('a_lr', 'b_lr', 'lr'), show_density=T, fill_color = 'skyblue')
plot_trace <- stan_trace(fit_rl, pars=c('k_tau','theta_tau', 'a_lr', 'b_lr'), inc_warmup = F)

# save model and summary table
save(fit_rl, file='modelling/fitted_models/model_1.RData')
