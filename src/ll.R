## load packages and data
if(!require(pacman)) install.packages(pacman)
pacman::p_load(lubridate, tidyverse, rio, magrittr, remotes, epicontacts, glue, here)
remotes::install_github("finlaycampbell/outbreaker2@updated_hosp")
load(here("data/eps0.999.RData"))
source(here("R/functions.R"))

## load data
data <- d_fix0.999$data
config <- d_fix0.999$config
data <- outbreaker2:::add_convolutions(data, config)
res <- d_fix0.999$res
wards <- d_fix0.999$wards
meta <- tibble(
  id = rownames(data$dna),
  date_onset = d_fix0.999$raw$dateofonset_foranalysis
)

## beta prior on mu (using parameters calculated above)
prior_mu <- function(param) {
  dbeta(param$mu, 89.9, 809.1, log = TRUE)
}

## poisson ll on genetic distance
ll_genetic <- function(data, param, i = NULL) {
  if(is.null(i)) i <- seq_along(param$alpha)
  sum(
    dpois(
      data$D[cbind(param$alpha[i], i)],
      param$mu,
      log = TRUE
    ),
    na.rm = TRUE
  )
}

## create ll and prior objects
priors <- custom_priors(mu = prior_mu)
likelihoods <- custom_likelihoods(genetic = ll_genetic)

## set up config
config <- create_config(
  ## define iterations and sampling (this is a short exploratory run)
  n_iter = 1e3, sample_every = 50,
  ## prior on the reporting probability pi
  move_pi = TRUE, prior_pi = c(5, 5),
  ## eps is the probability of infections occuring between cases registered on the same ward
  move_eps = TRUE, prior_eps = c(1, 1),
  ## tau is the probability of an unobserved cases being moved to a different ward
  move_tau = TRUE, prior_tau = c(1, 1),
  ## set initial mu values
  move_mu = TRUE, init_mu = 0.1,
  ## leave these options as is for now
  move_joint = TRUE, move_model = TRUE, move_kappa = FALSE,
  ## increase sd of proposal for infection times
  sd_t_inf = 6, sd_t_onw = 15,
  ## increase sd of proposal for mu
  init_mu = 0.2, sd_mu = 0.05
)

## identify initial tree
config <- get_initial_tree(
  data = data,
  config = config,
  priors = priors,
  likelihoods = likelihoods,
  n_iter = 1e2,
  max_dist = 1e3,
  days_after_admission = 2,
  scale = 1
)

## run outbreaker once for likelihood calculations
pre_run <- outbreaker(data, config, priors = priors, likelihoods = likelihoods)

## make empty param
param <- create_param(data)$current

## get the average likelihood per case (make sure to set the correct burnin)
ll <- get_average_likelihood(pre_run, data, param, config, likelihoods, burnin = 500)

## assign imports for 100 runs with an 60% import probability
configs <- assign_imports(n_chains = 10, p_import_mu = 0.6, p_import_sd = 0.1, res, ll, config)

## run and save outbreaker (in outputs folder) and return as list
results <- imap(configs, run_and_save, data, priors, likelihoods)
