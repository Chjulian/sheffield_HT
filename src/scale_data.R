if(!require(pacman)) install.packages(pacman)
pacman::p_load(lubridate, tidyverse, rio, magrittr, remotes, epicontacts, glue,
               here)
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

## load raw ward data
wards_raw <- import(here("data/patient-ward-movements-20210207.csv")) %>%
  rename_all(tolower) %>%
  rename(start = `in`, end = out, id = barcode) %>%
  mutate(
    ## convert these to POSIXct dates
    across(
      c(start, end),
      as.POSIXct,
      format = "%Y-%m-%d %H:%M:%S",
      tz = "GMT"
    ),
    across(
      c(id, ward),
      ~ tolower(str_replace_all(.x, "-", "_"))
    )
  ) %>%
  drop_na(start, end) %>%
  arrange(id, start, end)

## define minimum date
origin <- floor_date(min(wards_raw$start, na.rm = TRUE), "day")

## define number of sections to break one day into
scale <- 6

## scale the dates up by a factor of 6 - note that the actual dates won't be
## correct anymore, but you can scale them back down to their correct, rounded
## values using function below. You want to use this object for the outbreaker run.
wards_scaled <- scale_ward_dates(wards_raw, origin, to = "up", scale = scale)

## convert back to normal dates - you can see that the times have been rounded to
## 1/6th days i.e. 4 hours. You DON'T use this for the outbreaker run.
wards_unscaled <- scale_ward_dates(wards_scaled, origin, to = "down", scale = scale)

## scale the onset dates
onset <- d_fix0.999$raw$dateofonset_foranalysis
scaled_onset <- scale_onset_dates(onset, origin, scale = scale)
scaled_onset %<>% setNames(rep("local", length(.)))

## scale the incubation period and generation times
f_dens_scaled <- apply(d_fix0.999$data$f_dens, 2, scale_distribution, scale = scale)
w_dens_scaled <- scale_distribution(d_fix0.999$data$w_dens, scale = scale)

## created scaled data object

## THIS WILL NOT WORK as the dna / id / wards don't quite match up with the ward
## data I have, but the script should work in your pipeline hopefully.
data_scaled <- outbreaker_data(
  dates = scaled_onset,
  dna = dna,
  ctd_timed = wards_scaled,
  ids = meta$id,
  w_dens = w_dens_scaled,
  f_dens = f_dens_scaled,
  p_trans = p_trans
)

## set up config
config <- create_config(
  ## define iterations and sampling (this is a short exploratory run)
  n_iter = 1e4, sample_every = 50,
  ## prior on the reporting probability pi
  move_pi = TRUE, prior_pi = c(5, 5),
  ## eps is the probability of infections occuring between cases registered on the same ward
  move_eps = TRUE, prior_eps = c(1, 1),
  ## tau is the probability of an unobserved cases being moved to a different ward
  move_tau = TRUE, prior_tau = c(1, 1),
  ## set initial mu values
  move_mu = TRUE, init_mu = 0.1,
  ## leave these options as is for now
  move_joint = TRUE, move_model = TRUE,
  ## increase sd of proposal for infection times
  sd_t_inf = 6, sd_t_onw = 15
)

## get beta prior values from mean and sd
get_beta_par(0.98, 0.001)

## beta prior on mu (using parameters calculated above)
prior_mu <- function(param) {
  dbeta(param$mu, 89.9, 809.1, log = TRUE)
}

## poisson ll on genetic distance
ll_genetic <- function(data, param) {
  sum(
    dpois(
      data$D[cbind(param$alpha, seq_along(param$alpha))],
      param$mu,
      log = TRUE
    ),
    na.rm = TRUE
  )
}

## create ll and prior objects
priors <- custom_priors(mu = prior_mu)
likelihoods <- custom_likelihoods(genetic = ll_genetic)

## identify initial tree (run this for now, it will identify a possible starting point)
config <- get_initial_tree(
  data = data,
  config = config,
  priors = priors,
  likelihoods = likelihoods,
  n_iter = 1e2,
  max_dist = 2,
  scale = scale
)

## run outbreaker
res <- outbreaker(data, config, prior = priors, likelihoods = likelihoods)
