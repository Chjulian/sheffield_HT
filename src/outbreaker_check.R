## load packages and data
if(!require(pacman)) install.packages(pacman)
pacman::p_load(here, tidyverse, rio, magrittr, remotes, epicontacts, glue)
devtools::load_all("~/Dropbox/dev/outbreaker2")
load(here("data/eps0.999.RData"))
source(here("R/functions.R"))

data <- d_fix0.999$data
config <- d_fix0.999$config
data <- add_convolutions(data, config)
res <- d_fix0.999$res

## ## set up config
## config <- create_config(
##   ## define iterations and sampling (this is a short exploratory run)
##   n_iter = 1e4, sample_every = 50,
##   ## prior on the reporting probability pi
##   move_pi = TRUE, prior_pi = c(5, 5),
##   ## eps is the probability of infections occuring between cases registered on the same ward
##   move_eps = TRUE, prior_eps = c(1, 1),
##   ## tau is the probability of an unobserved cases being moved to a different ward
##   move_tau = TRUE, prior_tau = c(1, 1),
##   ## leave these options as is for now
##   move_joint = TRUE, move_model = TRUE
## )

## ## identify initial tree (run this for now, it will identify a possible starting point)
## config <- get_initial_tree(data, config, n_iter = 1e3, max_dist = 2)

## ## run outbreaker
## res <- outbreaker(data, config)

## burnin
res <- res[20:nrow(res),]

## get param object
param <- create_param(data, config)$current

## update param object
param <- update_par(param, nrow(res), res, data, config)

## visualise results
vis_ward(res, meta, wards)

## extract tpair data with burnin
tpairs <- get_tpairs(res, data)

length(unique(paste(tpairs$from,tpairs$to,sep="_")))
pp <- as.numeric(table(paste(tpairs$from,tpairs$to, sep="_"))/nrow(res))
hist(pp);summary(pp)

## calculate proportion of ancestries that are on a different ward
## we see that only ~1% of links are between different, known wards
## we do see that ~40% of links are from an unknown -> known ward
tpairs %>%
  filter(kappa == 1) %>%
  summarise(
    ## these are both on a known ward which is different
    diff_ward_known = mean(
      ward_from != ward_to &
      ward_from != 0 &
      ward_to != 0
    ),
    ## wards are different but one is unknown
    diff_ward_unknown = mean(
      ward_from != ward_to &
      (ward_from == 0 | ward_to == 0)
    ),
    ## on the same ward, both known
    same_ward_known = mean(
      ward_from == ward_to &
      ward_from != 0
    ),
    ## both on an unknown ward
    same_ward_unknown = mean(
      ward_from == ward_to &
      ward_from == 0
    )
  )

## find links between different, known wards that are accepted more than 10% of
## the time (i.e. regularly)
diff_ward_known <- filter(
  tpairs,
  kappa == 1,
  ward_from != ward_to,
  ward_from != 0,
  ward_to != 0
) %>%
  group_by(to) %>%
  summarise(prop = n()/nrow(res)) %>%
  filter(prop > 0.1) %>%
  arrange(desc(prop))

## search for potential ancestors on the same ward for these "wrong" ancestries
potential_known <- lapply(unique(diff_ward_known$to), get_potential, data)

## we see there are no plausible ancestors for 95% of cases for which outbreaker
## has assigned "wrong" ancestries between wards. the remaining wrong links are
## likely due to the epi data constraining the tree in a manner that forces
## these links to be accepted.
mean(sapply(potential_known, nrow) == 0)

## find links between unknown -> known or unknown -> unknown wards that are
## accepted more than 10% of the time (i.e. regularly)
diff_ward_unknown <- filter(
  tpairs,
  kappa == 1 &
  (
    (ward_from != ward_to & (ward_from == 0 | ward_to == 0)) |
    (ward_from == ward_to & ward_from == 0)
  )
) %>%
  group_by(from, to) %>%
  summarise(prop = n()/nrow(res)) %>%
  group_by(to) %>%
  mutate(total_prop = sum(prop)) %>%
  filter(total_prop > 0.1) %>%
  arrange(desc(total_prop))

## search for potential ancestors on the same ward for these "wrong" ancestries
potential_unknown <- unique(diff_ward_unknown$to) %>%
  setNames(., .) %>%
  lapply(get_potential, data)

## we can see that for 82% of cases for which an unknown -> known link is
## assigned, this is because there is no possible ancestor on the same ward
mean(sapply(potential_unknown, nrow) == 0)

## for "wrong" ancestries with potentially correct ones, check to see if they
## have been explored in the MCMC, or check if the potentially correct one is
## actually the infectee already (i.e. it's a loop where only one can be correct)
potential_check <- names(potential_unknown)[sapply(potential_unknown, nrow) > 0] %>%
  setNames(., .) %>%
  lapply(
    function(i) {
      potents <- unique(potential_unknown[[i]]$from)
      tibble(
        from = potents,
        to = i,
        ## this checks if the proposed ancestor has been accepted as ancestor in the MCMC
        already_ancestor = sapply(
          potents,
          function(potent) sum(res[[glue("alpha_{i}")]] == potent, na.rm = TRUE)/nrow(res)
        ),
        ## this checks if the proposed ancestor has been accepted as the infectee in the MCMC
        ancestor_is_infectee = sapply(
          potents,
          function(potent) sum(res[[glue("alpha_{potent}")]] == i, na.rm = TRUE)/nrow(res)
        ),
        either = already_ancestor + ancestor_is_infectee
      )
    }
  ) %>%
  bind_rows()
