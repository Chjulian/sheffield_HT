## load packages and data
if(!require(pacman)) install.packages(pacman)
pacman::p_load(ape, linelist, epitrix, fitdistrplus,
               outbreaker2, tidyverse, rio,
               magrittr, remotes, epicontacts,
               visNetwork, lubridate, pairsnp,
               ggplot2, chron,glue, here)

source("src/functions_wards.R")
source("src/functions.R")
source("src/onset_distribution.R")
source("src/functions_analyse.R")

scale <- 6
reporting.probability <- 0.5
number.import.chains <- 5
output.name <- "w1_i0.7_pi0.5"

## make unique file output name.
output.name <- paste0(output.name, "_",
                      format(Sys.time(), "%m%d%H%M%S"), "_",
                      sample(1:99999, 1))

#################################
#load linelist data
#################################
# load data, clean dates and filter-out observations
mydata <- read.csv("data/hoci-phylo-metadata-onset-2021-02-07.csv",
                   stringsAsFactors =F)
mydata$HealthcareAssociation <- recode(mydata$HealthcareAssociation, "Hospital Onset-Indetermite Healthcare Associated" = "Hospital Onset-Indeterminate Healthcare Associated")
mydata$DateOfCollection <- as.Date(mydata$DateOfCollection, format="%d/%m/%y")
mydata$DateOfOnset <- as.Date(mydata$DateOfOnset, format="%d/%m/%y")
mydata$DateOfAdmission <- as.Date(mydata$DateOfAdmission, format="%d/%m/%y")
mydata <- onsetdistribution(mydata)
mydata<- clean_data(mydata, guess_dates = which(names(mydata)=='DateOfOnset_forAnalysis'))
mydata <- mydata[mydata$category%in%c('inpatient', 'outpatient','staff'),]
mydata <- mydata[mydata$loc1%!in%c('com'),] #rm community associated
mydata <- mydata[!is.na(mydata$dateofonset_foranalysis),]

#################################
#Consolidate ward data
#################################
wards.patients<- read.csv("data/patient-ward-movements-20210207.csv")
wards.patients<- wards.patients[,c('Barcode','Ward','In','Out')]
names(wards.patients) <- c('id','ward','adm','dis') #rename as examples files for consistency
wards<- ward_occupation(wards.patients=wards.patients)
#set nncl and multiple as NA
wards$ward[wards$ward%in%c('nncl','multiple')] <- NA

#################################
#remove barcodes with missing location data
#################################
barcodes.with.movements <- unique(wards$id[!is.na(wards$ward)])
wards <- wards[wards$id %in% barcodes.with.movements,]
mydata <- mydata[mydata$barcode %in% barcodes.with.movements,]

#################################
#process onset dates
#################################
date_onset<- mydata$dateofonset_foranalysis
names(date_onset) <- rep("local", length(date_onset))

#################################
#load dna data in DNAbin format,
#################################
#one sequence per individual
dna <- read.FASTA("data/hoci_phylo-2021-02-07_error-removed_masked.aln")
mylabels <- clean_data(as.data.frame(names(dna)))
names(dna) <- mylabels$names_dna
dna <- dna[mydata$barcode]; mylabels <- names(dna)
write.FASTA(dna, paste0("data/", output.name,  "error-removed_masked.ord.aln"))
dna.SNP <- import_fasta_sparse(paste0("data/", output.name,  "error-removed_masked.ord.aln")) %>% snp_dist(.)
rownames(dna.SNP) <- colnames(dna.SNP) <- mylabels; rm(mylabels)
unlink(paste0("data/", output.name,  "error-removed_masked.ord.aln"))


#################################
# scale dates
#################################
## reformat and set correct time zone
wards %<>%
  rename(start = adm, end = dis) %>%
  drop_na(start, end) %>%
  arrange(id, start, end) %>%
  mutate(across(c(start, end), with_tz, "GMT"))

## correct time zone
date_onset %<>% with_tz("GMT")

## define minimum date
origin <- floor_date(min(wards$start, na.rm = TRUE), "day")

## scale the dates up by a factor of 6 - note that the actual dates won't be
## correct anymore, but you can scale them back down to their correct, rounded
## values using function below. You want to use this object for the outbreaker run.
wards_scaled <- scale_ward_dates(wards, origin, to = "up", scale = scale)

## convert back to normal dates - you can see that the times have been rounded to
## 1/6th days i.e. 4 hours. You DON'T use this for the outbreaker run.
wards_unscaled <- scale_ward_dates(wards_scaled, origin, to = "down", scale = scale)

## scale the onset dates
scaled_date_onset <- date_onset %>%
  scale_onset_dates(origin, scale = scale)

names(scaled_date_onset) <- rep("local", length(scaled_date_onset))
wards_scaled <- wards_scaled %>% rename(adm =  start, dis = end)


#################################
#set incubation_period: a vector indexed at day = 1
#################################
n <- rlnorm(10000, meanlog = 1.621, sdlog = 0.418)
fit.gamma <- fitdist(n, distr = "gamma", method = "mle")
# build discretized gamma
incubation_period <- distcrete::distcrete("gamma",
                                          shape = fit.gamma$estimate['shape'],
                                          scale = fit.gamma$estimate['rate'],
                                          w = 0.5, interval = 1)
rm(n,fit.gamma)

#################################
#set generation_time: a vector indexed at day = 1
#################################
# build discretized gamma
generation_time <- distcrete::distcrete("gamma",
                                        shape = 4.83,
                                        scale = 1.051851,
                                        w = 0.5, interval = 1)
#We can specify multiple incubation periods but we using only one
f <- matrix(
  c(incubation_period$d(1:50), incubation_period$d(1:50)), ncol = 2,
  dimnames = list(NULL, c("local", "import"))
)


## scale the incubation period and generation time
w_dens_scaled <- scale_distribution(generation_time$d(1:50), scale = scale)
f_dens_scaled <- apply(f, 2, scale_distribution, scale = scale)


#################################
#prepare for run
#################################
## define a random set of transition matrices between wards
unq <- unique(wards$ward)[!is.na(unique(wards$ward))]
n <- length(unq)
transition <- matrix(runif(n*n), n, n, dimnames = list(unq, unq))
diag(transition) <- 0
transition <- t(apply(transition, 1, function(x) x/sum(x)))
rm(unq,n)

#create meta object
meta <- mydata[,c("barcode","dateofonset_foranalysis")]
names(meta) <- c("id", "date_onset")
row.names(meta) <- NULL


#################################
#set outbreaker run
#################################
data <- outbreaker_data(
  dates = scaled_date_onset,
  dna = dna,
  ctd_timed = wards_scaled,
  ids = meta$id,
  w_dens = w_dens_scaled,
  f_dens = f_dens_scaled,
  p_trans = list(transition)
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
  n_iter = 1100, sample_every = 50,
  ## prior on the reporting probability pi
  move_pi = FALSE, init_pi=reporting.probability,
  ## eps is the probability of infections occuring between cases registered on the same ward
  move_eps = FALSE, init_eps=0.99,
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
                      
saveRDS(pre_run, paste0("output/", output.name, "_prerun.rds"))

## make empty param
param <- create_param(data)$current

## we need this step for the likelihood calculations
data %<>% outbreaker2:::add_convolutions(config)

## get the average likelihood per case (make sure to set the correct burnin)
ll <- get_average_likelihood(pre_run, data, param, config, likelihoods, burnin = 500)

## assign imports for 100 runs with an 70% import probability
configs <- assign_imports(n_chains = 5, p_import_mu = 0.7, p_import_sd = 0.075, pre_run, ll, config)

## run and save outbreaker (in outputs folder) and return as list
results <- imap(configs, run_and_save, data, priors, likelihoods, output.name)

#################################
#summarise and output
#################################

##read in bay links
baylink.df <- read.csv("data/bay_link_all.csv", header = T, stringsAsFactors = F)
baylink.df$from <- tolower(baylink.df$from)
baylink.df$from <- gsub("-", "_", baylink.df$from)
baylink.df$to <- tolower(baylink.df$to)
baylink.df$to <- gsub("-", "_", baylink.df$to)
baylink.df$pair <- paste0(baylink.df$from, "_", baylink.df$to)
baylink.pairs <- baylink.df$pair


## output data and ward objects 
dataRDS <- list('raw'=mydata, 'wards'=wards,
              'data'=data, 'config'=config )
saveRDS(dataRDS, paste0("output/", output.name, "_data.rds"))

## make and output summary file
rds.files <- sapply(seq(1, number.import.chains), 
                    function(i){
                      paste0("output/", output.name, "_run_", i, ".rds")
                    })

res.list <- lapply(rds.files, function(i){
  fun.df <- readRDS(i)
  fun.df$run <- match(i, rds.files)
  return(fun.df)
}) 

res.all <- res.list[[1]]
for (i in seq(2, length(res.list))){
  res.all <- rbind(res.all, as.data.frame(res.list[[i]]))
} ; rm(i)

burnin <- 100
summary.df <- lapply(seq(1, length(res.list)), function(i){
  loop.res <- res.all[res.all$run==i,]
  loop.steps <- get_steps(loop.res, data)
  loop.steps <- add_meta(loop.steps, mydata, data, scaled_date_onset, origin, scale)
  loop.steps <- loop.steps %>% filter(step>burnin)
  loop.steps$run <- output.name
  saveRDS(loop.steps, paste0("output/", output.name, "_run", i, "_steps.rds"))
  steps <- unique(loop.steps$step)
  step.summary <- analyse_run(loop.steps[loop.steps$step==steps[1],], mydata, baylink.pairs, "wave1")
  for (step in steps[2:length(steps)]){
    summary.row <- analyse_run(loop.steps[loop.steps$step==step,], mydata, baylink.pairs, "wave1")   
    step.summary <- rbind(step.summary, summary.row)
  }
  return(step.summary)
})


summary.df <- summary.df %>% bind_rows()

saveRDS(summary.df, paste0("output/", output.name, "_summary.rds"))




