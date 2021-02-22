#   /.)
#  /)\|   
# // /     
# /'" " 
#################################
#Outbreaker2 for ward studies
#################################

#################################
#Load libraries and functions 
#################################
if(!require(pacman)) install.packages(pacman)
pacman::p_load(ape, linelist, epitrix, fitdistrplus,
               outbreaker2, tidyverse, rio, 
               magrittr, remotes, epicontacts, 
               visNetwork, lubridate, pairsnp,
               ggplot2, chron,glue) 
source("src/functions_wards.R")
source("src/functions.R")
source("src/onset_distribution.R")

max_dist <- 2

#################################
#load linelist data
#################################
# load data, clean dates and filter-out observations
mydata <- read.csv("data/hoci-phylo-metadata-onset-2021-02-07.csv",
                   stringsAsFactors =F)
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
write.FASTA(dna, "data/hoci_phylo-2021-02-07_error-removed_masked.ord.aln")
dna.SNP <- import_fasta_sparse("data/hoci_phylo-2021-02-07_error-removed_masked.ord.aln") %>% snp_dist(.)
rownames(dna.SNP) <- colnames(dna.SNP) <- mylabels; rm(mylabels)

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
        dates = date_onset,
        dna = dna,
        ctd_timed = wards,
        ids = meta$id,
        w_dens = generation_time$d(1:50),
        f_dens = f,
        p_trans = list(transition)
)
rm(transition)

## set up outbreaker config file
config <- create_config(
        ## define iterations and sampling (this is a short exploratory run)
        n_iter = 1e4, sample_every = 50,
        ## prior on the rate
        move_mu = TRUE, init_mu = 0.1, sd_mu = 0.01,
        ## prior on the reporting probability pi
        move_pi = TRUE, prior_pi = c(30, 10), init_pi=0.75,
        ## eps is the probability of infections occuring between cases registered on the same ward
        move_eps = TRUE, prior_eps = c(9800.01, 98.99), init_eps=0.999,
        ## tau is the probability of an unobserved cases being moved to a different ward
        move_tau = TRUE, prior_tau = c(1, 1), init_tau=0.1,
        ## leave these options as is for now
        move_joint = TRUE, move_model = TRUE
)

# get beta prior values from mean and sd
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

#################################
#find a starting point
#################################
## Identify initial tree (run this for now, it will identify a possible starting point)
#from time to time it fails but usually work after first try

config <- get_initial_tree(
  data = data,
  config = config,
  priors = priors,
  likelihoods = likelihoods,
  n_iter = 1e3,
  max_dist = 2
)
# table(data$D[cbind(seq_along(config$init_alpha), config$init_alpha)],useNA = 'always')

#################################
#run outbreaker
#################################
res <- outbreaker(data, config, prior = priors, likelihoods = likelihoods)

#################################
#summaries
#################################

tChains.df <- summary.res(res=res, burnin=1000, support=0)

#################################
#save data
#################################

myRDS <- list('raw'=mydata, 'wards'=wards,
              'data'=data, 'config'=config , 'res'=res, 
              'res.sum'= mydf,
              'cons_tree'=cons_tree,
              'clusters'=myclusters,
              'df'=tChains.df)
saveRDS(myRDS, file = paste("output/res_", "max_dist_", max_dist, "_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".rds", sep=""))
