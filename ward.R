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
               ggplot2) 
source("src/functions_wards.R")
source("src/functions_maxdist.R")
source("src/onset_distribution.R")

max_dist <- 0


#################################
#load linelist data
#################################
# load data, clean dates and filter-out observations
mydata <- read.csv("data/hoci-phylo-metadata-onset-2020-09-16.V1.2.121120.csv")
mydata <- onsetdistribution(mydata)
mydata<- clean_data(mydata, guess_dates = which(names(mydata)=='DateOfOnset_forAnalysis'))
mydata <- mydata[mydata$category%in%c('inpatient', 'outpatient','staff'),]
mydata <- mydata[mydata$loc1%!in%c('com'),] #rm community associated
mydata <- mydata[!is.na(mydata$dateofonset_foranalysis),]

#################################
#load dna data in DNAbin format, 
#################################
#one sequence per individual
dna <- read.FASTA("data/hoci_phylo-2020-09-03_error-removed_masked.aln")
mylabels <- clean_data(as.data.frame(names(dna)))
names(dna) <- mylabels$names_dna
dna <- dna[mydata$barcode]; mylabels <- names(dna)
write.FASTA(dna, "data/hoci_phylo-2020-09-03_error-removed_masked.ord.aln")
dna.SNP <- import_fasta_sparse("data/hoci_phylo-2020-09-03_error-removed_masked.ord.aln") %>% snp_dist(.)
rownames(dna.SNP) <- colnames(dna.SNP) <- mylabels; rm(mylabels)

#################################
#process onset dates
#################################
date_onset<- mydata$dateofonset_foranalysis
names(date_onset) <- rep("local", length(date_onset)) 

#################################
#Consolidate ward data
#################################
wards.patients<- read.csv("data/patient-ward-movements-20201207.csv")
wards.patients<- wards.patients[,c('Barcode','Ward','In','Out')]
names(wards.patients) <- c('id','ward','adm','dis') #rename as examples files for consistency
wards<- ward_occupation(wards.patients=wards.patients)
#set nncl and multiple as NA
wards$ward[wards$ward%in%c('nncl','multiple')] <- NA 

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
        n_iter = 5e3, sample_every = 50,
        ## prior on the reporting probability pi
        move_pi = TRUE, prior_pi = c(5, 5),
        ## eps is the probability of infections occuring between cases registered on the same ward
        move_eps = TRUE, prior_eps = c(1, 1),
        ## tau is the probability of an unobserved cases being moved to a different ward
        move_tau = TRUE, prior_tau = c(1, 1),
        ## leave these options as is for now
        move_joint = TRUE, move_model = TRUE
)


#################################
#find a starting point
#################################
## Identify initial tree (run this for now, it will identify a possible starting point)
#from time to time it fails but usually work after first try
config <- get_initial_tree(data, config, n_iter = 5e3, max_dist = max_dist)
table(data$D[cbind(seq_along(config$init_alpha), config$init_alpha)],useNA = 'always')

#################################
#run outbreaker
#################################
res <- outbreaker(data, config)

#################################
#summaries
#################################

tChains.df <- summary.res(res=res_maxdist_0, burnin=1000, support=0)

#################################
#save data
#################################

myRDS <- list('data'=data, 'config'=config , 'res'=res, 'df'=tChains.df)
saveRDS(myRDS, file = paste("output/res_", "max_dist_", max_dist, "_", format(Sys.time(), "%Y_%m_%d_%H_%M"), ".rds", sep=""))
