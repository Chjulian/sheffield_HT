## This script calculates the true distribution of secondary cases by using the
## observed number of secondary cases (X) and adjusting for the unobserved secondary cases (Y)
## These unobserved secondary cases are defined as those secondary cases who neither they nor their descendants
### are sampled

# read in data files
wavedata = list()
wavedata[[1]] <- readRDS("wave1-sec-outcomes.RDS")
wavedata[[2]] <- readRDS("wave2-sec-outcomes.RDS")

cnames <- c("n.0", "n.1", "n.2", "n.3", "n.4", "n.5", "n.6", "n.7", "n.8", "n.9", "n.10")

# number of samples to draw from empirical distributions
num_samples <- 5000
  
# vector of secondary cases
secondary_cases <- seq(0,dim(wavedata[[1]])[2]-2,1)

# # initialise data 
# prob_observed_waves <- data.frame(matrix(0, nrow=length(secondary_cases), ncol = 2))
wave1_adjustedsec_outcomes <- data.frame(matrix(0, nrow = dim(wavedata[[1]])[1], ncol = dim(wavedata[[1]])[2]-1))
wave2_adjustedsec_outcomes <- data.frame(matrix(0, nrow = dim(wavedata[[2]])[1], ncol = dim(wavedata[[2]])[2]-1))
colnames(wave1_adjustedsec_outcomes) <- cnames
colnames(wave2_adjustedsec_outcomes) <- cnames
# unknown range of values of probability of true secondary cases = 0
prob_true0 <- seq(0.1, 1, 0.01)
prob_true1 <- seq(0.1, 1, 0.01)

# X = observed secondary cases, P(X=i)
# prob_observed_waves[,1] <- c(0.807, 0.122, 0.045, 0.017, 0.006, 0.002, 0.001, 0, 0, 0, 0)
# prob_observed_waves[,2] <- c(0.752, 0.152, 0.061, 0.022, 0.008, 0.003, 0.002, 0, 0, 0, 0)
reporting_rate <- 0.5
# par(mfrow=c(2,1))

# Y = unobserved secondary cases, P(Y = i)
prob_unobserved <- function(p0, p1, i){
(p0*(1-reporting_rate) + p1*(1-reporting_rate)^2)^i
}

# function to sample from an empirial distribution
sample_function <- function(x, n, p){
  
  sample(x = x, n, replace = T, prob = p) 
  
}

# loop through waves and calculate true distribution of secondary cases

for (wave in c(1,2)){
  print(paste("Calculating wave ", wave))
  # grab the wave data
  dataset <- wavedata[[wave]]
  
      # loop through each iteration of the posterior distribution
  #for (iteration in 1:5){    
  for (iteration in 1:dim(dataset)[1]){
              # set wave data  
              prob_observed <- dataset[iteration,-1] / sum(dataset[iteration,-1], na.rm = TRUE)
              prob_observed[is.na(prob_observed)] <- 0
              
              if(sum(prob_observed)!=1){
                warning("prob_observed doesn't add to one")
              }
              
              
              # sample from each of the observed and unobserved secondary cases distributions
              sample_observed <- sample_function(x = secondary_cases, n = num_samples, p = prob_observed)
              sample_unobserved <- data.frame(purrr::map2(
                .x = rep(prob_true0,length(prob_true0)),
                .y = sort(rep(prob_true1,length(prob_true1))),
                .f = ~sample_function(x = secondary_cases, n = num_samples, p = prob_unobserved(.x, .y, secondary_cases))))
              
              # calculate the distribution of true secondary cases
              sample_true <- purrr::map(
                .x = sample_unobserved,
                .f = ~(.x + sample_observed))
              
              # define the break points for the density bins
              breaks <- c(0,seq(0.9,floor(max(unlist(sample_true))), 1), ceiling(max(unlist(sample_true))))
              
              # calculate - for each 0 case probability - what the 0 case probability of model is
              zerodensitymodel <- unlist(purrr::map(
                .x = sample_true,
                .f = ~hist(.x, breaks = breaks, plot = FALSE)$density[1] # extract density of 0 true cases 
              ))
              
              onedensitymodel <- unlist(purrr::map(
                .x = sample_true,
                .f = ~hist(.x, breaks = breaks, plot = FALSE)$density[2] # extract density of 1 true cases 
              ))
              
              # look at the model calculated zero case probability (points) vs the expected value (line)
              # plot(zerodensitymodel, rep(prob_true0,length(prob_true0)) , xlim=c(0,1), ylim=c(0,1))
              # par(new=TRUE)
              # plot(c(0,1), c(0,1), type = "l", col = "red")
              # 
              # plot(onedensitymodel, sort(rep(prob_true1,length(prob_true1))), xlim=c(0,1), ylim=c(0,1))
              # par(new=TRUE)
              # plot(c(0,1), c(0,1), type = "l", col = "red")
              
              # calculate the difference between the model and the expectation
              difference <- data.frame(
                p0 = rep(prob_true0,length(prob_true0)),
                p1 = sort(rep(prob_true1,length(prob_true1))),
                diff0 = abs(zerodensitymodel - rep(prob_true0,length(prob_true0))),
                diff1 = abs(onedensitymodel - sort(rep(prob_true1,length(prob_true1)))), row.names = NULL)
              
              # absolute sum of differences
              difference <- cbind(difference, diff_total = difference$diff0 + difference$diff1)
              
              # minimise the sum of absolute differences
              index_calibrated <- which.min(difference$diff_total)
              
              
              # from this we need the rest of the distribution, which is given by 
              sample_true_bestfit <- sample_true[[index_calibrated]]
              
              # calculate the full distribution 
              distribution_true_bestfit <- hist(sample_true_bestfit, breaks = breaks,
                                                main = paste("wave", wave),
                                                xlab = "True secondary cases")$density
              
              if(wave==1){
                wave1_adjustedsec_outcomes[iteration,] <- distribution_true_bestfit[1:length(secondary_cases)]
              }else{
                wave2_adjustedsec_outcomes[iteration,] <- distribution_true_bestfit[1:length(secondary_cases)]
              }
              # print results
              # print(paste("P(Z=0)=", difference[index_calibrated, "p0"]))
              # print(paste("P(Z=1)=", difference[index_calibrated, "p1"]))
              # 
              # print(paste("Distribution of True cases for Wave", wave))
              # print(distribution_true_bestfit)
          }
      }
  
  saveRDS(wave1_adjustedsec_outcomes, "wave1-adjustedsec-outcomes.rds")
  saveRDS(wave2_adjustedsec_outcomes, "wave2-adjustedsec-outcomes.rds")
