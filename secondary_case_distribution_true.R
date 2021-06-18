
# number of samples to draw from empirical distributions
num_samples <- 10000
  
# vector of secondary cases
secondary_cases <- seq(0,10,1)

# initialise data 
prob_observed_waves <- data.frame(matrix(0, nrow=length(secondary_cases), ncol = 2))

# unknown range of values of probability of true secondary cases = 0
prob_true0 <- seq(0.1, 1, 0.01)
prob_true1 <- seq(0.1, 1, 0.01)

# X = observed secondary cases, P(X=i)
prob_observed_waves[,1] <- c(0.807, 0.122, 0.045, 0.017, 0.006, 0.002, 0.001, 0, 0, 0, 0)
prob_observed_waves[,2] <- c(0.752, 0.152, 0.061, 0.022, 0.008, 0.003, 0.002, 0, 0, 0, 0)
reporting_rate <- 0.5
par(mfrow=c(2,1))

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
      
      # set wave data  
      prob_observed <- prob_observed_waves[,wave]
      
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

      # print results
      print(paste("P(Z=0)=", difference[index_calibrated, "p0"]))
      print(paste("P(Z=1)=", difference[index_calibrated, "p1"]))
      
      print(paste("Distribution of True cases for Wave", wave))
      print(distribution_true_bestfit)
}
