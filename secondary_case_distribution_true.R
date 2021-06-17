# number of samples to draw from empirical distributions
num_samples <- 10000
  
# vector of secondary cases
secondary_cases <- seq(0,10,1)

# initialise data 
prob_observed_waves <- data.frame(matrix(0, nrow=length(secondary_cases), ncol = 2))

# unknown range of values of probability of true secondary cases = 0
prob_true0 <- seq(0.1, 1, 0.01)

# X = observation distribution of secondary cases
prob_observed_waves[,1] <- c(0.807, 0.122, 0.045, 0.017, 0.006, 0.002, 0.001, 0, 0, 0, 0)
prob_observed_waves[,2] <- c(0.752, 0.152, 0.061, 0.022, 0.008, 0.003, 0.002, 0, 0, 0, 0)
reporting_rate <- 0.5
par(mfrow=c(2,1))

# loop through waves and calculate true distribution of secondary cases

for (wave in c(1,2)){
      
      # set wave data  
      prob_observed <- prob_observed_waves[,wave]
      
      if(sum(prob_observed)!=1){
        warning("prob_observed doesn't add to one")
      }
      # Y = unobserved distribution of secondary cases
      prob_unobserved <- function(p, i){
        (reporting_rate*p)^i
      }
      
      # sample from an empirial distribution
      sample_function <- function(x, n, p){
        
        sample(x = x, n, replace = T, prob = p) 
        
      }
      
      # sample from each of the observed and unobserved secondary cases distributions
      sample_observed <- data.frame(matrix(
        rep(sample_function(x = secondary_cases, n = num_samples, p = prob_observed), length(prob_true0)),
        ncol = length(prob_true0)))
      sample_unobserved <- data.frame(purrr::map(
        .x = prob_true0,
        .f = ~sample_function(x = secondary_cases, n = num_samples, p = prob_unobserved(.x, secondary_cases))))
      
      # calculate the distribution of true secondary cases
      sample_true <- sample_observed + sample_unobserved
      
      # define the break points for the density bins
      breaks <- c(0,seq(0.9,floor(max(sample_true)), 1), ceiling(max(sample_true)))
      
      # calculate - for each 0 case probability - what the 0 case probability of model is
      zerodensitymodel <- unlist(purrr::map(
        .x = sample_true,
        .f = ~hist(.x, breaks = breaks, plot = FALSE)$density[1] # extract density of 0 true cases 
      ))
      
      # look at the model calculated zero case probability (points) vs the expected value (line)
      # plot(zerodensitymodel, prob_true0, xlim=c(0,1), ylim=c(0,1))
      # par(new=TRUE)
      # plot(c(0,1), c(0,1), type = "l")
      
      # calculate the difference between the model and the expectation
      difference <- data.frame(
        prob_true0 = prob_true0,
        value = abs(zerodensitymodel - prob_true0))
      
      # output the best fit model probability of the zero secondary cases
        difference[difference$value==min(difference$value), prob_true0]
      
      # from this we need the rest of the distribution, which is given by 
      sample_true_bestfit <- sample_true[, which.min(difference$value)]
      
      # calculate the full distribution 
      distribution_true_bestfit <- hist(sample_true_bestfit, breaks = breaks, 
                                        main = paste("wave", wave),
                                        xlab = "True secondary cases")$density
      
      print(paste("Distribution of True cases for Wave", wave))
      print(distribution_true_bestfit)
}
