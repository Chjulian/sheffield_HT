# number of samples to draw from empirical distributions
num_samples <- 10000
  
# vector of secondary cases
secondary_cases <- seq(0,10,1)

# unknown range of values of probability of true secondary cases = 0
prob_true0 <- seq(0.1, 1, 0.01)

# X = observation distribution of secondary cases
prob_observed <- c(0.8, 0.1, 0.05, 0.03, 0.02, 0, 0, 0, 0, 0, 0)

# Y = unobserved distribution of secondary cases
prob_unobserved <- function(p, i){
  (p/2)^i
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
  .f = ~hist(.x, breaks = breaks)$density[1] # extract density of 0 true cases 
))

# look at the model calculated zero case probability (points) vs the expected value (line)
plot(zerodensitymodel, prob_true0, xlim=c(0,1), ylim=c(0,1))
par(new=TRUE)
plot(c(0,1), c(0,1), type = "line")

# calculate the difference between the model and the expectation
difference <- data.frame(
  prob_true0 = prob_true0,
  value = abs(zerodensitymodel - prob_true0))

# output the best fit model probability of the zero secondary cases
  difference[difference$value==min(difference$value), prob_true0]

# from this we need the rest of the distribution, which is given by 
sample_true_bestfit <- sample_true[, which.min(difference$value)]

# calculate the 
distribution_true_bestfit <- hist(sample_true_bestfit, breaks = breaks)$density
