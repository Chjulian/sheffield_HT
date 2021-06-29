## This script calculates the true distribution of secondary cases by using the
## observed number of secondary cases (X) and adjusting for the unobserved secondary cases (Y)
## These unobserved secondary cases are defined as those secondary cases who neither they nor their descendants
### are sampled



secondary_case_distribution <- function(wave_number, parallel = TRUE){
    
    # timer
    tictoc::tic() 
  
    # start multithread processing
  if(parallel==TRUE){
    future::plan("future::multicore")
    # plan(sequential)
  }
    
    # read in data files
    number_networks <- 100
    wavedata = list()
    if (wave_number == 1){
      wavedata <- readRDS("wave1-sec-outcomes.RDS")
    }else{wavedata <- readRDS("wave2-sec-outcomes.RDS")}
    
    cnames <- c("n.0", "n.1", "n.2", "n.3", "n.4", "n.5", "n.6", "n.7", "n.8", "n.9", "n.10")
    
    
    

    # loop through waves and calculate true distribution of secondary cases
    
    # for (wave in c(1,2)){
      # print(paste("Calculating wave ", wave))
      # grab the wave data
    
    # grab relevant dataset
    dataset <- wavedata[1:number_networks,]
    print(sprintf("analysing %i networks", number_networks))
    
    # convert every row into a list and resave
    dataset <- split(dataset, seq(nrow(dataset)))   
    
    start_time <- Sys.time()
    
          # loop through each iteration of the posterior distribution
     
       # for (iteration in 1:5){    
                  # set wave data  
        # print(iteration)          
      find_adjusted_dist <- function(data){
        
        # function to sample from an empirical distribution
        sample_function <- function(x, n, p){
          
          sample(x = x, n, replace = T, prob = p) 
          
        }
        
        # number of samples to draw from empirical distributions
        num_samples <- 5000
        
        # vector of secondary cases
        secondary_cases <- seq(0,length(dataset[[1]])-2,1)
        
        # unknown range of values of probability of true secondary cases = 0
        # prob_true0 <- seq(0.1, 1, 0.01)
        # prob_true1 <- seq(0.1, 1, 0.01)
        
        
        
       
        # prob_observed <- dataset[iteration,-1] / sum(dataset[iteration,-1], na.rm = TRUE)
        #           prob_observed[is.na(prob_observed)] <- 0
         
         prob_observed <- data[-1] / sum(data[-1], na.rm = TRUE)
         prob_observed[is.na(prob_observed)] <- 0
                  
                  # if(sum(prob_observed)!=1){
                  #   warning("prob_observed doesn't add to one")
                  # }
                  
         
                  # sample from each of the observed and unobserved secondary cases distributions
                  sample_observed <- sample_function(x = secondary_cases, n = num_samples, p = prob_observed)
                  
                  
                  
                  
                  # breaks <- c(0,0.9,1.9,2.9,3.9,4.9,5.9,6.9,7.9,8.9,9.9,10.9)
                  # calculate - for each 0 case probability - what the 0 case probability of model is
                  
                  # 1. ####### if using a mapping for brute force
                  # update value range based on previous calculations
                                  # prob_true0 <- seq(0.45, 0.7, 0.005)
                                  # prob_true1 <- seq(0.2, 0.35, 0.005)
                                  # 
                                  # # X = observed secondary cases, P(X=i)
                                  # reporting_rate <- 0.5
                                  # 
                                  # # Y = unobserved secondary cases, P(Y = i)
                                  # prob_unobserved <- function(p0, p1, i){
                                  #   (p0*(1-reporting_rate) + p1*(1-reporting_rate)^2)^i
                                  # }
                                  # 
                                  # sample_unobserved <- data.frame(purrr::map2(
                                  #   .x = rep(prob_true0,length(prob_true1)),
                                  #   .y = sort(rep(prob_true1,length(prob_true0))),
                                  #   .f = ~sample_function(x = secondary_cases, n = num_samples, p = prob_unobserved(.x, .y, secondary_cases))))
                                  # 
                                  # # calculate the distribution of true secondary cases
                                  # sample_true <- purrr::map(
                                  #   .x = sample_unobserved,
                                  #   .f = ~(.x + sample_observed))
                                  #                 # zerodensitymodel <- unlist(purrr::map(
                                  #   .x = sample_true,
                                  #   .f = ~hist(.x, breaks = breaks, plot = FALSE)$density[1] # extract density of 0 true cases 
                                  # ))
                                  # # define the break points for the density bins
                                  # breaks <- c(0,seq(0.9,floor(max(unlist(sample_true))), 1), ceiling(max(unlist(sample_true))))
                                  # onedensitymodel <- unlist(purrr::map(
                                  #   .x = sample_true,
                                  #   .f = ~hist(.x, breaks = breaks, plot = FALSE)$density[2] # extract density of 1 true cases 
                                  # ))
                                  # 
                                  # 
                                  # # calculate the difference between the model and the expectation
                                  # difference <- data.frame(
                                  #   p0 = rep(prob_true0,length(prob_true1)),
                                  #   p1 = sort(rep(prob_true1,length(prob_true0))),
                                  #   diff0 = abs(zerodensitymodel - rep(prob_true0,length(prob_true1))),
                                  #   diff1 = abs(onedensitymodel - sort(rep(prob_true1,length(prob_true0)))), row.names = NULL)
                                  # 
                                  # # absolute sum of differences
                                  # difference <- cbind(difference, diff_total = difference$diff0 + difference$diff1)
                                  # 
                                  # # minimise the sum of absolute differences
                                  # index_calibrated <- which.min(difference$diff_total)
                                  # # from this we need the rest of the distribution, which is given by 
                                  # sample_true_bestfit <- sample_true[[index_calibrated]]
                                  # # calculate the full distribution 
                                  # wave_adjustedsec_outcomes <- hist(sample_true_bestfit, breaks = breaks, plot = FALSE)$density[1:length(secondary_cases)]
                                  # # wave_adjustedsec_outcomes[iteration,] <- distribution_true_bestfit[1:length(secondary_cases)]
                                  # 
                  
                  # 2. ####### using optimisation
                  calculate_sqsum <- function(params, sample_observed, secondary_cases, num_samples){
                    
                    p0 <- params[1]
                    p1 <- params[2]
                    
                    reporting_rate <- 0.5
                    
                    
                    
                    
                    sample_unobserved <- function(p0, p1, secondary_cases, num_samples, reporting_rate){
                      
                      prob_unobserved <- function(p0, p1, i, reporting_rate){
                        (p0*(1-reporting_rate) + p1*(1-reporting_rate)^2)^i
                      }
                      
                      prob1up <- prob_unobserved(p0, p1, secondary_cases[-1])
                      probs <- c(1-sum(prob1up), prob1up)
                      su <- sample_function(x = secondary_cases, n = num_samples, p = probs)
                    }   
                    
                    sample_true <- sample_observed + sample_unobserved(p0, p1, secondary_cases, num_samples, reporting_rate)
                    
                    # define the break points for the density bins
                    breaks <- c(0,seq(0.9,floor(max(unlist(sample_true))), 1), ceiling(max(unlist(sample_true))))
                    
                    density_zero <- hist(sample_true, breaks = breaks, plot = FALSE)$density[1] # extract density of 0 true cases
                    density_one <- hist(sample_true, breaks = breaks, plot = FALSE)$density[2] # extract density of 1 true cases
                    
                    diff0 <- (density_zero - p0)^2
                    diff1 <- (density_one - p1)^2
                    
                    return(diff0+diff1)
                    
                    }
                  
                  param_out <- optim(c(0.6, 0.2), calculate_sqsum, 
                          sample_observed = sample_observed, 
                          secondary_cases = secondary_cases, 
                          num_samples = num_samples)$par
                  
                  
                  # calculate the best fit distribution
                  
                  sample_unobserved <- function(p0, p1, secondary_cases, num_samples, reporting_rate){
                    
                    prob_unobserved <- function(p0, p1, i, reporting_rate){
                      (p0*(1-reporting_rate) + p1*(1-reporting_rate)^2)^i
                    }
                    
                    
                    prob1up <- prob_unobserved(p0, p1, secondary_cases[-1])
                    probs <- c(1-sum(prob1up), prob1up)
                    su <- sample_function(x = secondary_cases, n = num_samples, p = probs)
                    }   
                  
                  sample_true <- sample_observed + sample_unobserved(param_out[1], param_out[2], secondary_cases, num_samples, reporting_rate)
                
                  breaks <- c(0,seq(0.9,floor(max(unlist(sample_true))), 1), ceiling(max(unlist(sample_true))))
                  
                  best_fit_distribution <- hist(sample_true, breaks = breaks, plot = FALSE)$density 
                  
                  return(best_fit_distribution)
                  # ))
                  # 
                  # 
                  # # calculate the difference between the model and the expectation
                  # difference <- data.frame(
                  #   p0 = rep(prob_true0,length(prob_true1)),
                  #   p1 = sort(rep(prob_true1,length(prob_true0))),
                  #   diff0 = abs(zerodensitymodel - rep(prob_true0,length(prob_true1))),
                  #   diff1 = abs(onedensitymodel - sort(rep(prob_true1,length(prob_true0)))), row.names = NULL)
                  # 
                  # # absolute sum of differences
                  # difference <- cbind(difference, diff_total = difference$diff0 + difference$diff1)
                  # 
                  # # minimise the sum of absolute differences
                  # index_calibrated <- which.min(difference$diff_total)
                  # 
                  
                  
                  
                  
                  
                  # if(wave==1){
                  #   wave1_adjustedsec_outcomes[iteration,] <- distribution_true_bestfit[1:length(secondary_cases)]
                  # # }else{
                  #   wave2_adjustedsec_outcomes[iteration,] <- distribution_true_bestfit[1:length(secondary_cases)]
                  # }
                  # print results
                  # print(paste("P(Z=0)=", difference[index_calibrated, "p0"]))
                  # print(paste("P(Z=1)=", difference[index_calibrated, "p1"]))
                  # 
                  # print(paste("Distribution of True cases for Wave", wave))
                  # print(distribution_true_bestfit)
      }
       # }
      
      if(parallel==TRUE){
      wave_adjustedsec_outcomes <- furrr::future_map(
        .x = dataset,
        .f = find_adjusted_dist,
        .options = furrr::future_options(seed = TRUE)
      )
      }else{
        wave_adjustedsec_outcomes <- purrr::map(
          .x = dataset,
          .f = find_adjusted_dist
        )
      }
      wave_adjustedsec_outcomes <- data.frame(t(sapply(wave_adjustedsec_outcomes,c)))
    
      
      # end_time <- Sys.time()
      # print(end_time - start_time)
      tictoc::toc()
      
      if(wave_number==1){
        print("saving wave 1")
        saveRDS(wave_adjustedsec_outcomes, "wave1-adjustedsec-outcomes.rds")
      }else{
        print("saving wave 2")
        saveRDS(wave_adjustedsec_outcomes, "wave2-adjustedsec-outcomes.rds")
      }
      
    # }
      future::plan("future::sequential")
}



