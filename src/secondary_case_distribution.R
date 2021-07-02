## This script calculates the true distribution of secondary cases by using the
## observed number of secondary cases (X) and adjusting for the unobserved secondary cases (Y)
## These unobserved secondary cases are defined as those secondary cases who neither they nor their descendants
### are sampled



secondary_case_distribution_rehash <- function(wave_number, number_networks = 100, parallel = TRUE){
  
  # timer
  tictoc::tic() 
  
  # start multithread processing
  if(parallel==TRUE){
    future::plan("future::multicore")
    # plan(sequential)
  }
  
  # read in data files
  
  wavedata = list()
  if (wave_number == 1){
    wavedata <- readRDS("./output/wave1-sec-outcomes.RDS")
  }else{wavedata <- readRDS("./output/wave2-sec-outcomes.RDS")}
  
  cnames <- c("n.0", "n.1", "n.2", "n.3", "n.4", "n.5", "n.6", "n.7", "n.8", "n.9", "n.10")
  
  
  # grab relevant dataset
  dataset <- wavedata[1:number_networks,]
  print(sprintf("analysing %i networks", number_networks))
  
  # convert every row into a list and resave
  dataset <- split(dataset, seq(nrow(dataset)))   

  find_adjusted_dist <- function(data){
    
    # function to sample from an empirical distribution
    # sample_function <- function(x, n, p){
    #   
    #   sample(x = x, n, replace = T, prob = p) 
    #   
    # }
    # 
    # number of samples to draw from empirical distributions
    # num_samples <- 5000
    
    # vector of possible true cases
    vector_true_cases <- seq(0,20,1)
    # vector over which secondary cases have been identified
    vector_observed_cases <- seq(0,10,1)
    
    # reporting rate of patients / staff
    reporting_rate <- 0.5
    
    # vector of secondary cases
    secondary_cases <- seq(0,length(dataset[[1]])-2,1)
    
    
    prob_observed <- data[-1] / sum(data[-1], na.rm = TRUE)
    prob_observed[is.na(prob_observed)] <- 0
    
    
  
    
    # sample from each of the observed and unobserved secondary cases distributions
    # sample_observed <- sample_function(x = secondary_cases, n = num_samples, p = prob_observed)
    
    
    
    

    
    # 2. ####### using optimisation
    # calculate_sqsum <- function(params, sample_observed, secondary_cases, num_samples, reporting_rate){
    calculate_sqsum <- function(params, observed_distribution, vector_true_cases, vector_observed_cases, reporting_rate){
      
      p0 <- params[1]
      p1 <- params[2]
      
      prob_unobserved <- function(p0, p1, i, reporting_rate){
            (p0*(1-reporting_rate) + p1*(1-reporting_rate)^2)^i
      }
      # code to calculate the P(Z = j) for some p0 and p1
      
      # vector over which to calculate probability of true secondary cases
      # vector_true_cases <- seq(0,20,1)
      # vector over which secondary cases have been identified
      # vector_observed_cases <- seq(0,10,1)
      
      # define reporting rate
      # reporting_rate <- 0.5
      
      # define observed probability distribution
      # observed_distribution <- c(0.9, 0.05, 0.04, 0.01, 0.00, 0.00, 0.00, 0, 0, 0, 0)
      
      #initialise conditional probability : each row is a vector of true cases
      conditional_probability <- data.frame(matrix(0, 
                                                   nrow = length(vector_true_cases), 
                                                   ncol = length(vector_observed_cases)))
      colnames(conditional_probability) <- vector_observed_cases
      rownames(conditional_probability) <- vector_true_cases
      
      # find the combinations for true vs observed cases
      for (observed_cases in vector_observed_cases){
        
        for (true_cases in vector_true_cases){
          
          if (true_cases > observed_cases){
            conditional_probability[true_cases+1,observed_cases+1] <- 
                            prob_unobserved(p0, p1, true_cases - observed_cases, reporting_rate)
              
               # (1 - reporting_rate)^(true_cases - observed_cases) # this is the probability if only direct cases form dataset
              
              
          }
          
        }
        # add on P(no unobserved cases i.e. TRUE cases = OBS cases) as complement of rest of column
        
        conditional_probability[observed_cases+1,observed_cases+1] = 1 - sum(conditional_probability[,observed_cases+1])
      }
      
      # now reweight by probability of observed cases
      # P(Z = j) = P(Z = j | X = i) * P(X = i)
      
      reweighted_true_case_distributions <- rowSums(
        sweep(conditional_probability, MARGIN=2, t(observed_distribution), `*`))
      
      true_case_distribution <- reweighted_true_case_distributions / sum(reweighted_true_case_distributions)
    
      sumsqout <- unname((true_case_distribution[1] - p0)^2 + (true_case_distribution[2] - p1)^2)
      
      return(sumsqout)
      
    }
    
    param_out <- optim(c(0, 0.5), calculate_sqsum, 
                       observed_distribution = prob_observed, 
                       vector_true_cases = vector_true_cases, 
                       vector_observed_cases = vector_observed_cases,
                       reporting_rate = reporting_rate)$par
    
  
      
    ##### calculate the best fit distribution
    
              p0 <- param_out[1]
              p1 <- param_out[2]
              
              prob_unobserved <- function(p0, p1, i, reporting_rate){
                (p0*(1-reporting_rate) + p1*(1-reporting_rate)^2)^i
              }
              # code to calculate the P(Z = j) for some p0 and p1
              
              # vector over which to calculate probability of true secondary cases
              # vector_true_cases <- seq(0,20,1)
              # vector over which secondary cases have been identified
              # vector_observed_cases <- seq(0,10,1)
              
              # define reporting rate
              # reporting_rate <- 0.5
              
              # define observed probability distribution
              # observed_distribution <- c(0.9, 0.05, 0.04, 0.01, 0.00, 0.00, 0.00, 0, 0, 0, 0)
              
              #initialise conditional probability : each row is a vector of true cases
              conditional_probability <- data.frame(matrix(0, 
                                                           nrow = length(vector_true_cases), 
                                                           ncol = length(vector_observed_cases)))
              colnames(conditional_probability) <- vector_observed_cases
              rownames(conditional_probability) <- vector_true_cases
              
              # find the combinations for true vs observed cases
              for (observed_cases in vector_observed_cases){
                
                for (true_cases in vector_true_cases){
                  
                  if (true_cases > observed_cases){
                    conditional_probability[true_cases+1,observed_cases+1] <- 
                      prob_unobserved(p0, p1, true_cases - observed_cases, reporting_rate)
                    
                    # (1 - reporting_rate)^(true_cases - observed_cases) # this is the probability if only direct cases form dataset
                    
                    
                  }
                  
                }
                # add on P(no unobserved cases i.e. TRUE cases = OBS cases) as complement of rest of column
                
                conditional_probability[observed_cases+1,observed_cases+1] = 1 - sum(conditional_probability[,observed_cases+1])
              }
              
              # now reweight by probability of observed cases
              # P(Z = j) = P(Z = j | X = i) * P(X = i)
              
              reweighted_true_case_distributions <- rowSums(
                sweep(conditional_probability, MARGIN=2, t(prob_observed), `*`))
              true_case_distribution <- reweighted_true_case_distributions / sum(reweighted_true_case_distributions)
              
              return(true_case_distribution)
              
              
          
            
  }
  
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
  
  tictoc::toc()
  
  if(wave_number==1){
    print("saving wave 1")
    saveRDS(wave_adjustedsec_outcomes, "./output/wave1-adjustedsec-outcomes.rds")
  }else{
    print("saving wave 2")
    saveRDS(wave_adjustedsec_outcomes, "./output/wave2-adjustedsec-outcomes.rds")
  }
  
  # }
  future::plan("future::sequential")
}






