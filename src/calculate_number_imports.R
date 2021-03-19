calculate_number_imports <- function(data=mydata){
  pacman::p_load(data.table, lubridate,
                 dplyr, magrittr, ggplot2,
                 tidyr)

  # file_name <- "dataset_imputed_2020-12-08.csv"
  # data <- data.table::fread(file_name)
  # message("Using ", file_name, " as input")
  
  # incubation period as Lognormal
  dist_setup <- function(dist_mean = NULL, dist_sd = NULL) {
    out <- purrr::partial(plnorm,
                          meanlog = dist_mean,
                          sdlog = dist_sd)
    return(out)
  }
  
  incubation_fn <- dist_setup(dist_mean = 1.621,
                              dist_sd = 0.418)
  
  # grab the inpatients only (as outpatients and staff are not subject to HOCI criteria)
  # then: calculate prob that exposure occurred before hospital admission using pre-assigned distribution
  inpatient_data <- mydata %>%
    dplyr::filter(category=="inpatient") %>%
    mutate(dateofonset_foranalysis = as_date(dateofonset_foranalysis, format = "%y-%m-%d")) %>%
    mutate(dateofadmission = as_date(dateofadmission, format = "%y-%m-%d")) %>%
    mutate(delay_admissiononset = as.numeric(dateofonset_foranalysis - dateofadmission))  %>%
    mutate(prob_communityqcquired = ifelse(delay_admissiononset > 0, 
                                           1 - incubation_fn(delay_admissiononset),
                                           1))
  
  p <- mean(inpatient_data$prob_communityqcquired) 
  N <- length(inpatient_data$prob_communityqcquired)
  number_imports_average <- round(p * N)
  number_imports_loCI <- qbinom(0.025,N,p)
  number_imports_hiCI <- qbinom(0.975,N,p)
  
  proportion_imports_average <- number_imports_average / N
  proportion_imports_loCI <- number_imports_loCI / N
  proportion_imports_hiCI <- number_imports_hiCI / N
  
  number_imports_average <<- number_imports_average
  
  # save(number_imports_average, number_imports_loCI, number_imports_hiCI,
  #      proportion_imports_average, proportion_imports_loCI, proportion_imports_hiCI,
  #      file = paste("NumberImports_", file_name,".Rda"))
}
