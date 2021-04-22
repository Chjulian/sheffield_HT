##Written by Katie
##impute symptom onset dates for those individuals without known times
onsetdistribution <- function(data="mydata") {
  if(!require(pacman)) install.packages(pacman)
  pacman::p_load(data.table, lubridate, dplyr,
                 magrittr,ggplot2,tidyr)
  data_read <- as.data.table(data)
  data <- data_read %>%
    mutate(DateOfCollection = as_date(DateOfCollection, format = "%d/%m/%Y"))%>%
    mutate(DateOfAdmission = as_date(DateOfAdmission, format = "%d/%m/%Y")) %>%
    mutate(DateOfOnset = as_date(DateOfOnset, format = "%d/%m/%Y")) %>%
    mutate(Delay_OnsetAdmission = DateOfAdmission - DateOfOnset) %>%
    mutate(Delay_OnsetSwab = DateOfCollection - DateOfOnset) %>%
    mutate(Delay_AdmissionSwab = DateOfCollection - DateOfAdmission) %>%
    mutate(DateOfOnset_forAnalysis = DateOfOnset) %>%
    mutate(Delay_OnsetSwab_forAnalysis = Delay_OnsetSwab) %>%
    mutate(HealthcareAssociation = ifelse(HealthcareAssociation == "", "undefined", HealthcareAssociation)) %>%
    mutate(HealthcareAssociation = ifelse(HealthcareAssociation == "Not admitted-Community Associated", "Not Admitted-Community Associated", HealthcareAssociation))
  
  #split up data table by type of person 
  cat_names <- unique(data$Category) 
  cat_names <- cat_names[!(cat_names=="EXTERNAL" | cat_names=="HHCON")]
  data_bycat <- purrr::map(
    .x = cat_names,
    ~dplyr::filter(data, Category == .x)) %>%
    setNames(cat_names)
  out_partial <- list()
  
  # just work with inpatients first
  loop.cats <- c("INPATIENT", "STAFF")
  if ("OUTPATIENT" %in% unique(data$Category)){
    loop.cats <- c(loop.cats, "OUTPATIENT")
  }
  
  for (pcat in loop.cats){
    
    # those data that need imputation
    to_impute <- data_bycat[[pcat]] %>%
      filter(is.na(DateOfOnset))
    
    # those data on which the imputation is calculated
    to_impute_comp <- data_bycat[[pcat]] %>%
      filter(!is.na(DateOfOnset))
    
    # check the type of infection that are imputed
    type_names_toimpute <- unique(to_impute$HealthcareAssociation)
    type_names_assigned <- unique(to_impute_comp$HealthcareAssociation)
    
    type_names_intersection <- intersect(type_names_toimpute, type_names_assigned)
    type_names_unabletobeimputated <- setdiff(type_names_toimpute, type_names_intersection)
    
    # reassign those unable to be imputed as 'undefined'
    to_impute <- to_impute %>%
      mutate(HealthcareAssociation = ifelse(HealthcareAssociation %in% type_names_unabletobeimputated, "undefined", HealthcareAssociation))
    
    # the total number of entries that need to be imputed
    number_to_impute <- to_impute %>%
      group_by(HealthcareAssociation) %>%
      summarise(n = n())
    
    # the total number of entries on which the data are imputed 
    number_assigned <- to_impute_comp %>%
      group_by(HealthcareAssociation) %>%
      summarise(n = n())
    
    # split up the data by patient type
    data_by_type_toimpute <- purrr::map(
      .x = type_names_intersection,
      ~dplyr::filter(to_impute, HealthcareAssociation == .x)) %>%
      setNames(type_names_intersection)
    
    data_by_type_toimpute_comp <- purrr::map(
      .x = type_names_intersection,
      ~dplyr::filter(to_impute_comp, HealthcareAssociation == .x)) %>%
      setNames(type_names_intersection)
    
    # data used for imputation
    data_withvals <- data_bycat[[pcat]] %>% 
      filter(Delay_OnsetSwab > -1)  
    
    data_by_type_withvals <- purrr::map(
      .x = type_names_intersection,
      ~dplyr::filter(data_withvals, HealthcareAssociation == .x)) %>%
      setNames(type_names_intersection)
    
    # sample the imputed entries from the entries that have a date of onset
    impute_vals <- purrr::map(
      .x = type_names_intersection,
      ~sample(x = data_by_type_withvals[[.x]]$Delay_OnsetSwab_forAnalysis,
              size = as.numeric(number_to_impute[number_to_impute$HealthcareAssociation==.x,"n"]),
              replace = TRUE)) %>%
      setNames(type_names_intersection)
    
    # impute_vals <- purrr::map(
    #   .x = "",
    #   ~sample(x = data_by_type_withvals[[.x]]$Delay_OnsetSwab_forAnalysis,
    #           size = as.numeric(number_to_impute[number_to_impute$HealthcareAssociation==.x,"n"]),
    #           replace = TRUE)) 
    # %>%
    #   setNames(type_names)
    # 
    # add the new values to a list with subset of original data that needed to be imputed
    imputed_data <- purrr::map(
      .x = type_names_intersection,
      ~dplyr::mutate(
        data_by_type_toimpute[[.x]], 
        Delay_OnsetSwab_forAnalysis = impute_vals[[.x]],
        DateOfOnset_forAnalysis = DateOfCollection - impute_vals[[.x]])) %>%
      setNames(type_names_intersection)
    
    # bind together the new vals and original
    out_partial_bycat <- purrr::map(
      .x = type_names_intersection,
      ~rbindlist(list(data_by_type_toimpute_comp[[.x]], imputed_data[[.x]]))
    )
    
    out_partial[[pcat]] <- rbindlist(out_partial_bycat)
  }
  
  out_data <- rbindlist(out_partial)
}



