# estimating R0 and contribution of different routes of transmission for HCW and inpatients

# input parameters:
# reporting: distribution samples of average reporting fractions (length n)
# patient2patient / patient2staff / staff2staff / staff2patient: distribution samples MEAN number of secondary cases of each type

# output:
# dataframe of size n x 4 giving overall R0 then R0 without pat2pat, staff2staff or staff-patient transmission


Rcontribution  <- function(reporting = NULL, 
                          patient2patient = NULL, patient2staff = NULL, 
                          staff2patient = NULL, staff2staff = NULL){
  

 
  len <- length(reporting)
  
  # calculate overall R0
  repro.values <- list(patient2patient / reporting,
                       staff2patient / reporting, 
                       patient2staff / reporting, 
                       staff2staff / reporting)
  
  R0 <- unlist(
          purrr::pmap(
                .l = repro.values,
                .f = calculate_repro_number))

  # calculate R0 without patient2patient transmission 
  repro.values_nop2p <- list(rep(0,len),
                             staff2patient / reporting, 
                             patient2staff / reporting, 
                             staff2staff / reporting)
  R0_nop2p <- unlist(
          purrr::pmap(
            .l = repro.values_nop2p,
            .f = calculate_repro_number))
  
  # calculate R0 without staff2staff transmission 
  repro.values_nos2s <- list(patient2patient / reporting,
                             staff2patient / reporting, 
                             patient2staff / reporting, 
                             rep(0,len))
  R0_nos2s <- unlist(
    purrr::pmap(
      .l = repro.values_nos2s,
      .f = calculate_repro_number))

  
  # calculate R0 without staff2patient / patient2staff transmission 
  repro.values_nos2p <- list(patient2patient / reporting,
                             rep(0,len), 
                             rep(0,len), 
                             staff2staff / reporting)
  R0_nos2p <- unlist(
    purrr::pmap(
      .l = repro.values_nos2p,
      .f = calculate_repro_number))

  
  repro.out <- dplyr::as_tibble(cbind("R0"=R0, 
                                   "R0 No patient-patient"=R0_nop2p, 
                                   "R0 No patient-staff"=R0_nos2s, 
                                   "R0 No staff-staff"=R0_nos2p))
  
  
  return(repro.out)

}


calculate_repro_number <- function(p2p, s2p, p2s, s2s){
  
  ngm <- matrix(c(p2p, s2p, p2s, s2s), nrow = 2, ncol = 2)
  return(max(eigen(ngm, only.values = TRUE)$values))
  
            
}
