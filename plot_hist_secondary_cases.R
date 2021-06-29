library(magrittr)
library(ggplot2)
rm(list=ls())
# data_type can be adjusted (adjusting for unobserved cases) or raw (only observed)
plot_hist_secondary_cases <- function(data_type = "adjusted"){


  
# read in data file for adjusted secondary cases
d = 10000
sec_case_vector <- seq(0,10,1)
wavedata = list()
if (data_type == "adjusted"){
  
        print("reading in adjusted secondary cases")
        wavedata[[1]] <- readRDS("wave1-adjustedsec-outcomes.RDS")
        wavedata[[2]] <- readRDS("wave2-adjustedsec-outcomes.RDS")
}else if(data_type == "raw"){
         print("reading in unadjusted secondary cases")
        wavedata[[1]] <- readRDS("wave1-sec-outcomes.RDS")
        wavedata[[2]] <- readRDS("wave2-sec-outcomes.RDS")
        
        wavedata[[1]] <- wavedata[[1]][-1] / rowSums(wavedata[[1]][-1], na.rm = TRUE)
        wavedata[[1]][is.na(wavedata[[1]])] <- 0
        
        wavedata[[2]] <- wavedata[[2]][-1] / rowSums(wavedata[[2]][-1], na.rm = TRUE)
        wavedata[[2]][is.na(wavedata[[2]])] <- 0
}else{
  error("no data read in")
}
todaysdate <- format(Sys.Date(), "%Y%m%d")

# function to calculate quantiles
my_quantile <- function(x, probs) {
  dplyr::tibble(x = quantile(x, probs), probs = probs)$x
}

###### RUN WAVE 1 CALCULATIONS

data <- wavedata[[1]][1:d,1:length(sec_case_vector)]
colnames(data) <- paste("n.", sec_case_vector, sep="")

R0data_wave1 <- data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(mean_of_sample = sum(sec_case_vector * dplyr::c_across())) %>%
  dplyr::mutate(R0_greater_than_one = mean_of_sample>1) %>%
  dplyr::select(mean_of_sample, R0_greater_than_one) %>%
  dplyr::group_by() %>%
  dplyr::summarise(R0_mean = mean(mean_of_sample),
                   R0_lowCI = my_quantile(mean_of_sample, 0.025),
                   R0_hiCI = my_quantile(mean_of_sample, 0.975),
                   probR0_greater_than_one = sum(R0_greater_than_one)/d)

probdata_wave1 <- data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(probover5cases = sum(dplyr::c_across(7:length(sec_case_vector)))) %>%
  dplyr::select(n.0, probover5cases) %>%
  dplyr::group_by() %>%
  dplyr::summarise(prob0cases_mean = mean(n.0),
                   prob0cases_lowCI = my_quantile(n.0, 0.025),
                   prob0cases_hiCI = my_quantile(n.0, 0.975),
                    probover5cases_mean = mean(probover5cases),
                    probover5cases_lowCI = my_quantile(probover5cases, 0.025),
                    probover5cases_hiCI = my_quantile(probover5cases, 0.975))


plotdata1 <- data %>% 
  tidyr::pivot_longer(everything()) %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(mean = mean(value),
                   lowCI = my_quantile(value, 0.025),
                   hiCI = my_quantile(value, 0.975)) %>%
  dplyr::mutate(name = strtoi(name)) %>%
  dplyr::arrange(name)



###### RUN WAVE 2 CALCULATIONS
data <- wavedata[[2]][1:d,1:length(sec_case_vector)]
colnames(data) <- paste("n.", sec_case_vector, sep="")

R0data_wave2 <- data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(mean_of_sample = sum(sec_case_vector * dplyr::c_across())) %>%
  dplyr::mutate(R0_greater_than_one = mean_of_sample>1) %>%
  dplyr::select(mean_of_sample, R0_greater_than_one) %>%
  dplyr::group_by() %>%
  dplyr::summarise(R0_mean = mean(mean_of_sample),
                   R0_lowCI = my_quantile(mean_of_sample, 0.025),
                   R0_hiCI = my_quantile(mean_of_sample, 0.975),
                   probR0_greater_than_one = sum(R0_greater_than_one)/d)


probdata_wave2 <- data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(probover5cases = sum(dplyr::c_across(7:length(sec_case_vector)))) %>%
  dplyr::select(n.0, probover5cases) %>%
  dplyr::group_by() %>%
  dplyr::summarise(prob0cases_mean = mean(n.0),
                   prob0cases_lowCI = my_quantile(n.0, 0.025),
                   prob0cases_hiCI = my_quantile(n.0, 0.975),
                   probover5cases_mean = mean(probover5cases),
                   probover5cases_lowCI = my_quantile(probover5cases, 0.025),
                   probover5cases_hiCI = my_quantile(probover5cases, 0.975))

plotdata2 <- data %>% 
  tidyr::pivot_longer(everything()) %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(mean = mean(value),
                   lowCI = my_quantile(value, 0.025),
                   hiCI = my_quantile(value, 0.975)) %>%
  dplyr::mutate(name = strtoi(name)) %>%
  dplyr::arrange(name)



# COMBINE WAVES and DRAW PLOT
all_plot <- dplyr::bind_rows(plotdata1, plotdata2, .id = "wave")
all_calc <- dplyr::bind_rows(R0data_wave1, R0data_wave2, .id = "wave")
all_prob <- dplyr::bind_rows(probdata_wave1, probdata_wave2, .id = "wave")


p <- ggplot(data = all_plot) +
  geom_col(aes(x = factor(name), y = mean, fill = wave)) +
  geom_linerange(aes(x = factor(name), ymin = lowCI, ymax = hiCI)) +
  facet_wrap(vars(wave)) + 
  theme_minimal() +
  xlab("Number of secondary cases") +
  ylab("")

# save output
ggsave(
  paste("adjusted_secondary_cases", "_", todaysdate, ".pdf", sep=""),
  plot = p,
  width = 11,
  height = 8.5,
  units = "in",
  dpi = 300)  
saveRDS(all_calc, "R0calculations.rds")

#output to screen
warning(sprintf("Only using the first %i iterations", d))

print(all_calc)
print(all_prob[,1:4])
print(all_prob[,5:7])

}


