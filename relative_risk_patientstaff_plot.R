
n <- 10
s = 2
number_of_staff_to_patient_wave1 <- rnorm(n, 61, s)
number_of_staff_to_patient_wave2 <- rnorm(n, 48, s)

number_of_staff_to_staff_wave1 <-rnorm(n, 125, s)
number_of_staff_to_staff_wave2 <- rnorm(n, 45, s)

number_of_patient_to_patient_wave1 <- rnorm(n, 106, s)
number_of_patient_to_patient_wave2 <- rnorm(n, 74, s)

number_of_patient_to_staff_wave1 <- rnorm(n, 100, s)
number_of_patient_to_staff_wave2 <- rnorm(n, 183, s)

rel_risk_s2p_vs_s2s_wave1 <- number_of_staff_to_patient_wave1 / number_of_staff_to_staff_wave1
rel_risk_s2p_vs_s2s_wave2 <- number_of_staff_to_patient_wave2 / number_of_staff_to_staff_wave2

rel_risk_p2p_vs_p2s_wave1 <- number_of_patient_to_patient_wave1 / number_of_patient_to_staff_wave1
rel_risk_p2p_vs_p2s_wave2 <- number_of_patient_to_patient_wave2 / number_of_patient_to_staff_wave2

data <- data.frame(
  'Staff transmission (Wave 1)' = 1/rel_risk_s2p_vs_s2s_wave1,
  'Staff transmission (Wave 2)' = 1/rel_risk_s2p_vs_s2s_wave2,
  'Patient transmission (Wave 1)' = 1/rel_risk_p2p_vs_p2s_wave1,
  'Patient transmission (Wave 2)' = 1/rel_risk_p2p_vs_p2s_wave2) %>%
  pivot_longer(everything())

maxratio <- 3
data_complement <- data.frame(
  'Staff transmission (Wave 1)' = maxratio,
  'Staff transmission (Wave 2)' = maxratio,
  'Patient transmission (Wave 1)' = maxratio,
  'Patient transmission (Wave 2)' = maxratio) %>%
  pivot_longer(everything())


p <- ggplot2::ggplot() +
  # geom_rect(aes(xmin = -Inf, xmax = +Inf, ymin = 0, ymax = 1),
  #           fill = "", alpha = 0.01) +
  geom_col(data = data_complement, aes(x = name, y = value, fill = "pink")) +
  geom_bar(data = data, aes(x = name, y = value, fill = "lightblue"), fun = "mean", stat = "summary") +
  geom_violin(data = data, aes(x = name, y = value)) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 40, hjust=1)) +
  labs(x = "", 
       y = "Ratio of staff to patients",
       title = "When is transmission to patients more likely than transmission to staff?") +
  ylim(0,maxratio) + 
  scale_fill_identity(name = "",
                      labels = c("More likely", "Less likely"),
                      guide = "legend")
  

p
