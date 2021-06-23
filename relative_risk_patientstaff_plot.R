
wave1 <- readRDS("wave1-outcomes.rds")
wave2 <- readRDS("wave2-outcomes.rds")

n <- 10
s = 2
# number_of_staff_to_patient_wave1 <- rnorm(n, 61, s)
# number_of_staff_to_patient_wave2 <- rnorm(n, 48, s)

number_of_staff_to_patient_wave1 <- wave1$staff.inpatient.n
number_of_staff_to_patient_wave2 <- wave2$staff.inpatient.n

# number_of_staff_to_staff_wave1 <-rnorm(n, 125, s)
# number_of_staff_to_staff_wave2 <- rnorm(n, 45, s)

number_of_staff_to_staff_wave1 <-wave1$staff.staff.n
number_of_staff_to_staff_wave2 <- wave2$staff.staff.n

# number_of_patient_to_patient_wave1 <- rnorm(n, 106, s)
# number_of_patient_to_patient_wave2 <- rnorm(n, 74, s)

number_of_patient_to_patient_wave1 <- wave1$inpatient.inpatient.n
number_of_patient_to_patient_wave2 <- wave2$inpatient.inpatient.n


# number_of_patient_to_staff_wave1 <- rnorm(n, 100, s)
# number_of_patient_to_staff_wave2 <- rnorm(n, 183, s)

number_of_patient_to_staff_wave1 <- wave1$inpatient.staff.n
number_of_patient_to_staff_wave2 <- wave2$inpatient.staff.n

# calculate relative risks without the multiplier of staff to patient ratio

rel_risk_s2p_vs_s2s_wave1 <- number_of_staff_to_patient_wave1 / number_of_staff_to_staff_wave1
rel_risk_s2p_vs_s2s_wave2 <- number_of_staff_to_patient_wave2 / number_of_staff_to_staff_wave2

rel_risk_p2p_vs_p2s_wave1 <- number_of_patient_to_patient_wave1 / number_of_patient_to_staff_wave1
rel_risk_p2p_vs_p2s_wave2 <- number_of_patient_to_patient_wave2 / number_of_patient_to_staff_wave2

data <- data.frame(
  'a Staff transmission (Wave 1)' = 1/rel_risk_s2p_vs_s2s_wave1,
  'c Staff transmission (Wave 2)' = 1/rel_risk_s2p_vs_s2s_wave2,
  'b Patient transmission (Wave 1)' = 1/rel_risk_p2p_vs_p2s_wave1,
  'd Patient transmission (Wave 2)' = 1/rel_risk_p2p_vs_p2s_wave2) %>%
  pivot_longer(everything())

maxratio <- ceiling(max(data$value))
data_complement <- data.frame(
  'a Staff transmission (Wave 1)' = maxratio,
  'c Staff transmission (Wave 2)' = maxratio,
  'b Patient transmission (Wave 1)' = maxratio,
  'd Patient transmission (Wave 2)' = maxratio) %>%
  pivot_longer(everything())

l <- c("Infection from staff (Wave 1)", "Infection from a patient (Wave 1)", "Infection from staff (Wave 2)", "Infection from a patient (Wave 2)")

p <- ggplot2::ggplot() +
  # geom_rect(aes(xmin = -Inf, xmax = +Inf, ymin = 0, ymax = 1),
  #           fill = "", alpha = 0.01) +
  geom_col(data = data_complement, aes(x = name, y = value, fill = "pink")) +
  geom_bar(data = data, aes(x = name, y = value, fill = "lightblue"), fun = "median", stat = "summary") +
  geom_violin(data = data, aes(x = name, y = value)) +
  geom_abline(intercept = 1, slope = 0, linetype = "dashed", color = "darkgray") + 
  scale_x_discrete(labels = l) + 
  labs(x = "", 
       y = "Ratio of staff to patients") +
       # title = "When is transmission to patients more likely than transmission to staff?") +
  ylim(0,maxratio) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  scale_fill_identity(name = "Individual patient infection risk",
                      labels = c("Lower than for a staff member", "Higher than for a staff member"),
                      guide = "legend") 

p
