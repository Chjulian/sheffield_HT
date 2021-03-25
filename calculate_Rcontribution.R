
# script to run and plot the R0 values

source("Rcontribution.R")
# make up some data

n.samples <- 1000
reporting <- runif(n.samples, 0.6, 0.7)
patient2patient <- rpois(n.samples, 2)
patient2staff <- rpois(n.samples, 1) 
staff2patient <- rpois(n.samples, 1) 
staff2staff <- rpois(n.samples, 2) 

# calculate table of R0 values
data.out <- Rcontribution(reporting = reporting,
                          patient2patient = patient2patient,
                          patient2staff = patient2staff,
                          staff2patient, staff2staff) 

data.to.plot.all <- data.out %>%
      tidyr::pivot_longer(names_to = "R0type", values_to = "Value", cols = 1:4)


# find the proportion of simulations that result in R < 1

data.to.plot.all %>% 
  group_by(R0type, .drop = FALSE) %>%
  filter(Value < 1) %>%
  summarise(prop = n() / n.samples)

# plot output
p <- ggplot2::ggplot() +
          geom_violin(data = data.to.plot.all, aes(x = R0type, y = Value)) + 
          geom_abline(slope = 0, intercept = 1, color = "red", lty = 2) +
          theme(axis.text.x = element_text(angle = 70, hjust=1),
                  axis.title = element_blank())
p

