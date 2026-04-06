# Simulate data from a 2-component bivariate normal mixture
# and make a labeled scatterplot of reaction time vs tremor amplitude

library(MASS)
library(ggplot2)

set.seed(4)

# sample size
n <- 31

# mixture weight: probability of Parkinson's disease
pi_pd <- 0.45

# latent class indicator: 0 = no Parkinson's, 1 = Parkinson's
z <- rbinom(n, size = 1, prob = pi_pd)

# component-specific means
# y1 = reaction time (ms)
# y2 = tremor amplitude (mm)
mu0 <- c(320, 1.8)   # no Parkinson's
mu1 <- c(470, 5.2)   # Parkinson's

# component-specific covariance matrices
Sigma0 <- matrix(c(35^2,   12,
                   12,   0.7^2), nrow = 2, byrow = TRUE)

Sigma1 <- matrix(c(55^2,   30,
                   30,   1.2^2), nrow = 2, byrow = TRUE)

# generate observations
Y <- matrix(NA, nrow = n, ncol = 2)

for (i in 1:n) {
  if (z[i] == 0) {
    Y[i, ] <- mvrnorm(1, mu = mu0, Sigma = Sigma0)
  } else {
    Y[i, ] <- mvrnorm(1, mu = mu1, Sigma = Sigma1)
  }
}

# put into data frame
dat <- data.frame(
  reaction_time = Y[, 1],
  tremor_amplitude = Y[, 2],
  status = factor(z, levels = c(0, 1),
                  labels = c("No Parkinson's", "Parkinson's"))
)

# optional: keep tremor amplitude positive if random noise produces tiny negatives
dat$tremor_amplitude <- pmax(dat$tremor_amplitude, 0.05)

# scatterplot
p<- ggplot(dat, aes(x = reaction_time, y = tremor_amplitude)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(
    x = "Reaction time (ms)",
    y = "Tremor amplitude (mm)",
    color = "Disease status",
    title = "Motor measurements"
  ) +
  theme_minimal(base_size = 14)

ggsave(
  "~/Downloads/pd_plot.pdf",
  plot = p,
  width = 6,
  height = 6   # make square for PDF
)