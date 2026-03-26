set.seed(123)

# -----------------------------
# Simulate data
# -----------------------------
m <- 30
sigma <- 1
tau <- 0.1          # small tau -> strong dependence in centered version
mu_true <- 0

theta_true <- rnorm(m, mean = mu_true, sd = tau)
y <- rnorm(m, mean = theta_true, sd = sigma)

# -----------------------------
# Gibbs sampler: centered
# -----------------------------
gibbs_centered <- function(y, sigma, tau, n_iter = 5000) {
  m <- length(y)
  
  theta <- rep(0, m)
  mu <- 0
  
  mu_save <- numeric(n_iter)
  theta1_save <- numeric(n_iter)  # save one representative theta
  
  var_theta <- 1 / (1/sigma^2 + 1/tau^2)
  sd_theta <- sqrt(var_theta)
  
  for (s in 1:n_iter) {
    # update theta_i | mu, y_i
    for (i in 1:m) {
      mean_theta <- var_theta * (y[i]/sigma^2 + mu/tau^2)
      theta[i] <- rnorm(1, mean_theta, sd_theta)
    }
    
    # update mu | theta
    mu_mean <- mean(theta)
    mu_sd <- tau / sqrt(m)
    mu <- rnorm(1, mu_mean, mu_sd)
    
    mu_save[s] <- mu
    theta1_save[s] <- theta[1]
  }
  
  data.frame(iter = 1:n_iter, mu = mu_save, theta1 = theta1_save)
}

# -----------------------------
# Gibbs sampler: non-centered
# theta_i = mu + tau * eta_i
# -----------------------------
gibbs_noncentered <- function(y, sigma, tau, n_iter = 5000) {
  m <- length(y)
  
  eta <- rep(0, m)
  mu <- 0
  
  mu_save <- numeric(n_iter)
  theta1_save <- numeric(n_iter)
  
  var_eta <- sigma^2 / (sigma^2 + tau^2)
  sd_eta <- sqrt(var_eta)
  
  for (s in 1:n_iter) {
    # update eta_i | mu, y_i
    for (i in 1:m) {
      mean_eta <- tau * (y[i] - mu) / (sigma^2 + tau^2)
      eta[i] <- rnorm(1, mean_eta, sd_eta)
    }
    
    # update mu | eta, y
    mu_mean <- mean(y - tau * eta)
    mu_sd <- sigma / sqrt(m)
    mu <- rnorm(1, mu_mean, mu_sd)
    
    mu_save[s] <- mu
    theta1_save[s] <- mu + tau * eta[1]
  }
  
  data.frame(iter = 1:n_iter, mu = mu_save, theta1 = theta1_save)
}

# -----------------------------
# Run both samplers
# -----------------------------
out_c <- gibbs_centered(y, sigma, tau, n_iter = 5000)
out_nc <- gibbs_noncentered(y, sigma, tau, n_iter = 5000)

# -----------------------------
# Compare traceplots
# -----------------------------
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

plot(out_c$iter, out_c$mu, type = "l",
     main = "Centered: trace of mu", xlab = "Iteration", ylab = expression(mu))

plot(out_nc$iter, out_nc$mu, type = "l",
     main = "Non-centered: trace of mu", xlab = "Iteration", ylab = expression(mu))

acf(out_c$mu, main = "Centered: ACF of mu")
acf(out_nc$mu, main = "Non-centered: ACF of mu")

# -----------------------------
# Compare posterior dependence
# -----------------------------
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

plot(out_c$mu, out_c$theta1, pch = 16, cex = 0.4,
     xlab = expression(mu), ylab = expression(theta[1]),
     main = "Centered: posterior draws")

plot(out_nc$mu, out_nc$theta1, pch = 16, cex = 0.4,
     xlab = expression(mu), ylab = expression(theta[1]),
     main = "Non-centered: posterior draws")

# numerical correlations
cor_centered <- cor(out_c$mu, out_c$theta1)
cor_noncentered <- cor(out_nc$mu, out_nc$theta1)

cat("Posterior correlation in centered sampler:", round(cor_centered, 3), "\n")
cat("Posterior correlation in non-centered sampler:", round(cor_noncentered, 3), "\n")