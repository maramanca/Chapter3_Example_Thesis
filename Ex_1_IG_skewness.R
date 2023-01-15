#### IG skewness example ####
library(mvtnorm)

#### Data input of the two populations and setting the population parameters ####

x1 <- c(1460, 4050, 3570, 2060, 1300, 1390, 1720, 6280, 1360, 7440)

x2 <- c(5320, 1400, 3240, 2710, 4520, 4840, 8320, 13900, 71500, 6250)

sk1 <- 3*(sd(x1)/mean(x1))
sk2 <- 3*(sd(x2)/mean(x2))
  
n1 <- length(x1)
n2 <- length(x2)

t1 <- sum(x1)
t2 <- sum(x2)

m1 <- sum(1/x1)
m2 <- sum(1/x2)

# Inverse Gaussian log-posterior distribution with parameter theta = (mu, lambda)
log_posterior_IG <- function(theta, n, t, m){
  
  pd.1 <- ((n + 1)/2)*log(theta[2])
  pd.2 <- (-1/2)*log(theta[1])
  pd.3 <- (-t/2)*(theta[2]/(theta[1]^2))
  pd.4 <- n*(theta[2]/theta[1])
  pd.5 <- -1/2*m*theta[2]
  
  out = pd.1 + pd.2 + pd.3 + pd.4 + pd.5
  
  return(out)
}



#### Metropolis-Hastings algorithm for the Inverse Gaussian ####

set.seed(5)
MH_IG <- function(S, B, theta0, n, t, m, tau, thin){
  
  R <- B + S
  mat <- matrix(NA, nrow = R, ncol = length(theta0))
  acc <- rep(0, R)
  
  theta_chain <- theta0
  
  for (i in 1:R) {
    # Sampling from the proposal
    l_theta_star <- rmvnorm(1, mean = log(theta_chain), sigma = tau)
    theta_star <- exp(l_theta_star)
    
    # Computation of r
    l1 <- log_posterior_IG(theta = theta_star, n = n, t = t, m = m)
    l2 <- log_posterior_IG(theta = theta_chain, n = n, t = t, m = m)
    r <- min(1, exp(l1 - l2))
    
    # Drawing U from the uniform distribution in (0,1)
    u <- runif(1)
    
    # Accepting or rejecting theta_star
    if(r >= u) {
      theta_chain <- theta_star
      acc[i] <- 1
    }
    mat[i,] <- theta_chain
  }
  
  # Thinning
  m1 <- mat[(B+1):R,]
  m2 <- m1[seq(1, nrow(m1), by = thin),]
  acc1 <- acc[(B + 1):R]
  acc2 <- acc1[seq(1, length(acc1), by = thin)]
  
  # Output from thinning and discarding burn-in values
  list(values = m2, acc_rate = sum(acc2)/length(acc2))
}



#### Application of the sampling algorithm to Population 1 ####

# theta_0 parameter initialization
maximization_1 <- optim(function(x) -log_posterior_IG(x, n = n1, t = t1, m = m1),
                       hessian = T, method = "BFGS", par = c(2800,6500))
theta0_1 <- maximization_1$par
theta0_1


# tau covariance matrix for the proposal distribution
# tau_1 <- 0.0000015*solve(maximization_1$hessian)
# tau_1
tau_1 <- matrix(c(1,1,1,1),2,2)

# Sampling procedure
S <- 600000 # number of samples
B <- 300000  # burn-in samples

set.seed(1)
posterior_sample_1 <- MH_IG(S = S, B = B, theta0 = theta0_1, n = n1, t = t1,
                        m = m1, tau = tau_1, thin = 5)

mu_1 <- posterior_sample_1$values[, 1]
lambda_1 <- posterior_sample_1$values[, 2]

### Convergence diagnostics ###
# Samples visualization, mu1 vs lambda1
plot(posterior_sample_1$values, xlab = "mu_1", ylab = "lambda_1")

# Acceptance rate
posterior_sample_1$acc_rate

# Traceplots
par(mfrow = c(1,2))
plot(mu_1, type = "l")
plot(lambda_1, type = "l")

# Istogrammi
means <- c(mean(mu_1), mean(lambda_1))
means
medians <- c(median(mu_1), median(lambda_1))
medians

par(mfrow = c(1,2))
{hist(mu_1, 3000, freq = F, xlab = "mu_1", main = "a) Histogram of mu_1")
  abline(v = means[1], col = "red", lwd = 2)
  abline(v = medians[1], col = "green", lwd = 2, lty = 3)}
{hist(lambda_1, 3000, freq = F, xlab = "lambda_1", main = "b) lambda_1") 
  abline(v = means[2], col = "red", lwd = 2)
  abline(v = medians[2], col = "green", lwd = 2, lty = 3)}

# ACF
par(mfrow = c(1,2))
acf(mu_1, main = "ACF for mu_1", lag.max = 10000)
acf(lambda_1, main = "ACF for lambda_1", lag.max = 10000)



#### Application of the sampling algorithm to Population 2 ####

# theta_0 parameter initialization
maximization_2 <- optim(function(x) -log_posterior_IG(x, n = n2, t = t2, m = m2),
                       hessian = T, method = "L-BFGS-B", par = c(2800,6500))
theta0_2 <- maximization_2$par

# tau covariance matrix for the proposal distribution
# tau_2 <- 0.0009*solve(maximization_2$hessian)
tau_2 <- matrix(c(1,1,1,1),2,2)

# Sampling procedure
S <- 600000 # number of samples
B <- 300000  # burn-in samples

set.seed(1)
posterior_sample_2 <- MH_IG(S = S, B = B, theta0 = theta0_2, 
                        n = n2, t = t2,
                        m = m2, tau = tau_2, thin = 5)

mu_2 <- posterior_sample_2$values[, 1]
lambda_2 <- posterior_sample_2$values[, 2]

### Convergence diagnostics ###
# Samples visualization, mu_2 vs lambda_2
plot(posterior_sample_2$values, xlab = "mu_2", ylab = "lambda_2")

# Acceptance rate
posterior_sample_2$acc_rate

# Traceplots
par(mfrow = c(1,2))
plot(mu_2, type = "l")
plot(lambda_2, type = "l")

# Istogrammi
means <- c(mean(mu_2), mean(lambda_2))
means
medians <- c(median(mu_2), median(lambda_2))
medians

par(mfrow = c(1,2))
{hist(mu_2, 3000, freq = F, xlab = "mu_2",main = "a) Histogram of mu_2")
  abline(v = means[1], col = "red", lwd = 2)
  abline(v = medians[1], col = "green", lwd = 2, lty = 3)}
{hist(lambda_2, 3000, freq = F, xlab = "lambda_2", main = "b) lambda_2") 
  abline(v = means[2], col = "red", lwd = 2)
  abline(v = medians[2], col = "green", lwd = 2, lty = 3)}

# ACF
par(mfrow = c(1,2))
acf(mu_2, main = "ACF for mu_2", lag.max = 1000)
acf(lambda_2, main = "ACF for lambda_2", lag.max = 1000)



#### BDM computation and hypothesis testing H: psi_1 = psi_2 ####

N <- length(mu_1)

psi_1 <- 3*sqrt(mu_1/lambda_1) 
psi_2 <- 3*sqrt(mu_2/lambda_2)

I <- sum(psi_1 < psi_2)/N
I

d_H <- function(Int){
  
  if(Int < 0.5){
    out <- 1 - 2*I
  }
  else {
    out <- 1 - 2*(1 - Int)
  }
}

# BDM
delta_H <- d_H(Int = I)
delta_H
