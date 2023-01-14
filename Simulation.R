set.seed(1)

# Posterior distribution
rnormgamma <- function(n, eta, nu, alpha, beta) {
  
  if (length(n) > 1) 
  n <- length(n)
  phi <- rgamma(n, alpha, beta)
  mu <- rnorm(n, eta, 1/sqrt(nu*phi))
  
  data.frame(mu = mu, phi = phi)
}

# Function that gives the hyperparameters of the posterior distribution
hyperparam <- function(n, mean, s_square){
  eta <- mean          
  nu <- n
  alpha <- (n - 1)/2
  beta <- (n*s_square)/2
  list(eta = eta, nu = nu, alpha = alpha, beta = beta)
} 

# Function that computes the BDM 
d_H <- function(Int){
  if(Int < 0.5){
    out <- 1 - 2*Int
  }
  else {
    out <- 1 - 2*(1 - Int)
  }
}

mu_1 <- mu_2 <-  3

std_1 <- std_2 <- 1
  
S <- 50000     # number of simulations
n_1 <- 1000    # sample size
n_2 <- 1000    # sample size
D <- 10000     # number of posterior draws


x1 <- matrix(NA, n_1, S)      # sample 1
x2 <- matrix(NA, n_2, S)      # sample 2

I <- rep(NA, S)               # Integral for the BDM
delta_H <- rep(NA, S)         # BDM

for (i in 1:S) {
  x1[,i] <- rnorm(n = n_1, mean = mu_1, sd=std_1)     
  x1_mean <- mean(x1[,i])
  s_square1 <- var(x1[,i])
  
  x2[,i] <- rnorm(n = n_2, mean = mu_2, sd=std_2)
  x2_mean <- mean(x2[,i])
  s_square2 <- var(x2[,i])
  
  
  # Hyperparameters of the two populations
  h1 <- hyperparam(n = n_1, mean = x1_mean, s_square = s_square1)
  h2 <- hyperparam(n = n_2, mean = x2_mean, s_square = s_square2)
  
  # Posterior draws 
  phi_mu_1 <- rnormgamma(n = D, 
                         eta = h1$eta, nu = h1$nu, 
                         alpha = h1$alpha, beta = h1$beta)

  phi_mu1_mu <- phi_mu_1$mu
  phi_mu1_phi <- phi_mu_1$phi
  
    phi_mu_2 <- rnormgamma(n = D, 
                         eta = h2$eta, nu = h2$nu, 
                         alpha = h2$alpha, beta = h2$beta)
  
  phi_mu2_mu <- phi_mu_2$mu
  phi_mu2_phi <- phi_mu_2$phi
  
  psi_1 <- 1/(abs(phi_mu1_mu)*(phi_mu1_phi)^(0.5))
  psi_2 <- 1/(abs(phi_mu2_mu)*(phi_mu2_phi)^(0.5))
  
  I[i] <- sum(psi_1<psi_2)/D
  
  delta_H[i] <- d_H(Int = I[i])
  
  print(i)
}

# Summary for each evidence measure
summary(delta_H)

# False positive rate for threshold = 0.90, 0.95, 0.99
  round(sum(delta_H>=0.9)/S,3)

  round(sum(delta_H>=0.95)/S,3)

  round(sum(delta_H>=0.99)/S,3)
  
  
# Comparisons cvs
  mean_1_vec <- apply(x1,2,mean)
  std1_vec <- apply(x1,2,sd)
  cv1 <- std1_vec/abs(mean_1_vec)
  
  mean_2_vec <- apply(x2,2,mean)
  std2_vec <- apply(x2,2,sd)
  cv2 <- std2_vec/abs(mean_2_vec)
  
  diff <- abs(cv1-cv2)
  summary(diff)
  hist(diff)

  
# Plot of a given chain
x11()
plot(density(theta[,5000]))

# AA AR
# RA RR
matrix(data = round(c(n_11,n_12,n_21,n_22),3),2,2, 
       byrow = TRUE)

# r prop Jeffreys
n_11 = sum(e_g<0.95 & d_H<0.95)/S
n_22 = sum(e_g>=0.95 & d_H>=0.95)/S
n_12 = sum(e_g<0.95 & d_H>=0.95)/S
n_21 = sum(e_g>=0.95 & d_H<0.95)/S

# AA AR
# RA RR
matrix(data = round(c(n_11,n_12,n_21,n_22),3),2,2, 
       byrow = TRUE)
