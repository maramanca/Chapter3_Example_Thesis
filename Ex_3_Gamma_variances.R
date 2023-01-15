library(mvtnorm)
library(coda)
library(robustbase)

# Comparison variances Gamma
Yield1 = c(0.001, 0.003, 0.007, 0.020,
                              0.030, 0.040, 0.041, 0.077,
                              0.100, 0.454, 0.490, 1.020)
Yield2 = c(0.020, 0.031, 0.086, 0.130,
                              0.160, 0.160, 0.180, 0.300,
                              0.400, 0.440, 0.510, 0.720, 
                              0.950)

n1 <- length(Yield1)
n2 <- length(Yield2)

s1 <- sum(Yield1)
s2 <- sum(Yield2)

p1 <- prod(Yield1)
p2 <- prod(Yield2)


c(n1,s1,p1)
c(n2,s2,p2)

# Posterior
posterior = function(theta,n,s,p){
  out = (((exp(theta[2]))^(n*exp(theta[1])))/((gamma(exp(theta[1])))^n))*
    p^(exp(theta[1])-1)*exp(theta[1]+theta[2]-s*exp(theta[2]))
  
  return(out)
}


# log-posterior
log_posterior = function(theta,n,s,p){
  out = exp(theta[1])*n*theta[2]-n*log(gamma(exp(theta[1])))+
    (exp(theta[1])-1)*log(p)+theta[1]+theta[2]-s*exp(theta[2])
  return(out)
}

# Metropolis-Hastings
MH = function(S, B, theta0, n, s, p, tau, thin){
  
  R = S+B
  mat = matrix(NA, nrow = R, ncol = length(theta0))
  acc = rep(0, R)
  
  Tau = matrix(0, length(tau), length(tau)); diag(Tau) = tau
  
  theta_chain = theta0
  
  for (m in 1:R) {
    # Sampling from the proposal
    theta_star = rmvnorm(1, mean = theta_chain, sigma = Tau)
    # Compute r
    r = min(1, exp(log_posterior(theta = theta_star, n = n, s = s, p = p)-
                     log_posterior(theta = theta_chain, n = n, s = s, p = p)))
    # Draw U from the uniform distribution in (0,1)
    u = runif(1)
    # Accept or reject theta.s
    if(r >= u) {
      theta_chain = theta_star
      acc[m] = 1
    }
    mat[m,] = theta_chain
  }
  
  m1 <- mat[(B+1):R,]
  m2 <- mat[seq(1,nrow(m1),by=thin),]
  acc1 <- acc[(B+1):R]
  acc2 <- acc1[seq(1,length(acc1),by=thin)]
  
  # Output
  list(values = m2, acc_rate=sum(acc2)/length(acc2))
}

#### Pop 1 ####

# theta_0
maximization_1 = optim(function(x) - log_posterior(x,n1,s1,p1),
                       hessian=T, method = "L-BFGS-B", par=c(0,0))
theta0_1 = maximization_1$par
theta0_1

# tau
tau_1 = 1.5*diag(solve(maximization_1$hessian))

# Simulation
S = 150000
B = 60000

set.seed(1)
posterior_sample_1 = MH(S = S, B = B, theta0 = theta0_1, n = n1, s = s1,
                        p = p1, tau = tau_1, thin = 3)

x11()
plot(posterior_sample_1$values)

# Accuracy rate
posterior_sample_1$acc_rate

# (alpha,beta)
alpha_1  = exp(posterior_sample_1$values[,1])
beta_1   = exp(posterior_sample_1$values[,2])

# Traceplots
par(mfrow=c(1,2))
plot(alpha_1, type="l")
plot(beta_1, type="l")

# Histograms
means=c(mean(alpha_1),mean(beta_1))
means
medians=c(median(alpha_1),median(beta_1))
medians

par(mfrow=c(1,2))
{hist(alpha_1,3000, freq=F,xlab="a",main="a) Histogram of alpha_1")
  abline(v=means[1], col="red",lwd=2)
  abline(v=medians[1], col="green",lwd=2,lty=3)}
{hist(beta_1,3000, freq=F,xlab="beta_1",main="b) beta_1") 
  abline(v=means[2], col="red",lwd=2)
  abline(v=medians[2], col="green",lwd=2,lty=3)}

# ACF
par(mfrow=c(1,2))
acf(alpha_1,main="ACF for alpha_1",lag.max = 100)
acf(beta_1,main="ACF for beta_1",lag.max = 100)

# (mu,sigma2)
mu_1     = alpha_1 / beta_1
sigma2_1 = alpha_1 / beta_1^2

####
x11()
plot(cumsum(mu_1) / seq_along(mu_1),type='l',col="red")
plot(cumsum(sigma2_1) / seq_along(sigma2_1),type='l',col="blue")



#### Pop 2 ####
# theta0
maximization_2 = optim(function(x) - log_posterior(x,n2,s2,p2),
                       hessian=T, method = "L-BFGS-B", par=c(1,1))
theta0_2 = maximization_2$par
theta0_2

# tau
tau_2 = 0.8*diag(solve(maximization_2$hessian))

set.seed(1)
posterior_sample_2 = MH(S = S, B = B, theta0 = theta0_2, 
                        n = n2, s = s2,
                        p = p2, tau = tau_2, thin = 3)

plot(posterior_sample_2$values)

posterior_sample_2$acc_rate

alpha_2  = exp(posterior_sample_2$values[,1])
beta_2   = exp(posterior_sample_2$values[,2])

x11()
par(mfrow=c(1,2))
plot(alpha_2, type="l")
plot(beta_2, type="l")

x11()
means=c(mean(alpha_2),mean(beta_2))
means
medians=c(median(alpha_2),median(beta_2))
medians

par(mfrow=c(1,2))
{hist(alpha_2,3000, freq=F,xlab="alpha_2",main="a) Histogram of alpha_2")
  abline(v=means[1], col="red",lwd=2)
  abline(v=medians[1], col="green",lwd=2,lty=3)}
{hist(beta_2,3000, freq=F,xlab="beta_2",main="b) beta_2") 
  abline(v=means[2], col="red",lwd=2)
  abline(v=medians[2], col="green",lwd=2,lty=3)}

x11()
par(mfrow=c(1,2))
acf(alpha_2,main="ACF for alpha_2",lag.max = 150)
acf(beta_2,main="ACF for beta_2",lag.max = 150)

mu_2     = alpha_2 / beta_2
sigma2_2 = alpha_2 / beta_2^2

x11()
plot(cumsum(mu_2) / seq_along(mu_2),type='l',col="red")
plot(cumsum(sigma2_2) / seq_along(sigma2_2),type='l',col="blue")

N = length(sigma2_1)
I = sum(sigma2_1<sigma2_2)/N
I

delta_H = 1-2*(1-I)
delta_H


