
library(actuar) ## distributions


# extra funs --------------------------------------------------------------

approx_log_fac <- function(n){ 
  a <- n*log(n) - n + log(pi*(1+4*pi*(1+2*n)))/6 + log(pi)/2
  a
}

gammabrob <- function(x) as.brob(gamma(x))


# Li-Chun Zhang models ----------------------------------------------------

## overall model

lichun_ll <- function(param, data) {
  alpha <- param[1]
  beta <- param[2]
  phi <- param[3]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  mu <- N^alpha*(n/N)^beta
  res <- m * log(mu) - (m + phi) * log(mu + phi) + (m + phi - 0.5) * log(m + phi) + 0.5 * log(phi)
  res
}

### model with one variable
#### sex
lichun_ll_sex <- function(param, data) {
  alpha_0 <- param[1]
  alpha_1 <- param[2]
  beta_0 <- param[3]
  phi <- param[4]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  X <- as.matrix(data$sex)
  mu <- N^(alpha_0+alpha_1*X)*(n/N)^(beta_0)
  res <- m * log(mu) - (m + phi) * log(mu + phi) + (m + phi - 0.5) * log(m + phi) + 0.5 * log(phi)
  res
}

#### ukr
lichun_ll_ukr <- function(param, data) {
  alpha_0 <- param[1]
  alpha_1 <- param[2]
  beta_0 <- param[3]
  phi <- param[4]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  X <- as.matrix(data$ukr)
  mu <- N^(alpha_0+alpha_1*X)*(n/N)^(beta_0)
  res <- m * log(mu) - (m + phi) * log(mu + phi) + (m + phi - 0.5) * log(m + phi) + 0.5 * log(phi)
  res
}

### model with two variable

lichun_ll_two <- function(param, data) {
  alpha_0 <- param[1]
  alpha_1 <- param[2]
  alpha_2 <- param[3]
  beta_0 <- param[4]
  phi <- param[5]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  X1 <- as.matrix(data$sex)
  X2 <- as.matrix(data$ukr)
  mu <- N^(alpha_0+alpha_1*X1+alpha_2*X2)*(n/N)^(beta_0)
  res <- m * log(mu) - (m + phi) * log(mu + phi) + (m + phi - 0.5) * log(m + phi) + 0.5 * log(phi)
  res
}



## grads
### overall

lichun_ll_grad <- function(param, data) {
  alpha <- param[1]
  beta <- param[2]
  phi <- param[3]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  mu <- N^alpha*(n/N)^beta
  x <- cbind(log(N),log(n/N))
  d_gamma <- as.numeric((m - mu) / (mu + phi)) * (phi* x )
  d_phi <- -log(mu + phi) - (m + phi)/(mu + phi) + log(m + phi) + (m + phi -0.5) / (m + phi) + 1/(2*phi)
  cbind(d_gamma, d_phi)
  #c(alpha=d_gamma[1], beta = d_gamma[2], phi=sum(d_phi))
}


## with sex
lichun_ll_sex_grad <- function(param, data) {
  alpha_0 <- param[1]
  alpha_1 <- param[2]
  beta_0 <- param[3]
  phi <- param[4]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  X <- as.matrix(data$sex)
  mu <- N^(alpha_0+alpha_1*X)*(n/N)^(beta_0)
  
  x <- cbind(log(N), log(N)*X, log(n/N))
  
  d_gamma <- as.numeric((m - mu) / (mu + phi)) * (phi* x )
  d_phi <- -log(mu + phi) - (m + phi)/(mu + phi) + log(m + phi) + (m + phi -0.5) / (m + phi) + 1/(2*phi)
  cbind(d_gamma, d_phi)
  #c(alpha=d_gamma[1], beta = d_gamma[2], phi=sum(d_phi))
}

#### with ukr

lichun_ll_ukr_grad <- function(param, data) {
  alpha_0 <- param[1]
  alpha_1 <- param[2]
  beta_0 <- param[3]
  phi <- param[4]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  X <- as.matrix(data$ukr)
  mu <- N^(alpha_0+alpha_1*X)*(n/N)^(beta_0)
  
  x <- cbind(log(N), log(N)*X, log(n/N))
  
  d_gamma <- as.numeric((m - mu) / (mu + phi)) * (phi* x )
  d_phi <- -log(mu + phi) - (m + phi)/(mu + phi) + log(m + phi) + (m + phi -0.5) / (m + phi) + 1/(2*phi)
  cbind(d_gamma, d_phi)
  #c(alpha=d_gamma[1], beta = d_gamma[2], phi=sum(d_phi))
}

#### both variables

lichun_ll_two_grad <- function(param, data) {
  alpha_0 <- param[1]
  alpha_1 <- param[2]
  alpha_2 <- param[3]
  beta_0 <- param[4]
  phi <- param[5]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  X1 <- as.matrix(data$sex)
  X2 <- as.matrix(data$ukr)
  mu <- N^(alpha_0+alpha_1*X1+alpha_2*X2)*(n/N)^(beta_0)
  
  x <- cbind(log(N), log(N)*X1, log(N)*X2, log(n/N))
  
  d_gamma <- as.numeric((m - mu) / (mu + phi)) * (phi* x )
  d_phi <- -log(mu + phi) - (m + phi)/(mu + phi) + log(m + phi) + (m + phi -0.5) / (m + phi) + 1/(2*phi)
  cbind(d_gamma, d_phi)
  #c(alpha=d_gamma[1], beta = d_gamma[2], phi=sum(d_phi))
}


## hessian

### overall

lichun_ll_hess <- function(param, data) {
  alpha <- param[1]
  beta <- param[2]
  phi <- param[3]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  mu <- N^alpha*(n/N)^beta
  x <- cbind(log(N),log(n/N))
  
  #  
  d2_gamma <-  t(x) %*% (x * -as.numeric( (m+phi)/(mu+phi)*phi ))
  d2_phi <- sum( -(2*mu+phi-m)/(mu+phi)^2 + (m+phi +0.5)/(m+phi)^2 - 1/(2*phi^2) )
  d2_gamma_phi <- colSums( -as.numeric( (mu - m) / (mu+phi)^2 * mu) * x)
  
  h <- as.matrix(Matrix::bdiag(d2_gamma, d2_phi))
  colnames(h) <- rownames(h) <- c("alpha","beta","phi")
  h[3,1:2] <- h[1:2,3] <- d2_gamma_phi
  h
}

### with extra vars 

lichun_ll_sex_hess <- function(param, data) {
  alpha_0 <- param[1]
  alpha_1 <- param[2]
  beta_0 <- param[3]
  phi <- param[4]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  X <- as.matrix(data$sex)
  mu <- N^(alpha_0+alpha_1*X)*(n/N)^(beta_0)
  
  x <- cbind(log(N), log(N)*X, log(n/N))
  
  #  
  d2_gamma <-  t(x) %*% (x * -as.numeric( (m+phi)/(mu+phi)*phi ))
  d2_phi <- sum( -(2*mu+phi-m)/(mu+phi)^2 + (m+phi +0.5)/(m+phi)^2 - 1/(2*phi^2) )
  d2_gamma_phi <- colSums( -as.numeric( (mu - m) / (mu+phi)^2 * mu) * x)
  
  h <- as.matrix(Matrix::bdiag(d2_gamma, d2_phi))
  colnames(h) <- rownames(h) <- c("alpha_0","alpha_1","beta_0", "phi")
  h[4,1:3] <- h[1:3,4] <- d2_gamma_phi
  h
}

### with ukr

lichun_ll_ukr_hess <- function(param, data) {
  alpha_0 <- param[1]
  alpha_1 <- param[2]
  beta_0 <- param[3]
  phi <- param[4]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  X <- as.matrix(data$ukr)
  mu <- N^(alpha_0+alpha_1*X)*(n/N)^(beta_0)
  
  x <- cbind(log(N), log(N)*X, log(n/N))
  
  #  
  d2_gamma <-  t(x) %*% (x * -as.numeric( (m+phi)/(mu+phi)*phi ))
  d2_phi <- sum( -(2*mu+phi-m)/(mu+phi)^2 + (m+phi +0.5)/(m+phi)^2 - 1/(2*phi^2) )
  d2_gamma_phi <- colSums( -as.numeric( (mu - m) / (mu+phi)^2 * mu) * x)
  
  h <- as.matrix(Matrix::bdiag(d2_gamma, d2_phi))
  colnames(h) <- rownames(h) <- c("alpha_0","alpha_1","beta_0", "phi")
  h[4,1:3] <- h[1:3,4] <- d2_gamma_phi
  h
}

### both

lichun_ll_two_hess <- function(param, data) {
  alpha_0 <- param[1]
  alpha_1 <- param[2]
  alpha_2 <- param[3]
  beta_0 <- param[4]
  phi <- param[5]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  X1 <- as.matrix(data$sex)
  X2 <- as.matrix(data$ukr)
  mu <- N^(alpha_0+alpha_1*X1+alpha_2*X2)*(n/N)^(beta_0)
  
  x <- cbind(log(N), log(N)*X1, log(N)*X2, log(n/N))
  
  #  
  d2_gamma <-  t(x) %*% (x * -as.numeric( (m+phi)/(mu+phi)*phi ))
  d2_phi <- sum( -(2*mu+phi-m)/(mu+phi)^2 + (m+phi +0.5)/(m+phi)^2 - 1/(2*phi^2) )
  d2_gamma_phi <- colSums( -as.numeric( (mu - m) / (mu+phi)^2 * mu) * x)
  
  h <- as.matrix(Matrix::bdiag(d2_gamma, d2_phi))
  colnames(h) <- rownames(h) <- c("alpha_0","alpha_1","alpha_2", "beta_0", "phi")
  h[5,1:4] <- h[1:4,5] <- d2_gamma_phi
  h
}

## poisson

poisson_ll <- function(param, data) {
  alpha <- param[1]
  beta <- param[2]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  mu <- N^alpha*(n/N)^beta
  res <- dpois(x = m, lambda = mu, log = T)
  #res <- -mu + m*log(mu)
  res
}

poisson_grad <- function(param, data) {
  alpha <- param[1]
  beta <- param[2]
  m <- as.matrix(data$m)
  n <- as.matrix(data$n)
  N <- as.matrix(data$N)
  n_N <- n/N
  res_alpha <- -log(N)*n_N^beta*N^alpha + m*log(N) 
  res_beta <- -log(n_N)*N^alpha*(n_N^beta) + m*log(n_N)
  res <- cbind(res_alpha, res_beta)
  res
}
