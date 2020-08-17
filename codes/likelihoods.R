tpoisson_ll <- function(p) {
  mu <- pop^(p[1])*(pol_pop)^(p[2])
  #res <- actuar::dztpois(x  = border, lambda = mu, log = TRUE)
  res <- extraDistr::dtpois(x = border, lambda = mu, log = T, a = 0)
  -sum(res)
}

tpoisson_ll_ukr <- function(p) {
  mu <- pop^(p[1]+p[2]*ukr)*(pol_pop)^(p[3])
  #res <- actuar::dztpois(x  = border, lambda = mu, log = TRUE)
  res <- extraDistr::dtpois(x = border, lambda = mu, log = T, a = 0)
  -sum(res)
}

tpoisson_ll_sex <- function(p) {
  mu <- pop^(p[1]+p[2]*sex)*(pol_pop)^(p[3])
  #res <- actuar::dztpois(x  = border, lambda = mu, log = TRUE)
  res <- extraDistr::dtpois(x = border, lambda = mu, log = T, a = 0)
  -sum(res)
}

tpoisson_ll_both <- function(p) {
  mu <- pop^(p[1]+p[2]*ukr+p[3]*sex)*(pol_pop)^(p[4])
  #res <- actuar::dztpois(x  = border, lambda = mu, log = TRUE)
  res <- extraDistr::dtpois(x = border, lambda = mu, log = T, a = 0)
  -sum(res)
}


tnegbin_ll <- function(p) {
  mu <- pop^(p[1])*(pol_pop)^(p[2])
  phi <- p[3]
  res <- actuar::dztnbinom(x  = border, size = phi, prob = phi/(phi+mu), log = TRUE)
  -sum(res)
}

tnegbin_ll_ukr <- function(p) {
  mu <- pop^(p[1]+p[2]*ukr)*(pol_pop)^(p[3])
  phi <- p[4]
  res <- actuar::dztnbinom(x  = border, size = phi, prob = phi/(phi+mu), log = TRUE)
  -sum(res)
}

tnegbin_ll_sex <- function(p) {
  mu <- pop^(p[1]+p[2]*sex)*(pol_pop)^(p[3])
  phi <- p[4]
  res <- actuar::dztnbinom(x  = border, size = phi, prob = phi/(phi+mu), log = TRUE)
  -sum(res)
}

tnegbin_ll_both <- function(p) {
  mu <- pop^(p[1]+p[2]*ukr+p[3]*sex)*(pol_pop)^(p[4])
  phi <- p[5]
  res <- actuar::dztnbinom(x  = border, size = phi, prob = phi/(phi+mu), log = TRUE)
  -sum(res)
}

