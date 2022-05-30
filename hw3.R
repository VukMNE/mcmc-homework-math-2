metropolis_hastings <- function(p, k, start_state, m, cov_matrix) {
  # p is the target distribution function
  # k is the proposal distribution function
  # start_state is starting state
  # m is the number of samples we want to draw
  
  delta <- rmvnorm(1, rep(0, 2), cov_matrix)
  x_new <- x + delta
  
  
}

cov_matrix <- matrix(c(1, 0, 0, 1), ncol=2, nrow=2)


# beta density M-H
rule_D <- function(x, alpha, beta, delta = 0.1) {
  x_new <- runif(1, x - delta, x + delta)
  if (x_new > 1 | x_new < 0) return (x)
  p <- exp(
    dbeta(x_new, alpha, beta, log = T) - 
      dbeta(x, alpha, beta, log = T))
  # okej, ovaj dio je isti
  if (runif(1) > p) x_new <- x # reject
  x_new
}
