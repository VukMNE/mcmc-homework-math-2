library(mvtnorm)
library(ggplot2)
library(gridExtra)

metropolis_hastings <- function(p, k, start_state, m, cov_matrix) {
  # p is the target distribution function
  # k is the proposal distribution function
  # start_state is starting state
  # m is the number of samples we want to draw
  x <- start_state
  ndims <- nrow(cov_matrix)
  samples <- x
  
  for (i in 1:(m-1)) {
    # k = rmvnorm
    # delta <- k(1, rep(0, nrow(cov_matrix)), cov_matrix)
    x_new <- k(1, x, cov_matrix)
    alpha <- min(1, p(x_new) / p(x))
    if(is.nan(alpha)) {
      alpha <- 0
    }
    u <- runif(1)
    if(u > alpha) {
      x_new <- x
    } 
    samples <- rbind(samples, x_new)
    x <- x_new
    
  }
  samples
  
}

# Rejection sampling

rejection_sampling <- function(m, target_dens_func, renvelope, denvelope, c) {
  # m is the number of samples that we want
  # target_dens_func needless to say, the pdf of target distribution
  # renvelope is a function for generating samples from envelope distribution
  samples <- vector()
  for (i in 1:m) {
    # print(i)
    # print('*******************')
    cand_valid <- F
    while(cand_valid == F) {
      x_cand <- renvelope(1)
      u <- runif(1, 0, denvelope(x_cand, c))
      # print(x_cand)
      # print(target_dens_func(x_cand))
      # print(u <= target_dens_func(x_cand))
      # print('-----------------------')
      cand_valid <- u <= target_dens_func(x_cand)
      if (cand_valid) {
        samples <- rbind(samples, x_cand)
      }
    }
  }
  
  samples
}


# CODE for HMC is in HMC.r file 

source("HMC.r")
source("Multiplot.r")
library(ggplot2)
library(numDeriv)
library(coda)
library(grid)


#HMC 
hmc_make_samples <- function(mlogU, grad_mlogU, epsilon, L, current_q, m) {
  samples <- NULL
  for (i in 1:m) {
    res = HMC(mlogU, grad_mlogU, epsilon, L, current_q)
    samples = rbind(samples, data.frame(Q1 = res$next_q[1], Q2 = res$next_q[2]))
    current_q = res$next_q
  }
  samples
}
