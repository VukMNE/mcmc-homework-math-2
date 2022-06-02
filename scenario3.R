library(mvtnorm)
library(ggplot2)
library(gridExtra)

dataset <- read.csv('datset.csv')


model <- glm(y ~ X2,family=binomial(link='logit'),data=dataset)
px <- 1 / (1 + exp(- 1* ( - 0.715 * dataset$X2  -0.387 * dataset$X3)))
y <- dataset$y
L <- px^y * (1 - px)^(1-y)
log(prod(L))


a <- seq(1, 5, 1)
b <- seq(1, 5, 1)
te <- a %*% t(b)

li_df <- data.frame()

pdf_log <- function(b, x) {
  exp(b * x) / (1 + exp(b * x)) 
}

likelihood_regression_2cols <- function(betas) {
  # print(betas)
  X <- as.matrix(dataset[, c(1,2)])
  Y <- dataset$y
  B <- matrix(rep(betas, nrow(X)), nrow=nrow(X), byrow = T)
  p_of_xs <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  prod(p_of_xs^Y * (1 - p_of_xs)^(1-Y))
}

likelihood_regression_2cols(c(1.42, -0.715))

minus_log_likelihood_2cols <- function(betas) {
  X <- as.matrix(dataset[, c(1,2)])
  Y <- dataset$y
  B <- matrix(rep(betas, nrow(X)), nrow=nrow(X), byrow = T)
  p_of_xs <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  -1 * sum(Y * log(p_of_xs) + (1-Y) * log(1 - p_of_xs))
}

cov_logreg <- matrix(c(0.1, 0, 0, 0.1), ncol=2, nrow=2)

df <- metropolis_hastings(likelihood_regression_2cols, rmvnorm,  c(0, 0), 1000, cov_logreg)
x1 <- df[,1]
x2 <- df[,2]
df <- data.frame(x1,x2)

ggplot(df, aes(x = x1, y = x2)) + 
  geom_density_2d_filled(alpha = 0.4) +
  geom_density_2d(color='black') +
  theme(
    axis.text.x = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.title.y = element_text(size=16)
  )
summary(df)

bn_df <- vector()
x <- c(0, 0)
for (i in 1:1000) {
  x_n <- banana_proposal(1, x, cov_banana)
  bn_df <- rbind(bn_df, x_n)
  x <- x_n
}
bn_df <- as.data.frame(bn_df)

# rejection sampling

rejection_sampling <- function(m, target_dens_func, renvelope, denvelope, c) {
  # m is the number of samples that we want
  # target_dens_func needless to say, the pdf of target distribution
  # renvelope is a function for generating samples from envelope distribution
  samples <- vector()
  for (i in 1:m) {
    print(i)
    # print('*******************')
    cand_valid <- F
    while(cand_valid == F) {
      x_cand <- renvelope(1)
      u <- runif(1, 0, denvelope(x_cand, c))
      print(paste('dx: ', denvelope(x_cand, c)))
      print(paste('u: ', u))
      print(paste('p(x): ', target_dens_func(x_cand)))
      cand_valid <- u <= target_dens_func(x_cand)
      if (cand_valid) {
        samples <- rbind(samples, x_cand)
      }
    }
  }
  
  samples
}


renvelope_logistic_2col <- function(n) {
  x <- rnorm(n, 1, 0.5)
  y <- rnorm(n)
  c(x,y)
}

denvelope_logistic_2col <- function(x, c) {
  dx <- dnorm(x[1], 1, 0.5)
  dy <- dnorm(x[2])
  dx * dy * c
}

rj <- rejection_sampling(1000, likelihood_regression_2cols, renvelope_logistic_2col, denvelope_logistic_2col, 1e-220)
rj <- as.data.frame(rj)
ggplot(rj, aes(x = V1, y = V2)) +
  geom_density_2d_filled(alpha = 0.4) +
  geom_density_2d(color='black') +
  theme(
    axis.text.x = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.title.y = element_text(size=16)
  )
summary(rj)

rej_sampl_chains_un <- vector(mode = "list", length = 5)
names(rej_sampl_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
rj_times <- rep(0, 5)
for (i in 1:5) {
  rj_times[i] <- system.time(rej_sampl_chains_un[[paste('chain_', i, sep = '')]] <- rejection_sampling(1000, dmvnorm, renvelope_for_bivariate, denvelope_for_bivariate, 2.5))[3]
}

rjdf <- data.frame(mu1 = c(), mu2=c(), samples=c(), chain=c())
for (i in 1:5) {
  tv1 <- traceplot(rej_sampl_chains_un[[paste('chain_', i, sep = '')]][,1], mean)
  tv2 <- traceplot(rej_sampl_chains_un[[paste('chain_', i, sep = '')]][,2], mean)
  rjdf <- rbind(rjdf, data.frame(mu1= tv1, mu2=tv2, samples=1:1000, chain=rep(as.character(i), length(tv1))))
}

p3 <- ggplot(rjdf) + geom_line(aes(x=samples, y=mu1, color=chain))
p4 <- ggplot(rjdf) + geom_line(aes(x=samples, y=mu2, color=chain))
grid.arrange(p1, p2, p3, p4, ncol=2, nrow=2)


# HMC
source("HMC.r")
source("Multiplot.r")
library(ggplot2)
library(numDeriv)
library(coda)
library(grid)




pdf_log <- function(b, x) {
  exp(b * x) / (1 + exp(b * x)) 
}


minus_log_likelihood_2cols <- function(betas) {
  X <- as.matrix(dataset[, c(1,2)])
  Y <- dataset$y
  B <- matrix(rep(betas, nrow(X)), nrow=nrow(X), byrow = T)
  p_of_xs <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  -1 * sum(Y * log(p_of_xs) + (1-Y) * log(1 - p_of_xs))
}

minus_log_like_grad <- function(betas) {
  X <- as.matrix(dataset[, c(1,2)])
  Y <- dataset$y
  B <- matrix(rep(betas, nrow(X)), nrow=nrow(X), byrow = T)
  p_of_xs <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  grads <- vector()
  for(j in ncol(X)) {
    grads <- c(grads, -1 * sum((Y - p_of_xs) * X[,j]))
  }
  grads
}

## HMC
L = 50
epsilon = 0.01
current_q = c(0,0)

samples <- hmc_make_samples(minus_log_likelihood_2cols, minus_log_like_grad, epsilon, L, current_q, 1000)
minus_log_likelihood_2cols(c(-66.75890, -80.85702))
# example 2 - banana-shaped distribution
B <- 0.05
minus_log_likelihood(c(0.1, 0.2))


hmc_make_samples <- function(mlogU, grad_mlogU, epsilon, L, current_q, m) {
  samples <- NULL
  for (i in 1:m) {
    print(i)
    res = HMC(mlogU, grad_mlogU, epsilon, L, current_q)
    samples = rbind(samples, data.frame(Q1 = res$next_q[1], Q2 = res$next_q[2]))
    current_q = res$next_q
  }
  samples
}


ggplot(samples, aes(x = Q1, y = Q2)) + geom_point() +
  geom_density_2d_filled(alpha = 0.4) +
  geom_density_2d(color='black') +
  theme(
    axis.text.x = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.title.y = element_text(size=16)
  )

scenario3_run_chains <- function() {
  mh_chains_un <- vector(mode = "list", length = 5)
  mh_times <- rep(0, 5)
  names(mh_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
  cov_logreg <- matrix(c(0.1, 0, 0, 0.1), ncol=2, nrow=2)
  df <- vector()
  for (i in 1:5) {
    mh_times[i] <- system.time(mh_chains_un[[paste('chain_', i, sep = '')]] <- metropolis_hastings(likelihood_regression_2cols, rmvnorm,  c(0, 0), 1000, cov_logreg))[3]
    df <- rbind(df, cbind(as.data.frame(mh_chains_un[[paste('chain_', i, sep = '')]]), samples=1:1000,chain=rep(as.character(i), 1000)  ))
  }
  
  df <- as.data.frame(df)
  p1 <- ggplot(df) + geom_line(aes(x=samples, y=V1, color=chain)) + theme(legend.position = "none") 
  p2 <- ggplot(df) + geom_line(aes(x=samples, y=V2, color=chain)) + theme(legend.position = "none")
  
  L = 50
  epsilon = 0.01
  current_q = c(0,0)
  m = 1000
  
  hmc_make_samples <- function(mlogU, grad_mlogU, epsilon, L, current_q, m) {
    samples <- NULL
    for (i in 1:m) {
      print(i)
      res = HMC(mlogU, grad_mlogU, epsilon, L, current_q)
      samples = rbind(samples, data.frame(Q1 = res$next_q[1], Q2 = res$next_q[2]))
      current_q = res$next_q
    }
    samples
  }
  
  minus_log_likelihood_2cols <- function(betas) {
    X <- as.matrix(dataset[, c(1,2)])
    Y <- dataset$y
    B <- matrix(rep(betas, nrow(X)), nrow=nrow(X), byrow = T)
    p_of_xs <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
    -1 * sum(Y * log(p_of_xs) + (1-Y) * log(1 - p_of_xs))
  }
  
  minus_log_like_grad <- function(betas) {
    X <- as.matrix(dataset[, c(1,2)])
    Y <- dataset$y
    B <- matrix(rep(betas, nrow(X)), nrow=nrow(X), byrow = T)
    p_of_xs <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
    grads <- vector()
    for(j in ncol(X)) {
      grads <- c(grads, -1 * sum((Y - p_of_xs) * X[,j]))
    }
    grads
  }
  
  hmc_chains_un <- vector(mode = "list", length = 5)
  names(hmc_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
  hmc_times <- rep(0, 5)
  hmc_df <- vector()
  for (i in 1:5) {
    hmc_times[i] <- system.time(hmc_chains_un[[paste('chain_', i, sep = '')]] <- hmc_make_samples(minus_log_likelihood_2cols, minus_log_like_grad, epsilon, L, current_q, 1000))[3]
    hmc_df <- rbind(hmc_df, cbind(as.data.frame(hmc_chains_un[[paste('chain_', i, sep = '')]]), samples=1:1000,chain=rep(as.character(i), 1000)  )) 
  }
  
  hmc_df <- as.data.frame(hmc_df)
  p3 <- ggplot(hmc_df) + geom_line(aes(x=samples, y=Q1, color=chain)) + theme(legend.position = "none")
  p4 <- ggplot(hmc_df) + geom_line(aes(x=samples, y=Q2, color=chain)) + theme(legend.position = "none")
  
  
  
  rej_sampl_chains_un <- vector(mode = "list", length = 5)
  names(rej_sampl_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
  rj_times <- rep(0, 5)
  rjdf <- vector()
  for (i in 1:5) {
    rj_times[i] <- system.time(rej_sampl_chains_un[[paste('chain_', i, sep = '')]] <- rejection_sampling(1000, likelihood_regression_2cols, renvelope_logistic_2col, denvelope_logistic_2col, 1e-220))[3]
    rjdf <- rbind(rjdf, cbind(as.data.frame(rej_sampl_chains_un[[paste('chain_', i, sep = '')]]), samples=1:1000,chain=rep(as.character(i), 1000)  )) 
  }
  
  rjdf <- as.data.frame(rjdf)
  p5 <- ggplot(rjdf) + geom_line(aes(x=samples, y=V1, color=chain)) + theme(legend.position = "none")
  p6 <- ggplot(rjdf) + geom_line(aes(x=samples, y=V2, color=chain)) + theme(legend.position = "none")
  
  
  grid.arrange(p1, p3, p5, p2, p4, p6, ncol=3, nrow=2)
}

mh_df <- df
mh_df$samples <- NULL
mh_df$chain <- NULL
plot(acf(mh_df))

ac_df <- hmc_df
ac_df$samples <- NULL
ac_df$chain <- NULL
plot(acf(ac_df))


ac_df <- rjdf
ac_df$samples <- NULL
ac_df$chain <- NULL
plot(acf(ac_df))


mh_df <- df
mh_df$samples <- NULL
mh_df$chain <- NULL
plot(acf(mh_df))

ac_df <- hmc_df
ac_df$b1 <- ac_df$mu1
ac_df$b2 <- ac_df$mu2
ac_df$mu1 <- NULL
ac_df$mu2 <- NULL
ac_df$samples <- NULL
ac_df$chain <- NULL
plot(acf(ac_df))


ac_df <- rjdf
ac_df$b1 <- ac_df$mu1
ac_df$b2 <- ac_df$mu2
ac_df$mu1 <- NULL
ac_df$mu2 <- NULL
ac_df$samples <- NULL
ac_df$chain <- NULL
plot(acf(ac_df))


library(mcmcse)

df$alg <- 'Metropolis-Hastings'
rjdf$alg <- 'Rejection-Sampling'
hmc_df$alg <- 'HMC'

hmc_df$V1 <- hmc_df$Q1
hmc_df$V2 <- hmc_df$Q2
hmc_df$Q1 <- NULL
hmc_df$Q2 <- NULL

df$ess_mu1 <- ess(df$V1)
rjdf$ess_mu1 <- ess(rjdf$V1)
hmc_df$ess_mu1 <- ess(hmc_df$V1)
ess_df <- rbind(df, rjdf, hmc_df)
scen2_ess1 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_mu1, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS')

df$ess_mu2 <- ess(df$V2)
rjdf$ess_mu2 <- ess(rjdf$V2)
hmc_df$ess_mu2 <- ess(hmc_df$V2)
ess_df <- rbind(df, rjdf, hmc_df)
scen2_ess2 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_mu2, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS')

summary(rjdf)

ess_df$ess_per_sec_mu1 <- 0
ess_df$ess_per_sec_mu2 <- 0
for (j in 1:5) {
  ess_df[ess_df$alg == 'Metropolis-Hastings' & ess_df$chain == as.character(j),]$ess_per_sec_mu1 <- ess_df[ess_df$alg == 'Metropolis-Hastings' & ess_df$chain == as.character(j),]$ess_mu1 / mh_times[j]
  ess_df[ess_df$alg == 'Metropolis-Hastings' & ess_df$chain == as.character(j),]$ess_per_sec_mu2 <- ess_df[ess_df$alg == 'Metropolis-Hastings' & ess_df$chain == as.character(j),]$ess_mu2 / mh_times[j]
  ess_df[ess_df$alg == 'HMC' & ess_df$chain == as.character(j),]$ess_per_sec_mu1 <- ess_df[ess_df$alg == 'HMC' & ess_df$chain == as.character(j),]$ess_mu1 / hmc_times[j]
  ess_df[ess_df$alg == 'HMC' & ess_df$chain == as.character(j),]$ess_per_sec_mu2 <- ess_df[ess_df$alg == 'HMC' & ess_df$chain == as.character(j),]$ess_mu2 / hmc_times[j]
  ess_df[ess_df$alg == 'Rejection-Sampling' & ess_df$chain == as.character(j),]$ess_per_sec_mu1 <- ess_df[ess_df$alg == 'Rejection-Sampling' & ess_df$chain == as.character(j),]$ess_mu1 / rj_times[j]
  ess_df[ess_df$alg == 'Rejection-Sampling' & ess_df$chain == as.character(j),]$ess_per_sec_mu2 <- ess_df[ess_df$alg == 'Rejection-Sampling' & ess_df$chain == as.character(j),]$ess_mu2 / rj_times[j]
  
}
scen2_ess2_sec1 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_per_sec_mu1, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS')

scen2_ess2_sec2 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_per_sec_mu2, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS')

grid.arrange(scen2_ess1, scen2_ess1, scen2_ess2_sec1, scen2_ess2_sec2, ncol=2, nrow=2)
summary(rjdf)

ess_df$ess_per_sec_b1 <- 0
ess_df$ess_per_sec_b2 <- 0
for (j in 1:5) {
  ess_df[ess_df$alg == 'Metropolis-Hastings' & ess_df$chain == as.character(j),]$ess_per_sec_b1 <- ess_df[ess_df$alg == 'Metropolis-Hastings' & ess_df$chain == as.character(j),]$ess_b1 / mh_times[j]
  ess_df[ess_df$alg == 'Metropolis-Hastings' & ess_df$chain == as.character(j),]$ess_per_sec_b2 <- ess_df[ess_df$alg == 'Metropolis-Hastings' & ess_df$chain == as.character(j),]$ess_b2 / mh_times[j]
  ess_df[ess_df$alg == 'HMC' & ess_df$chain == as.character(j),]$ess_per_sec_b1 <- ess_df[ess_df$alg == 'HMC' & ess_df$chain == as.character(j),]$ess_b1 / hmc_times[j]
  ess_df[ess_df$alg == 'HMC' & ess_df$chain == as.character(j),]$ess_per_sec_b2 <- ess_df[ess_df$alg == 'HMC' & ess_df$chain == as.character(j),]$ess_b2 / hmc_times[j]
  ess_df[ess_df$alg == 'Rejection-Sampling' & ess_df$chain == as.character(j),]$ess_per_sec_b1 <- ess_df[ess_df$alg == 'Rejection-Sampling' & ess_df$chain == as.character(j),]$ess_b1 / rj_times[j]
  ess_df[ess_df$alg == 'Rejection-Sampling' & ess_df$chain == as.character(j),]$ess_per_sec_b2 <- ess_df[ess_df$alg == 'Rejection-Sampling' & ess_df$chain == as.character(j),]$ess_b2 / rj_times[j]
  
}
p_ess_sec_b1 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_per_sec_b1, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS')

p_ess_sec_b2 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_per_sec_b2, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS')

mcse(df$V1)
mcse(df$V2)

mcse(hmc_df$V1)
mcse(hmc_df$V2)

mcse(rjdf$V1)
mcse(rjdf$V2)