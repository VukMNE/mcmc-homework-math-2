library(mvtnorm)
library(ggplot2)
library(gridExtra)
source('algorithms.R')


dataset <- read.csv('datset.csv')


# fitting the model with glm
model <- glm(y ~ X2,family=binomial(link='logit'),data=dataset)
px <- 1 / (1 + exp(- 1* ( - 0.715 * dataset$X2  -0.387 * dataset$X3)))
y <- dataset$y
L <- px^y * (1 - px)^(1-y)
log(prod(L))


likelihood_regression_2cols <- function(betas) {
  # print(betas)
  X <- as.matrix(dataset[, c(1,2)])
  Y <- dataset$y
  B <- matrix(rep(betas, nrow(X)), nrow=nrow(X), byrow = T)
  p_of_xs <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  prod(p_of_xs^Y * (1 - p_of_xs)^(1-Y))
}

# checking the value of maximum likelihood
likelihood_regression_2cols(c(1.42, -0.715))

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


# REJECTION SAMPLING

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