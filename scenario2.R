source("HMC.r")
source("Multiplot.r")
library(ggplot2)
library(numDeriv)
library(coda)
library(grid)
source('algorithms.R')

B <- 0.05

banana_func <- function(x) {
  exp(-(x[1]^2)/200 - 0.5 * (x[2] + 0.05 * x[1] - 5)^2)
}

# example 2 - banana-shaped distribution


banana_proposal <- function(n, x, cov_mat) {
  z <- rmvnorm(n, c(0, 0), cov_mat)
  x[1] <- z[1]
  x[2] <- z[2] + B * z[1]^2 - 100 * B
  x
}

cov_banana <- matrix(c(100, 0, 0, 1), ncol=2, nrow=2)


mh_chains_un <- vector(mode = "list", length = 5)
mh_times <- rep(0, 5)
names(mh_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
df <- vector()
for (i in 1:5) {
  mh_times[i] <- system.time(mh_chains_un[[paste('chain_', i, sep = '')]] <- metropolis_hastings(banana_func, banana_proposal,  c(0, 0), 1000, cov_banana))[3]
  df <- rbind(df, cbind(as.data.frame(mh_chains_un[[paste('chain_', i, sep = '')]]), samples=1:1000,chain=rep(as.character(i), 1000)  )) 
}

df <- as.data.frame(df)
p1 <- ggplot(df) + geom_line(aes(x=samples, y=V1, color=chain)) + theme(legend.position = "none")
p2 <- ggplot(df) + geom_line(aes(x=samples, y=V2, color=chain)) + theme(legend.position = "none")

# HMC 
L = 27
epsilon = 0.6
current_q = c(0,0)
m = 1000

minus_logf <- function(x) {
  -(-(x[1]^2)/200- 0.5 * (x[2]+ B * x[1]^2 - 100*B)^2 )
}

minus_logf_grad <- function(x) {
  g1 <- -(x[1])/100- 1.0 * (2* B * x[1]) * (x[2]+ B * x[1]^2 - 100*B)
  g2 <- - 1.0 * (x[2]+ B * x[1]^2 - 100*B)
  -c(g1,g2)
}


hmc_chains_un <- vector(mode = "list", length = 5)
names(hmc_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
hmc_times <- rep(0, 5)
hmc_df <- vector()
for (i in 1:5) {
  hmc_times[i] <- system.time(hmc_chains_un[[paste('chain_', i, sep = '')]] <-  hmc_make_samples(minus_logf, minus_logf_grad, epsilon, L, current_q, 1000))[3]
  hmc_df <- rbind(hmc_df, cbind(as.data.frame(hmc_chains_un[[paste('chain_', i, sep = '')]]), samples=1:1000,chain=rep(as.character(i), 1000)  )) 
}
hmc_df <- as.data.frame(hmc_df)

p3 <- ggplot(hmc_df) + geom_line(aes(x=samples, y=Q1, color=chain)) + theme(legend.position = "none")
p4 <- ggplot(hmc_df) + geom_line(aes(x=samples, y=Q2, color=chain)) + theme(legend.position = "none")


# REJECTION SAMPLING


renvelope_banana <- function(n) {
  # x <- rbeta(n, 2, 2) * 40 - 20
  # y <- rbeta(n, 2, 2) * 50 - 40
  x <- runif(n, -30, 30)
  y <- runif(n, -40, 10)
  c(x,y)
}

denvelope_banana <- function(x, c) {
  # x - an instance we need to evaulate
  # c - constant multiplier
  # dx <- dbeta((x[1] + 20) / 40, 2, 2)
  # dy <- dbeta((x[2] + 40) / 50, 2, 2)
  dx <- dunif(x, -30, 30) / 60
  dy <- dunif(x, -40, 10) / 50
  dx * dy * c
}

rej_sampl_chains_un <- vector(mode = "list", length = 5)
names(rej_sampl_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
rj_times <- rep(0, 5)
rjdf <- vector()
for (i in 1:5) {
  rj_times[i] <- system.time(rej_sampl_chains_un[[paste('chain_', i, sep = '')]] <- rejection_sampling(1000, banana_func, renvelope_banana, denvelope_banana, 1e-5))[3]
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
grid.arrange(p1, p3, p5, p2, p4, p6, ncol=3, nrow=2, top="Main Title")


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
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS') + theme(legend.position = "none")

df$ess_mu2 <- ess(df$V2)
rjdf$ess_mu2 <- ess(rjdf$V2)
hmc_df$ess_mu2 <- ess(hmc_df$V2)
ess_df <- rbind(df, rjdf, hmc_df)
scen2_ess2 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_mu2, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS') + theme(legend.position = "none")

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
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS') + theme(legend.position = "none")

scen2_ess2_sec2 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_per_sec_mu2, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS') + theme(legend.position = "none")

grid.arrange(scen2_ess1, scen2_ess1, scen2_ess2_sec1, scen2_ess2_sec2, ncol=2, nrow=2)

mcse(df$V1)
mcse(df$V2)

mcse(hmc_df$V1)
mcse(hmc_df$V2)

mcse(rjdf$V1)
mcse(rjdf$V2)