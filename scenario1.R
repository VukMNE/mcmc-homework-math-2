library(mvtnorm)
library(ggplot2)
library(gridExtra)
source('algorithms.R')

mh_chains_un <- vector(mode = "list", length = 5)
names(mh_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
cov_matrix_prop <- matrix(c(1, 0, 0, 1), ncol=2, nrow=2)
mh_times <- rep(0, 5)

df <- vector()
for (i in 1:5) {
  mh_times[i] <- system.time(mh_chains_un[[paste('chain_', i, sep = '')]] <- metropolis_hastings(dmvnorm, rmvnorm,  c(0, 0), 1000, cov_matrix_prop))[3]
  df <- rbind(df, cbind(as.data.frame(mh_chains_un[[paste('chain_', i, sep = '')]]), samples=1:1000,chain=rep(as.character(i), 1000)  )) 
  
}
df <- as.data.frame(df)
p1 <- ggplot(df) + geom_line(aes(x=samples, y=V1, color=chain)) + theme(legend.position = "none") 
p2 <- ggplot(df) + geom_line(aes(x=samples, y=V2, color=chain)) + theme(legend.position = "none")



# HMC 
L = 1
epsilon = 0.25
current_q = c(0,0)
m = 1000

minus_log_bivar_st_norm <- function(x) {
  log(2 * pi) + 0.5 * sum(x * t(x))
}

grad_minus_log_bivar_st_norm <- function(x) {
  x
}



hmc_chains_un <- vector(mode = "list", length = 5)
names(hmc_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
hmc_times <- rep(0, 5)
hmc_df <- vector()
for (i in 1:5) {
  hmc_times[i] <- system.time(hmc_chains_un[[paste('chain_', i, sep = '')]] <-  hmc_make_samples(minus_log_bivar_st_norm, grad_minus_log_bivar_st_norm, epsilon, L, current_q, 1000))[3]
  hmc_df <- rbind(hmc_df, cbind(as.data.frame(hmc_chains_un[[paste('chain_', i, sep = '')]]), samples=1:1000,chain=rep(as.character(i), 1000)  )) 
}
hmc_df <- as.data.frame(df)
p3 <- ggplot(hmc_df) + geom_line(aes(x=samples, y=V1, color=chain)) + theme(legend.position = "none")
p4 <- ggplot(hmc_df) + geom_line(aes(x=samples, y=V2, color=chain)) + theme(legend.position = "none")

# rejection sampling

renvelope_for_bivariate <- function(n) {
  x <- rbeta(n, 2, 2) * 8 - 4
  y <- rbeta(n, 2, 2) * 8 - 4
  c(x,y)
}

denvelope_for_bivariate <- function(x, c) {
  # x - an instance we need to evaulate
  # c - constant multiplier
  dx <- dbeta((x[1] + 4) / 8, 2, 2)
  dy <- dbeta((x[2] + 4) / 8, 2, 2)
  dx * dy * c
}

rej_sampl_chains_un <- vector(mode = "list", length = 5)
names(rej_sampl_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
rj_times <- rep(0, 5)
rjdf <- vector()
for (i in 1:5) {
  rj_times[i] <- system.time(rej_sampl_chains_un[[paste('chain_', i, sep = '')]] <- rejection_sampling(1000, dmvnorm, renvelope_for_bivariate, denvelope_for_bivariate, 2.5))[3]
  rjdf <- rbind(rjdf, cbind(as.data.frame(rej_sampl_chains_un[[paste('chain_', i, sep = '')]]), samples=1:1000,chain=rep(as.character(i), 1000)  )) 
}
rjdf <- as.data.frame(rjdf)
p5 <- ggplot(rjdf) + geom_line(aes(x=samples, y=V1, color=chain)) + theme(legend.position = "none")
p6 <- ggplot(rjdf) + geom_line(aes(x=samples, y=V2, color=chain)) + theme(legend.position = "none") 


grid.arrange(p1, p3, p5, p2, p4, p6, ncol=3, nrow=2)

# DRAWING PLOTS FOR AUTO-CORRELATION, ESS and ESS/second

acf(df$V1)
acf(ac_df)
acf(rjdf$V1)

ac_df <- data.frame(V1=hmc_df$V1, V2=hmc_df$V2)
library(mcmcse)

df$alg <- 'Metropolis-Hastings'
rjdf$alg <- 'Rejection-Sampling'
hmc_df$alg <- 'HMC'

df$ess_mu1 <- ess(df$V1)
rjdf$ess_mu1 <- ess(rjdf$V1)
hmc_df$ess_mu1 <- ess(hmc_df$V1)
ess_df <- rbind(df, rjdf, hmc_df)
s1_ess <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_mu1, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS') + theme(legend.position = "none")

df$ess_mu2 <- ess(df$V2)
rjdf$ess_mu2 <- ess(rjdf$V2)
hmc_df$ess_mu2 <- ess(hmc_df$V2)
ess_df <- rbind(df, rjdf, hmc_df)
s1_ess2 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_mu2, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS')  + theme(legend.position = "none")

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
ess_sec1 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_per_sec_mu1, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS')  + theme(legend.position = "none")

ess_sec2 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_per_sec_mu2, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS')  + theme(legend.position = "none")

grid.arrange(s1_ess, s1_ess2, ess_sec1, ess_sec2, ncol=2, nrow=2)