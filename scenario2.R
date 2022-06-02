library(mvtnorm)
library(ggplot2)
library(gridExtra)
B <- 0.05

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


mh_chains_un <- vector(mode = "list", length = 5)
mh_times <- rep(0, 5)
names(mh_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
for (i in 1:5) {
  mh_times[i] <- system.time(mh_chains_un[[paste('chain_', i, sep = '')]] <- metropolis_hastings(banana_func, banana_proposal,  c(0.2, 0.3), 1000, cov_banana))[3]
}


traceplot <- function(x, func) {
  tp <- vector()
  print(paste('Length is:', length(x)))
  for (i in 1:length(x)) {
    tp <- c(tp, func(x[1:i]))
  }
  tp
}


cov_banana <- matrix(c(100, 0, 0, 1), ncol=2, nrow=2)

banana_proposal <- function(n, x, cov_mat) {
  z <- rmvnorm(n, c(0, 0), cov_mat)
  x[1] <- z[1]
  x[2] <- z[2] + B * z[1]^2 - 100 * B
  x
}
df <- metropolis_hastings(banana_func, banana_proposal,  c(0, 0), 1000, cov_banana)
x <- df[,1]
y <- df[,2]
df <- data.frame(x,y)

ggplot(df, aes(x = x, y = y)) + 
  geom_density_2d_filled(alpha = 0.4) +
  geom_density_2d(color='black') +
  theme(
    axis.text.x = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.title.y = element_text(size=16)
  )
mean(df$x)
mean(df$y)


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
    # print(i)
    # print('*******************')
    cand_valid <- F
    while(cand_valid == F) {
      x_cand <- renvelope(1)
      u <- runif(1, 0, denvelope(x_cand, c))
      cand_valid <- u <= target_dens_func(x_cand)
      if (cand_valid) {
        samples <- rbind(samples, x_cand)
      }
    }
  }
  
  samples
}


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


rj <- rejection_sampling(1000, banana_func, renvelope_banana, denvelope_banana, 1e-5)
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







u <- rbeta(1000, 2, 2) * 8 - 4
n <- rnorm(1000)
x <- seq(-4, 4, length.out = 512)
du <- density(u)$y * 2.5
dn <- density(n)$y
comp_df <- data.frame(x, du, dn)
cdf <- data.frame(u, n)
ggplot(comp_df) + geom_line(aes(x=x, y=du), color='blue') + geom_density(data=cdf, aes(x=n), color='red')

comp_df <- data.frame(u, xu, du, n)
ggplot(comp_df) + geom_density(aes(x=xu), color='blue') + geom_density(aes(x=n), color='red')

u <- rbeta(1, 2, 2) * 8 - 4
du <- pbeta((u + 4) / 8, 2, 2)
summary(rbeta(100, 2, 2))


x <- rbeta(1000, 2, 2) * 8 - 4
y <- rbeta(1000, 2, 2) * 8 - 4
df <- data.frame(x, y)

ggplot(df, aes(x = x, y = y)) +
  geom_density_2d_filled(alpha = 0.4) +
  geom_density_2d(color='black') +
  theme(
    axis.text.x = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.title.y = element_text(size=16)
  )

xd <- rmvnorm(1000)
df$x <- xd[,1]
df$x <- xd[,y]
ggplot(df, aes(x = x, y = y)) +
  geom_density_2d_filled(alpha = 0.4) +
  geom_density_2d(color='black') +
  theme(
    axis.text.x = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.title.y = element_text(size=16)
  )



# HMC
source("HMC.r")
source("Multiplot.r")
library(ggplot2)
library(numDeriv)
library(coda)
library(grid)

minus_log_bivar_st_norm <- function(x) {
  log(2 * pi) + 0.5 * sum(x * t(x))
}

grad_minus_log_bivar_st_norm <- function(x) {
  x
}

banana_func <- function(x) {
  exp(-(x[1]^2)/200 - 0.5 * (x[2] + 0.05 * x[1] - 5)^2)
}

minus_logf <- function(x) {
  -(-(x[1]^2)/200- 0.5 * (x[2]+ B * x[1]^2 - 100*B)^2 )
}

minus_logf_grad <- function(x) {
  g1 <- -(x[1])/100- 1.0 * (2* B * x[1]) * (x[2]+ B * x[1]^2 - 100*B)
  g2 <- - 1.0 * (x[2]+ B * x[1]^2 - 100*B)
  -c(g1,g2)
}


# example 2 - banana-shaped distribution
B <- 0.05


## HMC
L = 27
epsilon = 0.6
current_q = c(0.33,0.33)
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

samples <- hmc_make_samples(minus_logf, minus_logf_grad, epsilon, L, current_q, 1000)

ggplot(samples, aes(x = Q1, y = Q2)) + geom_point() +
geom_density_2d_filled(alpha = 0.4) +
  geom_density_2d(color='black') +
  theme(
    axis.text.x = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.text.y = element_text(size=14),
    axis.title.y = element_text(size=16)
  )

scenario2_run_chains <- function() {
  cov_banana <- matrix(c(100, 0, 0, 1), ncol=2, nrow=2)
  
  banana_proposal <- function(n, x, cov_mat) {
    z <- rmvnorm(n, c(0, 0), cov_mat)
    x[1] <- z[1]
    x[2] <- z[2] + B * z[1]^2 - 100 * B
    x
  }
  
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
  
  
  L = 27
  epsilon = 0.6
  current_q = c(0,0)
  m = 1000
  
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