library(mvtnorm)
library(ggplot2)
library(gridExtra)
source('algorithms.R')


dataset <- read.csv('datset.csv')
# fitting the model to see best parameters
model <- glm(y ~ . - X1,family=binomial(link='logit'),data=dataset)

likelihood_dataset <- function(betas) {
  # print(betas)
  X <- as.matrix(dataset[, 1:ncol(dataset)-1])
  Y <- dataset$y
  B <- matrix(rep(betas, nrow(X)), nrow=nrow(X), byrow = T)
  p_of_xs <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  prod(p_of_xs^Y * (1 - p_of_xs)^(1-Y))
}

mh_chains_un <- vector(mode = "list", length = 5)
mh_times <- rep(0, 5)
names(mh_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
cov_dataset <- diag(11) * 0.1

df <- vector()
for (i in 1:5) {
  mh_times[i] <- system.time(mh_chains_un[[paste('chain_', i, sep = '')]] <- metropolis_hastings(likelihood_dataset, rmvnorm,  rep(0, 11), 1000, cov_dataset))[3]
  df <- rbind(df, cbind(as.data.frame(mh_chains_un[[paste('chain_', i, sep = '')]]), samples=1:1000,chain=rep(as.character(i), 1000)  )) 
}
df <- as.data.frame(df)


p1 <- ggplot(df) + geom_line(aes(x=samples, y=V1, color=chain)) + theme(legend.position = "none")
p2 <- ggplot(df) + geom_line(aes(x=samples, y=V2, color=chain)) + theme(legend.position = "none")


# HMC 

minus_log_likelihood_dataset <- function(betas) {
  X <- as.matrix(dataset[, 1:ncol(dataset)-1])
  Y <- dataset$y
  B <- matrix(rep(betas, nrow(X)), nrow=nrow(X), byrow = T)
  p_of_xs <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  -1 * sum(Y * log(p_of_xs) + (1-Y) * log(1 - p_of_xs))
}

minus_log_like_grad <- function(betas) {
  X <- as.matrix(dataset[, 1:ncol(dataset)-1])
  Y <- dataset$y
  B <- matrix(rep(betas, nrow(X)), nrow=nrow(X), byrow = T)
  p_of_xs <-  exp(rowSums(B * X)) / (1 + exp(rowSums(B * X)))
  grads <- vector()
  for(j in 1:ncol(X)) {
    grads <- c(grads, -1 * sum((Y - p_of_xs) * X[,j]))
  }
  grads
}



hmc_make_samples <- function(mlogU, grad_mlogU, epsilon, L, current_q, m) {
  samples <- NULL
  for (i in 1:m) {
    print(i)
    res = HMC(mlogU, grad_mlogU, epsilon, L, current_q)
    samples = rbind(samples, data.frame(Q1 = res$next_q[1],
                                        Q2 = res$next_q[2],
                                        Q3 = res$next_q[3],
                                        Q4 = res$next_q[4],
                                        Q5 = res$next_q[5],
                                        Q6 = res$next_q[6],
                                        Q7 = res$next_q[7],
                                        Q8 = res$next_q[8],
                                        Q9 = res$next_q[9],
                                        Q10 = res$next_q[10],
                                        Q11 = res$next_q[11]
    ))
    #print(summary(res$next_q))
    current_q = res$next_q
  }
  samples
}


L = 50
epsilon = 0.01
current_q = rep(0, 11)
m = 1000

hmc_chains_un <- vector(mode = "list", length = 5)
names(hmc_chains_un) <- c('chain_1', 'chain_2', 'chain_3', 'chain_4', 'chain_5')
hmc_times <- rep(0, 5)
hmc_df <- vector()
for (i in 1:5) {
  hmc_times[i] <- system.time(hmc_chains_un[[paste('chain_', i, sep = '')]] <-   hmc_make_samples(minus_log_likelihood_dataset, minus_log_like_grad, epsilon, L, current_q, 1000))[3]
  hmc_df <- rbind(hmc_df, cbind(as.data.frame(hmc_chains_un[[paste('chain_', i, sep = '')]]), samples=1:1000,chain=rep(as.character(i), 1000)  )) 
}
hmc_df <- as.data.frame(hmc_df)
names(hmc_df) <- names(df)


p3 <- ggplot(hmc_df) + geom_line(aes(x=samples, y=V1, color=chain)) + theme(legend.position = "none")
p4 <- ggplot(hmc_df) + geom_line(aes(x=samples, y=V2, color=chain)) + theme(legend.position = "none")





acf(df$V1)
acf(hmc_df$V2)

library(mcmcse)

df$alg <- 'Metropolis-Hastings'
hmc_df$alg <- 'HMC'

df$ess_b1 <- ess(df$V1)
hmc_df$ess_b1 <- ess(hmc_df$V1)

ess_df <- rbind(df, hmc_df)
p_ess1 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_b1, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS') + theme(legend.position = "none")



df$ess_b2 <- ess(df$V2)
hmc_df$ess_b2 <- ess(hmc_df$V2)
ess_df <- rbind(df, hmc_df)
p_ess2 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_b2, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS') + theme(legend.position = "none")

summary(rjdf)

ess_df$ess_per_sec_b1 <- 0
ess_df$ess_per_sec_b2 <- 0
for (j in 1:5) {
  ess_df[ess_df$alg == 'Metropolis-Hastings' & ess_df$chain == as.character(j),]$ess_per_sec_b1 <- ess_df[ess_df$alg == 'Metropolis-Hastings' & ess_df$chain == as.character(j),]$ess_b1 / mh_times[j]
  ess_df[ess_df$alg == 'Metropolis-Hastings' & ess_df$chain == as.character(j),]$ess_per_sec_b2 <- ess_df[ess_df$alg == 'Metropolis-Hastings' & ess_df$chain == as.character(j),]$ess_b2 / mh_times[j]
  ess_df[ess_df$alg == 'HMC' & ess_df$chain == as.character(j),]$ess_per_sec_b1 <- ess_df[ess_df$alg == 'HMC' & ess_df$chain == as.character(j),]$ess_b1 / hmc_times[j]
  ess_df[ess_df$alg == 'HMC' & ess_df$chain == as.character(j),]$ess_per_sec_b2 <- ess_df[ess_df$alg == 'HMC' & ess_df$chain == as.character(j),]$ess_b2 / hmc_times[j]

}
p_ess_sec_b1 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_per_sec_b1, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS') + theme(legend.position = "none")

p_ess_sec_b2 <- ggplot(ess_df) + 
  geom_bar(aes(factor(alg), ess_per_sec_b2, fill = chain), stat="identity", position = "dodge") +
  scale_fill_brewer(palette = "Set1") + xlab('') + ylab('ESS') + theme(legend.position = "none")

grid.arrange(p_ess1, p_ess2, p_ess_sec_b1, p_ess_sec_b2, ncol=2, nrow=2)

