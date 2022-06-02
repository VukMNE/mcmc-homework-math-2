source("HMC.r")
source("Multiplot.r")
library(ggplot2)
library(numDeriv)
library(coda)
library(grid)


# example 2 - banana-shaped distribution
B <- 0.05

minus_logf <- function(x) {
  -(-(x[1]^2)/200- 0.5 * (x[2]+ B * x[1]^2 - 100*B)^2 )
}

minus_logf_grad <- function(x) {
  g1 <- -(x[1])/100- 1.0 * (2* B * x[1]) * (x[2]+ B * x[1]^2 - 100*B)
  g2 <- - 1.0 * (x[2]+ B * x[1]^2 - 100*B)
  -c(g1,g2)
}


## HMC
L = 27
epsilon = 0.6
current_q = c(0,0)
m = 1000

samples <- NULL
pdf(paste("trajectories-ex02.pdf",sep=""), width = 9, height = 3)
for (i in 1:m) {
  print(i)
  res = HMC(minus_logf, minus_logf_grad, epsilon, L, current_q)
  samples = rbind(samples, data.frame(Q1 = res$next_q[1], Q2 = res$next_q[2]))
  current_q = res$next_q
  

}
dev.off()


pdf(paste("densities-ex02.pdf",sep=""), width = 8, height = 4)
# true density
x <- seq(-25,25,0.2)
x0 <- expand.grid(x,x)
y <- apply(x0,1,minus_logf)
df <- data.frame(x0,y = exp(-y))

gA <- ggplot(df, aes(Var1, Var2, z = y)) + geom_contour(colour="black") + 
  theme_bw()  + coord_cartesian(xlim=c(-25, 25), ylim=c(-20,10)) + ggtitle("true density")
# sample density
gB <- ggplot(samples, aes(x=Q1, y=Q2)) + geom_density2d(colour="black") + 
  theme_bw()  + coord_cartesian(xlim=c(-25, 25), ylim=c(-20,10)) +
  geom_path() + ggtitle("sample density and path")
multiplot(gA,gB,cols=2)
dev.off()

