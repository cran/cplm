if (FALSE){
library(cplm)
library(tweedie)
library(pscl)

fit1 <- cpglm(RLD ~ factor(Stock) ,
  data = fineroot)
fit2 <- cpglm(RLD ~ factor(Zone) ,
  data = fineroot)
setwd("C:\\Documents and Settings\\CAB2007\\My Documents\\2011\\cplm")
#setwd("~/2011/cplm")
source("./R/classMethods.R")
source("./R/utilities.R")

# simulate data for use 
gamma <- c(-1, 0.2)
beta <- coef(fit1)
beta[1] <- beta[1] + 1
p <- fit1$p
phi <- fit1$phi
Z <- model.matrix(fit2)
dd <- exp(Z %*% gamma)
q <- dd / (1 + dd)
n <- nrow(fineroot)
X <- model.matrix(fit1)
s <- rbinom(n, 1, q)
mu<- as.numeric(exp(X %*%beta))
fineroot$y <- rtweedie(n, power = p, phi = phi, mu = mu)
fineroot$y[s == 1] <- 0


link="log"
control=list()
trace=T
contrasts = NULL

zcpglm <- function(formula, link = "log", data, weights, offset, 
                  subset, na.action = NULL, contrasts = NULL, 
                  control = list(), chunksize = 0, ...) {}

call <- match.call(zcpglm,call("zcpglm", y ~ factor(Stock) || factor(Zone),
  data=fineroot))
formula <- y ~ factor(Stock) || factor(Zone)

system.time(f1 <- zcpglm(y ~ factor(Stock) || factor(Zone), data=fineroot, control = list(trace =1)))
coef(f1)



#########################################
# run simulations
#########################################
library(cplm)
library(pscl)

fit1 <- cpglm(RLD ~ factor(Stock) ,
  data = fineroot)
fit2 <- cpglm(RLD ~ factor(Zone) ,
  data = fineroot)

# simulate data for use 
gamma <- c(-1, 0.2)
beta <- coef(fit1)
beta[1] <- beta[1] + 1
p <- fit1$p
phi <- fit1$phi
Z <- model.matrix(fit2)
dd <- exp(Z %*% gamma)
q <- dd / (1 + dd)
n <- nrow(fineroot)
X <- model.matrix(fit1)

nsim <- 100 
sims<- matrix(0, nsim, length(gamma)+ length(beta) + 2)

set.seed(10)
for ( i in 1:nsim){
  s <- runif(n) 
  s[s > q] <- 0
  s[s != 0] <- 1
  fineroot$y <- tweedie::rtweedie(n, power = p, phi = phi, mu = as.numeric(exp(X %*%beta)))
  fineroot$y[s == 1] <- 0
  f1 <- zcpglm(y ~ factor(Stock) || factor(Zone), data=fineroot)
  sims[i,] <- as.numeric(c(unlist(coef(f1)), f1$phi, f1$p))
}

apply(sims,2,median)
c(gamma, beta, phi, p)
}