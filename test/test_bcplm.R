if (FALSE){

# load package
library(tweedie)
library(rbenchmark)
library(ggplot2)
library(coda)
library(lme4)
library(amer)
library(cplm)


options(error = recover)
#setwd("~/2011/cplm")
setwd("C:\\Documents and Settings\\CAB2007\\My Documents\\2012\\cplm")
load("./data/fineroot.RData")
source("./R/classMethods.R")
source("./R/cpglm.R")
source("./R/utilities.R")


#dyn.load("./src/cplm.so")
link="log"
control=list()
trace=T
n.chains=3
n.iter = 200
n.burnin = 0
n.sims=600
n.thin=1
n.report=2
inits=NULL
bound.p=c(1.2,1.9)
phi.shape=0.001
phi.scale=0.001
tune.iter=3000
n.tune= 15 
tune.weight= 0.25
bound.phi =100
prior.beta.var=NULL
prior.beta.mean=NULL
contrasts = NULL
block.update.beta = FALSE

dyn.load("src/bcplm.dll")  
.Call("init")

bcplm <- function(formula, link = "log", data, inits = NULL,
                  weights, offset, subset, na.action, contrasts = NULL, 
                  n.chains = 3, n.iter = 2000, n.burnin = floor(n.iter / 2),
                  n.thin = max(1, floor(n.chains * (n.iter - n.burnin) / n.sims)),
                  n.sims = 1000, n.report = 1000, prior.beta.mean = NULL, 
                  prior.beta.var = NULL, bound.phi = 100, bound.p = c(1.01, 1.99), 
                  tune.iter = 4000, n.tune = 10, tune.weight = 0.25, 
                  basisGenerators = c("tp","tpU","bsp","sp2d"), 
                  block.update.beta = FALSE, ...) {} 

set.seed(10)

# test glm 
call <- amer:::expand.call(bcplm,call("bcplm", RLD~ Stock * Zone ,
                               data= FineRoot))

formula <- RLD ~ factor(Zone) * factor(Stock)


# test glmm
call <- amer:::expand.call(bcplm,call("bcplm",RLD ~  Stock * Zone +  (1|Plant) , 
                                       data = FineRoot))
formula <- RLD ~  Stock * Zone +  (1|Plant)


fineroot$a <- rnorm(nrow(fineroot))
formula <- RLD~ Zone+(1+a|Plant)+(1|Stock) 
contrasts = NULL
call <- amer:::expand.call(bcplm,call("bcplm", RLD~ Zone+(1+a|Plant)+(1|Stock),data=fineroot))


dyn.unload("src/bcplm.dll")             

# fit the fineroot data with Bayesian models
# first use the latent variable approach
set.seed(10)
fit1 <- bcplm(RLD~ factor(Zone)*factor(Stock), data=fineroot,
               n.iter=20000, n.burnin=10000,
               n.thin=20, n.report=5000, method = "gibbs")
gelman.diag(fit1$sims.list)
# diagnostic plots                             
acfplot(fit1$sims.list,lag.max=20)
xyplot(fit1$sims.list)                              
densityplot(fit1$sims.list)               

# summary results
summary(fit1$sims.list)                             

# now try direct density evaluation (This is much slower)
set.seed(10)
fit2 <- bcpglm(RLD~ factor(Zone)*factor(Stock), data=fineroot,
               n.iter=11000, n.burnin=1000,
               n.thin=10, n.report=5000)
gelman.diag(fit2$sims.list)
summary(fit2)                             


# now fit the Bayesian model to an insurance loss triangle 
# (see Peters et al. 2009)
fit3 <- bcplm(increLoss~ factor(year)+factor(lag), data=InsReserve,
               tune.iter=5000, n.tune = 10,
               n.iter=11000, n.burnin=1000,
               n.thin=10, n.report=5000, bound.p=c(1.1,1.95))
gelman.diag(fit3$sims.list)
summary(fit3$sims.list)                             


fit0<- cpglmm(RLD ~  Stock +  (1|Plant), data=fineroot)


set.seed(10)
fit1 <- bcplm(RLD ~  Stock +   (1|Plant), data=fineroot,
              n.iter=60000, n.burnin=10000,
              n.thin=50, n.report=5000, method = "gibbs")
gelman.diag(fit1$sims.list)
# diagnostic plots                             
acfplot(fit1$sims.list,lag.max=20)
xyplot(fit1$sims.list)                              
densityplot(fit1$sims.list)               

# summary results
summary(fit1$sims.list)


set.seed(12)
fit2<- bcplm(RLD ~  Stock +   (1|Plant), data=fineroot,
              n.iter=20000, n.burnin=10000,
              n.thin=20, n.report=5000, method = "gibbs")
summary(fit2$sims.list)
xyplot(fit2$sims.list)



mu1 <- as.numeric(exp(input@X %*% input@fixef + t(input@Zt) %*% input@u))

.Call("bcplm_update_mu", input)

sum(mu1 - input@mu)
#####################
# check the computation of full conditionals

# 1. post of beta_k
post_betak <- function(x, da){
  beta <- da@fixef 
  p2 <- 2 - da@p 
  p1 <- da@p - 1
  k <- da@k + 1
  beta[k] <- x
  mu <- as.numeric(exp(da@X %*% beta + t(da@Zt) %*% da@u))
  -(sum(mu^p2)/p2 + sum(da@y * mu^(-p1))/p1)/da@phi - 
    0.5 * (x - da@pbeta.mean[k])^2/da@pbeta.var[k]  
}

input@k <- as.integer(1)
post_betak(1.0, input)
.Call("bcplm_update_mu", input)
.Call("bcplm_post_betak", 1.0, input)


# 2. post of phi:
post_phi <- function(x, da){
  sum(log(dtweedie(da@y, mu = da@mu, power = da@p, phi = x)))  
}


post_phi(1.3, input)
.Call("bcplm_post_phi", 1.3, input)


# 2. post of p:
post_p <- function(x, da){
  sum(log(dtweedie(da@y, mu = da@mu, power = x, phi = da@phi)))  
}

post_p(1.3, input)
.Call("bcplm_post_p", 1.3, input)


# 4. post of u_k: assume "bcplm_update_mu" is called
post_uk <- function(x, da){
  u <- da@u 
  p2 <- 2 - da@p 
  p1 <- da@p - 1
  k <- da@k + 1
  u[k] <- x
  mu <- as.numeric(exp(da@X %*% da@fixef + t(da@Zt) %*% u))
  -(sum(mu^p2)/p2 + sum(da@y * mu^(-p1))/p1)/da@phi - 
    0.5 * x^2/as.numeric(da@Sigma[[1]]) 
}

# uniform prior
post_sigma <- function(x, da){
  dm <- da@dims
  as.numeric(-dm["n.u"]/2 * log(x)  - sum(da@u^2)/(2 * x^2))
}

post_sigma(0.08, input)

xx = seq(0.01, 0.2, by = 0.01)

yy <- sapply(xx, function(t) post_sigma(t, input))
yy2 <- yy - max(yy)
plot(xx, exp(yy2), type = "l")


input@k <- as.integer(1)
post_uk(1.0, input)
.Call("bcplm_update_mu", input)
.Call("bcplm_post_uk", 1.0, input)




post_phi(1.0, input)
.Call("bcplm_post_phi", 1.0, input)

pv <- c(0.01802352, 0.05255602, 0.04535603, 0.05694167, 0.11490249, 0.12897299,
0.002066605, 0.0006284337,
0.03479402, 0.0292, 0.0310, 0.0880, 0.02506620, 0.0257, 0.0347, 0.0387,
0.5)

input@mh.var <- pv

fit <- cpglmm(RLD ~ Stock * Zone + (1|Plant), data = FineRoot)
input@u <- fit@ranef
input@Sigma <- list(matrix(as.vector(VarCorr(fit)[[1]])))

do_mcmc <- function(input){
  psd <- sqrt(input@mh.var)
  ig.shape <- 0.0001
  ig.scale <- 0.0001
  dm <- input@dims
  acc <- rep(0, dm["n.beta"] + dm["n.u"] + 2 + 1)
  # update xb, zu, eta and mu slots etc
  .Call("bcplm_update_mu", input)
  
  sim <- matrix(0, dm["n.keep"], dm["n.beta"] + dm["n.u"] + 3)
  ns <- 0
  for (i in 1:dm["n.iter"]){
    
    # sample beta
    for (k in 1:dm["n.beta"]){
      input@k <- as.integer(k - 1)
      bk <- metrop_rw(1, input@fixef[k], psd[k], post_betak, input)
      input@fixef[k] <- as.vector(bk)
      acc[k] <- acc[k] + attr(bk, "accept")
    }
        
    # sample phi and p
    input@mu <- as.numeric(exp(input@X %*% input@fixef + t(input@Zt) %*% input@u)) 
    pos <- dm["n.beta"] + 1
    phi <- metrop_rw(1, input@phi, psd[pos], post_phi, input, lower = 0, upper = 100)
    input@phi <- phi
    acc[pos] <- acc[pos] + attr(phi, "accept")
    p <- metrop_rw(1, input@p, psd[pos + 1], post_p, input, lower = 1.01, upper = 1.99)
    input@p <- p
    acc[pos + 1] <- acc[pos + 1] + attr(p, "accept")
    
    # sample u
    for (k in 1:dm["n.u"]){
      input@k <- as.integer(k - 1)
      pos <- dm["n.beta"] + 2 + k
      uk <- metrop_rw(1, input@u[k], psd[pos], post_uk, input)
      input@u[k] <- as.vector(uk)
      acc[pos] <- acc[pos] + attr(uk, "accept")
      
    }
    
    # simulate Sigma
     #input@Sigma[[1]] <- matrix(1/rgamma(1, dm["n.u"]/2 + ig.shape, sum(input@u^2)/2 + ig.scale))
    pos <- dm["n.beta"] + 2 + dm["n.u"] + 1
    sigma <- metrop_rw(1, as.numeric(sqrt(input@Sigma[[1]])), 
                       psd[pos], post_sigma, input, lower = 0, upper = 100)
    
    input@Sigma <- list(matrix(sigma^2))
    acc[pos] <- acc[pos] + attr(sigma, "accept")
    
    
    # set parameters
    if (i > dm["n.burnin"] && ((i - dm["n.burnin"])%%dm["n.thin"] == 0)){
      ns <- ns + 1
      sim[ns, ] <- c(input@fixef, input@phi, input@p, input@u, as.numeric(input@Sigma[[1]]))
    }
  }
  
  return(list(accept = acc, sim = sim))
  
}

dm <- input@dims
dm["n.iter"] <- as.integer(500)
dm["n.burnin"] <- as.integer(0)
dm["n.thin"] <- as.integer(1)
dm["n.keep"] <- as.integer(500)
input@dims <- dm



dm <- input@dims
dm["n.iter"] <- as.integer(6000)
dm["n.burnin"] <- as.integer(1000)
dm["n.thin"] <- as.integer(5)
dm["n.keep"] <- as.integer(3000)
input@dims <- dm

a <- do_mcmc(input)
a[[1]]/500


out <- vector("list", dm["n.chains"])
for (i in 1:dm["n.chains"]){
  init <- input@inits[[i]]
  
  input@fixef <- init[1:dm["n.beta"]]
  input@phi <- init[dm["n.beta"] + 1]
  input@p <- init[dm["n.beta"] + 2]
  input@u <- init[(dm["n.beta"] + 3):(dm["n.beta"] + dm["n.u"] + 2)]
  input@Sigma[[1]] <- matrix(init[length(init)])
  .Call("bcplm_update_mu", input)
  cat("chain:")
  cat(i)
  out[[i]] <- do_mcmc(input)
  
}

a = lapply(out, "[[", 2)
b= lapply(a, as.mcmc)
b <- as.mcmc.list(b)
summary(b)



set.seed(10)
ss1 <- mh(500, postLog, xinit=5,lb=0,sigma=0.1)
ss2 <- mh(500, postLog, xinit=5,lb=0,sigma=50)

fun <- function(x, sd2 = 5){
  dnorm(x, 0, sd = sd2, log = TRUE)
}
n <- 100000
set.seed(10)
s1 <- metrop_rw(n, par = 0, sd = 7, fun, sd2 = 5)
set.seed(10)
s2 <- metrop_rw1(n, par = 0, v = 50, fun, sd = 5)

c(mean(s1), sd(s1))
c(mean(s2), sd(s2))

attr(s1, "accept")/n
for (i in 1:20){
m1 <- MCMCmetrop1R(fun, 0, burnin = 0, seed =i,mcmc = n, V = matrix(50), sd = 5)
print(c(mean(m1), sd(m1)))
}


library(rbenchmark)
f1 <- function() metrop_rw1(n, par = 0, v = 50, fun, sd = 5)
f2 <- function() MCMCmetrop1R(fun, 0, burnin = 0, mcmc = n, 
                              V = matrix(50), sd = 5, verbose = FALSE)

benchmark(f1(), f2(), replications = 50)

metrop_rw1(10, 0, 0.1, post_uk, input)

dyn.load("src/bcplm.dll")  
.Call("init")


dyn.unload("src/bcplm.dll")  






}
