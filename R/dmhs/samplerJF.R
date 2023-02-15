rm(list=ls())
set.seed(121)

library(rstan)
library(parallel)
library(doParallel)
library(doRNG)

#options(mc.cores = parallel::detectCores())
registerDoParallel(cores = 50)
rstan_options(auto_write = TRUE)


source("R/gen_nl.R")
#1a
simdata <- genmech_het_nl(npred = 25, nset = 50, overlap = 1.0)
##1b
#simdata <- genmech_het_nl(npred = 50, nset = 50, overlap = 1.0)
##3a
#simdata <- genmech_het_nl(npred = 25, nset = 50, overlap = 0.8)
##3b
#simdata <- genmech_het_nl(npred = 50, nset = 50, overlap = .8)

npc2 <- function(output, trtsgn, myoutot){
  K <- dim(output)[3]
  n <- dim(output)[1]
  myctut <- array(0, dim = c(3, 3, K))
  myctutSum <- NULL
  for (i in 1:K) {
    mycurdata <- output[, , i]
    mypre <- NULL
    pretrt1 <- apply(mycurdata[, 1:3], 1, which.max)
    pretrt2 <- apply(mycurdata[, 4:6], 1, which.max)
    mypreTall <- cbind(pretrt1, pretrt2)
    for (j in 1:n) {
      mypre[j] <- mypreTall[j, trtsgn[j]]
    }
    sts <- table(factor(mypre, levels = 1:3), factor(myoutot, levels = 1:3))
    mysdls <- as.numeric(rownames(sts))
    str1 <- matrix(0, nrow = 3, ncol = 3)
    str1[mysdls, ] <- sts
    myctut[, , i] <- str1 * diag(3)
    myctutSum[i] <- sum(str1 * diag(3))
  }
  res <- cbind(myctutSum)
  return(res)
}

Yl <- simdata$ymat
Xl <- simdata$pred
Zl <- simdata$prog
trtl <- simdata$trtsgn

R <- length(Yl)
n <- dim(Yl[[1]])[1]
J <- dim(Yl[[1]])[2]
p <- dim(Xl[[1]])[2]
q <- dim(Zl[[1]])[2]

K <- 10#repliche
npat_pred <- 28

Y <- array(0, dim = c(n + npat_pred, J, R))
X <- array(0, dim = c(n + npat_pred, p, R))
Z <- array(0, dim = c(n + npat_pred, q, R))
Xdis <- array(0, dim = c(n + npat_pred, (p + p + q + 1), R))
trt <- matrix(0, n + npat_pred, R)


for(r in 1:R){
  Y[1:152,,r] <- Yl[[r]]
  Y[153:180,,r] <- Y[125:152,,r]
  X[1:152,,r] <- Xl[[r]]
  X[153:180,,r] <- X[125:152,,r]
  Z[1:152,,r] <- Zl[[r]]
  Z[153:180,,r] <- Z[125:152,,r]
  trt[1:124,r] <- trtl[[r]][1:124]-1
  trt[125:152,r] <- 0
  trt[153:180,r] <- 1
  Xdis[,,r] <- cbind(X[,,r], Z[,,r], trt[,r], (trt[,r]*X[,,r]))
}

#params1 <- params2 <- array(npat_pred, J, K)

predAPT_all <- array(0, dim = c(npat_pred, 9, K))
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)

wk <- c(0, 40, 100)

res <- foreach(k = 1:K) %dorng%
  {
    X_train <- Xdis[1:124, ,k]
    #X_train <- X_train[which(trt[1:124,k] == 0),]
    Y_train <- Y[1:124, ,k]
    #Y_train <- Y_train[which(trt[1:124,k] == 0),]

    X_test <- Xdis[125:180,,k]
    Y_test <- Y[125:180,,k]

    ss_data = list(N = dim(X_train)[1], M = dim(X_test)[1], J = dim(Y_train)[2],
                    P = dim(X_train)[2], X = X_train, Xp = X_test, Y = Y_train,
                    sd_prior = 1.0, psi = .25)

    fit <- rstan::stan(file = "R/dmhs/model.stan", data = ss_data, cores = 1, iter = 1000,
                        chains = 1, verbose = T, warmup = 200, seed = 121,
                        control = list(max_treedepth = 15, adapt_delta = 0.995))#995))


    #check_hmc_diagnostics(fit)
    params <- rstan::extract(fit)
    pp <- apply(params$pipred, c(2,3), median)

    out <- cbind(pp1 = pp[1:28,], pp2 = pp[29:56,])
    return(out)
  }

#summ <- rstan::summary(fit)$summary

for(k in 1:K){
  pp1 <- res[[k]][,1:3]
  A1 <- pp1%*%wk
  pp2 <- res[[k]][,4:6]
  A2 <- pp2%*%wk

  predAPT_all[, 1, k] <- A1
  predAPT_all[, 2, k] <- A2
  myt <- as.numeric(A1 < A2) + 1
  predAPT_all[, 3, k] <- myt
  predAPT_all[, 4:9, k] <- cbind(pp1, pp2)

  myprob <- simdata$prob[[k]]
}
mywk1 <- myprob[[1]][125:152,]%*%wk
mywk2 <- myprob[[2]][125:152,]%*%wk
optrt <- as.numeric(mywk1 < mywk2) + 1
utsum <- sum(abs(mywk2 - mywk1))
utdiff <- abs(as.numeric(mywk2 - mywk1))

#MOT
PPMXCT <- c()
for(k in 1:K){
  PPMXCT[k] <-  sum(abs(predAPT_all[, 3, k] - optrt))
}

MOT <- c(round(mean(PPMXCT), 4), round(sd(PPMXCT), 4))

#MTUg
PPMXpp <- c()
for(k in 1:K){
  PPMXpp[k] <- -(2*sum(abs((predAPT_all[, 3, k] - optrt)) * utdiff) - utsum);
}

MTUg <- c(round(mean(PPMXpp/utsum), 4), round(sd(PPMXpp/utsum), 4))

#NPC
PPMXCUT <- c()
for(k in 1:K){
  trtsgn_test <- simdata$trtsgn[[k]][125:152]
  temp <- array(0, dim = c(28, 6, 1))
  temp[,,1] <- predAPT_all[, 4:9,k]
  myoutot <- simdata$yord[[k]][125:152,]
  PPMXCUT[k] <- npc2(temp, trtsgn_test, myoutot)
}
NPC <- c(round(mean(PPMXCUT), 4), round(sd(PPMXCUT), 4))

#results
resDMHS <- rbind(MOT, MTUg, NPC)#, WAIC, lpml)
#resDMHS <- rbind(MOT, MTUg)#, WAIC, lpml)
colnames(resDMHS) <- c("mean", "sd")
#save(resDMHS, file="output/simulation-study/main/scen3b_dmhs_res.RData")
resDMHS