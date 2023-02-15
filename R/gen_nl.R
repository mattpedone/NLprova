genoutcome <- function(nobs, alpha, beta1, beta2, beta3, metx, x2, x3){
  prob <- matrix(0, nrow = nobs, ncol = 3)
  myy <- rep(0, nobs)
  
  for(i in 1:nobs){
    temp <- NULL
    temp <- beta1*(metx[i])^2 + beta2*x2[i] + beta3*x3[i]
    eta1 <- alpha[1] + temp[1]
    eta2 <- alpha[2] + temp[2]
    
    prob[i,1] <- 1/((exp(eta1)+1)*(exp(eta2)+1))
    prob[i,2] <- exp(eta1)/((exp(eta1)+1)*(exp(eta2)+1))
    prob[i,3] <- exp(eta2)/(exp(eta2)+1)
    myy[i] <- sample(c(0,1,2), size = 1, prob = prob[i,])
  }
  outcome <- round(cbind(myy,prob), 4)
  return(outcome)
}

getdata <- function(...)
{
  e <- new.env()
  name <- data(..., envir = e)[1]
  e[[name]]
}

pred_sample <- function(p = 25, o = .80){
  ov <- floor(p*o)
  all <- c(1:p)
  comm <- sample(x = all, size = ov, replace = F)
  all <- all[-comm]
  res <- length(all)
  all <- sample(all, res)
  
  g1 <- all[1:floor(res/2)]
  g2 <- all[(floor(res/2)+1):res]
  g1 <- c(g1, comm)
  g2 <- c(g2, comm)
  return(list(g1 = sort(g1), g2=sort(g2)))
}

tt <- function(pred, prog){
  genenorm <- pred
  mypca <- prcomp(genenorm)
  
  nobs <- length(mypca$x[,1])
  
  # Predictive Markers
  ## use the combination of the pca to generate the exponential
  #metx1 <- scale(mypca$x[,1]) + scale(mypca$x[,2]) + 0.50
  #metx <- sign(metx1)*(sign(metx1)*metx1)^(0.2) * 0.45
  x1 <- scale(mypca$x[,1])
  x2 <- scale(mypca$x[,2])
  metx <- cos(x1) + cos(x2)
  #rgl::plot3d(x1, x2, metx, col = "blue", size = 3)
  # Prognostic Markers
  ## transformation
  z2 <- prog[,1]
  z3 <- prog[,2]
  z2 <- sign(z2)*(sign(z2)*z2)^(0.5)
  z3 <- sign(z3)*(sign(z3)*z3)^(0.2)
  
  # pmts probabilities for treatment 1
  alpha1 <- c(-0.5, -1)
  beta11 <- c(1.5, 2)
  #beta11 <- c(1, 1)
  # pmts probabilities for treatment 2
  alpha2 <- c(0.7, -1)
  beta21 <- c(.5, 1.0)
  #beta21 <- c(1, 1)
  # pmts probabilities with prognostic only
  alpha3 <- c(1, -0.5); beta2 <- c(1, 0.5); beta3 <- c(0.7,1)
  
  #probabilities for treatment 1
  prob1 <- genoutcome(nobs, alpha1, beta11 ,c(0,0), c(0,0), metx, z2, z3)
  #probabilities for treatment 2
  prob2 <- genoutcome(nobs, alpha2, beta21, c(0,0), c(0,0), metx, z2, z3)
  probprog <- genoutcome(nobs, alpha3, c(0,0), beta2, beta3, metx, z2, z3)
  prog <- cbind(z2, z3)
  
  # Now we construct prob with both prog and pred features
  myprob1 <- myprob2 <- matrix( 0, nrow = nobs, ncol = 3)
  
  for (i in 1:nobs){
    myprob1[i,] <- prob1[i,2:4]*probprog[i,2:4]/sum(prob1[i,2:4]*probprog[i,2:4])
    myprob2[i,] <- prob2[i,2:4]*probprog[i,2:4]/sum(prob2[i,2:4]*probprog[i,2:4])
  }
  
  myprob <- list(myprob1, myprob2)
  trtsgn <- rep(c(1,2), nobs/2)
  
  myy <- matrix(0, nrow = nobs, ncol = 3)
  myyout <- matrix(0, nrow = nobs, ncol = 1)
  for(k in 1:nobs){
    trtemp <- trtsgn[k];
    if(trtemp == 1){
      myy[k,1:3] <- t(rmultinom(n = 1, size = 1, prob = myprob[[1]][k,]))
    }
    if(trtemp == 2){
      myy[k,1:3] <- t(rmultinom(n = 1, size = 1, prob = myprob[[2]][k,]))
    }
    myyout[k] <- match(1,myy[k,])
    trtemp <- NULL
  }
  
  return(list(myoutot = myyout, mytot = myy, pred = pred, prog = prog,
              trtsgn = trtsgn, myprob = myprob))
}

genmech_het_nl <- function(npred = 10, nset = 30, overlap = 0.8){
  set.seed(121)
  mydata <- getdata("simupats")
  
  #for(i in 1:nset){
  genenorm <- scale(as.matrix(mydata))
  #genenorm <- genenorm[sample(nrow(genenorm)),]
  if(npred > 90) stop("Using the simupats dataset the maximum number of predictive covariates is 90.")
  pred <- genenorm[,c(1:npred)]#restituisco questi, ma riordinati
  prog <- genenorm[,c(91:92)]#restituisco questi, ma riordinati
  
  groups <- pred_sample(p = npred, o = overlap)
  id_train <- c(1:124)
  id_test <- c(125:152)
  
  yord <- ymat <- predmk <- progmk <- trtsgn <- prob <- vector("list", length = nset)
  
  for(i in 1:nset){
    train <- tt(pred = pred[id_train,groups$g1], prog = prog[id_train,])
    test <- tt(pred[id_test,groups$g2], prog[id_test,])
    
    yord[[i]] <- rbind(train$myoutot, test$myoutot)
    ymat[[i]] <- rbind(train$mytot, test$mytot)
    predmk[[i]] <- pred[c(id_train, id_test),]
    progmk[[i]] <- prog[c(id_train, id_test),]
    trtsgn[[i]] <- c(train$trtsgn, test$trtsgn)#[c(id_train[,i], id_test[,i])]
    prob[[i]] <- list(rbind(train$myprob[[1]], test$myprob[[1]]),
                      rbind(train$myprob[[2]], test$myprob[[2]]))
  }
  
  return(list(yord = yord, ymat = ymat, pred = predmk, prog = progmk, trtsgn = trtsgn,
              prob = prob))
}