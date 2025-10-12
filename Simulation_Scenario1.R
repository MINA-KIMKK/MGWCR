############################################################
############################################################
### Simulation design 1
### Different degrees of spatial heterogeneity
### Simulation
### Mina Kim
############################################################
############################################################

library(MASS)
library(foreach)
library(doParallel)
library(cvTools)
library(dplyr)
library(survival)
library(sp)
library(Rcpp)
library(RcppParallel)
select <- dplyr::select
mutate <- dplyr::mutate
rename <- dplyr::rename 
summarise <- dplyr::summarise 

source("~/bw.sel.r")
source("~/gw.weight.r")
source("~/gw.dist.r")
source("~/mgwcr_func.R")
sourceCpp("~/cox_derivatives_parallel4.cpp")
sourceCpp("~/gw_reg.cpp")
sourceCpp("~/MGWmodel.cpp")


############################################################

results <- data.frame()
betahat <- list()
truebetas <- list()
res.se <-  list()

iter=1000
for (r in 1:iter){
  print(r)
  set.seed(r)
  data1 <- data.generator1(N = 1250,
                           lam = 1,
                           probC = 0,
                           tau = 10)
  

  data.ori <- data1
  coordinates <- data.ori[, c("loc1", "loc2")]
  data <- SpatialPointsDataFrame(coords = coordinates, data = data.ori)
  varnames = c("estop","estatus","X1","X2","X3")
  formula = Surv(estop,estatus==1) ~ X1 + X2 + X3
  longlat=FALSE
  hatmatrix=TRUE
  approach = "AIC"
  kernel="gaussian"
  adaptive=FALSE
  verbose=F
  nlower = 10
  parallel.method=F
  parallel.arg=NULL
  bws.reOpts=5
  criterion="socf"
  max.iterations=100
  threshold1=10^-6
  threshold2=10^-5
  
  
  kernel.id <- 0
  kernel.id <- switch(kernel,
                      gaussian = 0,
                      exponential = 1,
                      bisquare = 2,
                      tricube = 3,
                      boxcar = 4)
  
  
  ###################################################################### 
  ## Data
  if (is(data, "Spatial")){
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    regression.points <- data
    data <- as(data, "data.frame")
  }else{
    stop("Given regression data must be Spatial*DataFrame")
  }
  
  dp.n <- nrow(dp.locat) 
  x <- data[,varnames[-(1:2)]]
  time <- data[,varnames[1]]
  status <- as.integer(data[,varnames[2]])
  var.n <- ncol(x)
  
  
  #################################################
  # Centering and scaling predictors
  predictor.centered <- rep(TRUE, length.out=var.n) # Centering and scaling predictors
  n.cent <- length(which(predictor.centered))
  x1 <- x 
  
  if (is.null(nrow(x1))) x1 <- matrix(x1, nrow=length(x1))
  if(n.cent>1){
    predictors.centered.means <- colMeans(x1[,predictor.centered])
  } else {
    predictors.centered.means <- mean(x1[,predictor.centered])
  }
  
  if(n.cent>1){
    x1[,predictor.centered] <- scale(x1[,predictor.centered], scale=FALSE)
  }
  
  colnames(x1) <-  colnames(x) # x : original, x1 : centered x
  allvars <- all.vars(formula)
  InDevars <- colnames(x) 
  
  
  #################################################
  ## Distance matrix
  dMats <- list()
  dMat <- gw.dist(dp.locat,longlat=longlat)
  dMats[[1]] <- dMat
  var.dMat.indx <- rep(1, var.n)
  
  
  #################################################
  ## initial setting
  beta0 <- c(0,0,0)# initial beta
  betas <- matrix(rep(beta0, nrow(dp.locat)), nrow = nrow(dp.locat), byrow = TRUE)
  bws0 <- rep(max(dMat), var.n)
  
  # initial adjusted dependent variable z
  eta.i <- rowSums(betas*x1)
  lpl <- cox_derivatives_cpp(eta.i, time, status, num_threads = 5)
  dd <- -lpl$second_derivative
  d.i <- diag(dd)
  u.i <- lpl$first_derivative
  z.i <- eta.i + solve(d.i)%*%u.i 
  gc()
  
  # Initialize the additive term fj
  f.i <- betas*x1 
  f.i2 <-rowSums(f.i)
  
  
  
  #################################################
  ## Hatmatrix for the whole process
  if(hatmatrix){
    Shat <- matrix(nrow=dp.n, ncol=dp.n)
    S.arrays <- array(dim=c(var.n, dp.n, dp.n))
    C <- array(dim=c(dp.n,var.n,dp.n))
    ## SEs
    Beta_SE <- matrix(nrow=dp.n, ncol=var.n)
  }
  
  
  
  #######################################################
  ## Backfitting algorithm
  iteration2 <- 0 
  criteria.outer <- 10000000
  
  while((iteration2 < max.iterations) && criteria.outer > threshold2){ 
    eta.i.se <- eta.i
    f.i_old <- f.i
    f.i <- matrix(rep(beta0, nrow(dp.locat)), nrow = nrow(dp.locat), byrow = TRUE) * x1
    f.i2 <-rowSums(f.i)
    socf_nu <- socf_de <- 0
    
    res.ori <- gwr.q2(as.matrix(x1), z.i,d.i, dp.locat, adaptive=adaptive, hatmatrix = hatmatrix,bw=bws0[1], kernel=kernel.id,dMat=dMat)
    Shat <- res.ori[[2]]
    C <- res.ori[[3]]
    idm <- diag(var.n)
    for(i in 1:var.n){
      for(j in 1:dp.n){
        S.arrays[i,j,] <- x1[j,i]*(idm[i,]%*%C[j,,])
      }
    }
    iteration1 <- 0
    bws <- bws0
    bws.vars <- bws0
    bw.seled <- rep(F, var.n)
    bws.thresholds <- rep(0.1, var.n)
    bws.change.NO <- numeric(var.n)
    criteria.inner <- 10000000
    loglik.old <- 10000000
    
    while((iteration1 < max.iterations) && criteria.inner > threshold1){
      for(i in 1:var.n){
        dMat <- dMats[[var.dMat.indx[i]]]  
        z.res <- z.i - f.i2 + f.i[,i] 
        if(bw.seled[i]){
          bw.i <- bws[i]
        }else{
          bw.i <- bw.gwr2(matrix(x1[, i], ncol = 1), z.res, dp.locat, dd, approach = approach, kernel=kernel.id, adaptive = adaptive, dMat, verbose = verbose, nlower = nlower,parallel.method=parallel.method,parallel.arg=parallel.arg)
          if (abs(bw.i - bws[i]) > bws.thresholds[i]) {
            bws.change.NO[i] <- 0
          } else {
            bws.change.NO[i] <- bws.change.NO[i] + 1
            if(bws.change.NO[i] < bws.reOpts){
            }else{
              bw.seled[i] <- T
            }
          }
        }
        bws[i] <- bw.i 
        res <- gwr.q2(matrix(x1[,i], ncol=1), z.res, d.i, dp.locat, adaptive=adaptive, hatmatrix = hatmatrix, bw=bw.i, kernel=kernel.id,dMat=dMat)
        betai <- res[[1]]
        betas[,i] <- betai
        f.i[,i] <- betas[,i]*x1[,i]
        f.i2 <- rowSums(f.i)
        Si <- res[[2]] 
        S.arrayi <- S.arrays[i,,] 
        S.arrays[i,,] <- Si%*%S.arrayi + Si - Si%*%Shat 
        Shat <- Shat- S.arrayi + S.arrays[i,,]
      }
      bws.vars <- rbind(bws.vars, bws) 

      # Convergence criteria for inner loop
      betas.df <- as.data.frame(betas) 
      colnames(betas.df) <- c("V1","V2","V3")
      tmp <- cbind(betas.df,x1,time,status) %>%
        mutate(X1beta1 = X1*V1,
               X2beta2 = X2*V2,
               X3beta3 = X3*V3) 
      m1 <- coxph(Surv(time, status)~X1beta1 + X2beta2 + X3beta3,
                  data=tmp, init=c(1,1,1), control=coxph.control(iter.max=0))
      loglik.new <- m1$loglik[2]
      criteria.inner <- abs(1-loglik.old/loglik.new)
      loglik.old <- loglik.new
      iteration1 <- iteration1 + 1
    }
    
    eta.i <- rowSums(betas*x1)
    lpl <- cox_derivatives_cpp(eta.i, time, status, num_threads = 5)
    dd <- -lpl$second_derivative
    d.i <- diag(dd)
    u.i <- lpl$first_derivative
    z.i <- eta.i + solve(d.i)%*%u.i
    gc()
    
    # Convergence criteria for outer loop 
    socf_nu <- sum(colMeans((f.i-f.i_old)^2))
    socf_de <- sum(rowSums((f.i)^2))
    socf <- sqrt(socf_nu/socf_de) 
    criteria.outer <- socf
    
    iteration2 <- iteration2+1 
    print(paste0("iteration2: ", iteration2," / ","criteria.outer: ",criteria.outer))
  }
  
  last_row <- bws.vars[nrow(bws.vars), ]
  results <- rbind(results, c(r, last_row))
  betas_df <- as.data.frame(betas) %>% mutate(loc=data1$loc)
  betahat[r] <- list(betas_df)
  truebetas[r] <- list(data1[,c("beta1","beta2","beta3")])
  
  ## Calculate the SEs
  lpl2 <- cox_derivatives_cpp(eta.i.se, time, status, num_threads = 5)
  dd <- -lpl2$second_derivative
  S_arrays_Cpp <- aperm(S.arrays, c(2, 3, 1)) 
  beta.se <- calculateBetaSE_Cpp(as.matrix(x1), dd, S_arrays_Cpp)
  res.se[r] <- list(cbind(beta.se,data1$loc)) 
  
}