 
## For bandwidth selection
bw.gwr2<-function(x, y, dp.locat,di, approach="CV",kernel="gaussian",adaptive=FALSE,dMat, verbose=F, parallel.method=F,parallel.arg=NULL, nlower = 10){
  dp.n <-  nrow(dp.locat)  # data point 개수
  if(adaptive){ # adaptive
    upper <- dp.n
    lower <- nlower
  } else{ # fixed
    upper<-max(dMat) # maximum of distance
    lower<-upper/5000
  }
  ########################## Now the problem for the golden selection is too computationally heavy
  #Select the bandwidth by golden selection
  bw<-NA
  # ## call for functions, # gold, gwr.bic, gwr.aic: GWmodel package
  if(approach == "bic" || approach == "BIC")
    bw <- gold(gwr.bic, lower, upper, adapt.bw = adaptive, x, y, di, kernel, adaptive, dp.locat, p=2, theta=0, longlat=F, dMat, verbose, parallel.method, parallel.arg)
  else if(approach == "aic" || approach == "AIC" || approach == "AICc")
    bw <- gold(gwr.aic, lower, upper, adapt.bw = adaptive, x, y, di, kernel, adaptive, dp.locat, p=2, theta=0, longlat=F, dMat, verbose, parallel.method, parallel.arg)    
  else 
    bw <- gold(gwr.cv, lower, upper, adapt.bw = adaptive, x, y, di, kernel, adaptive, dp.locat, p=2, theta=0, longlat=F, dMat, verbose, parallel.method, parallel.arg)
  # ## stop cluster
  if (parallel.method == "cluster") {
    if (missing(parallel.arg)) stopCluster(parallel.arg)
  }
  bw
}


## Geographically weighted regression betas

gwr.q2 <- function(x, y, d, loc, adaptive=F, hatmatrix,bw=sqrt(var(loc[,1])+var(loc[,2])),
                   kernel, p, theta, longlat,dMat){
  dp.n <- nrow(loc) # location 개수 
  var.n <- ncol(x)
  betas <- matrix(nrow=dp.n, ncol=var.n)
  S <- matrix(nrow=dp.n, ncol=dp.n)
  C <- array(dim=c(dp.n, var.n, dp.n))
  for (j in 1:dp.n){ 
    dist.vi <- dMat[,j]
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    gw.resi<-gw_reg(x,y,W.i, d, hatmatrix=hatmatrix,j)
    betas[j,]<-gw.resi[[1]]
    S[j,]<-gw.resi[[2]]
    C[j,,]<-gw.resi[[3]]
  }
  colnames(betas) <- colnames(x)
  res <- list(betas, S,C)
  res
}

gwr.q3 <- function(x, a, adaptive=F, bw, kernel, dMat, num_threads) {
  registerDoParallel(num_threads)
  dp.n <- nrow(x)  # 관측치 개수
  se_list <- foreach(i = 1:dp.n, .combine = 'rbind', .packages = c("RcppArmadillo")) %dopar% {
    dist.vi <- dMat[, i] 
    W.i <- gw.weight(dist.vi, bw, kernel, adaptive)
    as.vector(gw_reg_se(x, W.i, a))
  }
  stopImplicitCluster()
  return(se_list)
}