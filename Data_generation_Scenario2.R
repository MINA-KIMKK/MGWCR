############################################################
############################################################
### Simulation design 2
### Equal Degrees of Spatial Heterogeneity
### Data generation
### Mina Kim
############################################################
############################################################


data.generator2 <- function(N, lam, probC, tau){
  getdata.f <- function(id, tau, lam, X1, X2, coor) {
    
    tbeta1 <- 1/24 * (coor[1]+coor[2])
    tbeta2 <- 1/24 * (coor[1]-coor[2])
    
    lam <- lam * exp(tbeta1 * X1 + tbeta2 * X2)
    
    cur.t <- rexp(1, rate=lam)
    if (cur.t >= tau) {
      estart <- 0
      estop <- tau
      estatus <- 0
    } else{
      estart <- 0
      estop <- cur.t
      estatus <- 1
    }
    
    tmp <- data.frame(id = id, loc1 = coor[1], loc2 = coor[2], loc = coor[3],
                      estart = estart, estop = estop, estatus = estatus,
                      tau = tau, X1 = X1, X2 = X2, beta1=tbeta1, beta2=tbeta2)
    return(tmp)
  }
  
  
  if (probC == 0) {
    CC <- rep(tau, N)
  } else {
    CC <- rexp(N, rate = ((-1) * log(1 - probC)))
    CC <- ifelse(CC > tau, tau, CC)
  }
  
  # covariates
  X1 <- rnorm(N, 0, 1)
  X2 <- rnorm(N, 0, 1)
  
  
  # coordinates representing the locations of a dataset
  griddf <- cbind(as.matrix(expand.grid(latcoords = seq(1, 25, 1),
                                        lngcoords = seq(1, 25, 1))), loc = 1:625)
  griddf <- griddf[rep(seq_len(nrow(griddf)), each = 2), ]
  idx <- sample(1:N,N,replace = FALSE)
  coor <- griddf[idx,]
  
  
  event <- lapply(1:N, function(i) getdata.f(id = i, coor = coor[i,], X1 = X1[i], X2 = X2[i],
                                             tau = CC[i], lam = lam))
  data <- do.call(rbind, event)
  
  return(data)
}
