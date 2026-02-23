###################################################################
## MGWCR: Multiscale Geographically Weighted Cox Regression
## Simulation Scenario 2 - Varied Spatial Heterogeneity (Constant, Medium, High)
###################################################################

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

# Source core MGWCR R-functions and C++ scripts for high-performance computing
source("~/bw.sel.r")
source("~/gw.weight.r")
source("~/gw.dist.r")
source("~/mgwcr_func.R")
sourceCpp("~/cox_derivatives_parallel.cpp")
sourceCpp("~/gw_reg.cpp")
sourceCpp("~/MGWmodel.cpp")


###################################################################
### Data generation function
### Generates survival data with three covariates 
###################################################################

data.generator1 <- function(N, lam, probC, tau){
  getdata.f <- function(id, tau, lam, X1, X2, X3, coor) {
    
    # Define true beta surfaces with varying spatial scales
    tbeta1 <- 1
    tbeta2 <- 1/24 * (coor[1]+coor[2])
    tbeta3 <- 1/360 * ((36-(6-coor[1]/2)^2)*(36-(6-coor[2]/2)^2))

    # Hazard function following Cox Proportional Hazards model
    lam <- lam * exp(tbeta1 * X1 + tbeta2 * X2 + tbeta3 * X3)
    
    # Generate survival time using exponential distribution
    cur.t <- rexp(1, rate=lam)
    
    # Handle right-censoring at time 'tau'
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
                      tau = tau, X1 = X1, X2 = X2, X3 = X3,
                      beta1=tbeta1, beta2=tbeta2, beta3=tbeta3)
    return(tmp)
  }
  
  # Censoring time logic
  if (probC == 0) {
    CC <- rep(tau, N)
  } else {
    CC <- rexp(N, rate = ((-1) * log(1 - probC)))
    CC <- ifelse(CC > tau, tau, CC)
  }
  
  # Multivariate normal covariates with specific correlation structure
  mean_vec <- c(0, 0, 0)
  cov_matrix <- matrix(c(1, 0.8, 0.4,
                         0.8, 1, 0.6,
                         0.4, 0.6, 1), nrow = 3, ncol = 3)
  # Generate multivariate normal covariates
  covariates <- mvrnorm(n = N, mu = mean_vec, Sigma = cov_matrix)
  X1 <- covariates[, 1]
  X2 <- covariates[, 2]
  X3 <- covariates[, 3]
  
  # Coordinates representing the locations of a dataset
  griddf <- cbind(as.matrix(expand.grid(latcoords = seq(1, 25, 1),
                                        lngcoords = seq(1, 25, 1))), loc = 1:625)
  griddf <- griddf[rep(seq_len(nrow(griddf)), each = 2), ]
  idx <- sample(1:N,N,replace = FALSE)
  coor <- griddf[idx,]
  
  
  event <- lapply(1:N, function(i) getdata.f(id = i, coor = coor[i,], X1 = X1[i], X2 = X2[i], X3 = X3[i],
                                             tau = CC[i], lam = lam))
  data <- do.call(rbind, event)
  
  return(data)
}

# Main simulation results containers
results <- data.frame()
betahat <- list()
truebetas <- list()
res.se <-  list()

# Run simulation loop (1000 iterations)
iter=1000
for (r in 1:iter){
  print(r)
  set.seed(r)
  data1 <- data.generator1(N = 1250,
                           lam = 1,
                           probC = 0,
                           tau = 10)

###################################################################### 
# Initial Settings for Calibration
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
max.iterations=30
threshold1=10^-6
threshold2=10^-5
  
kernel.id <- switch(kernel,
                      gaussian = 0,
                      exponential = 1,
                      bisquare = 2,
                      tricube = 3,
                      boxcar = 4)
  
###################################################################### 
## Data Pre-processing
######################################################################

if (is(data, "Spatial")){
  p4s <- proj4string(data)
  dp.locat<-coordinates(data)
  regression.points <- data
  data <- as(data, "data.frame")
}else{
  stop("Given regression data must be Spatial*DataFrame")
}
  
dp.n <- nrow(dp.locat) # sample size
x <- data[,varnames[-(1:2)]]
time <- data[,varnames[1]]
status <- as.integer(data[,varnames[2]])
var.n <- ncol(x)  
  
  
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
InDevars <- colnames(x)  # x names


# Distance matrix calculation
dMats <- list()
dMat <- gw.dist(dp.locat,longlat=longlat)
dMats[[1]] <- dMat
var.dMat.indx <- rep(1, var.n)


###################################################################### 
## Initialization with GW Cox estimates 
######################################################################

load("~/GWCR_betas_mul.RData")
load("~/GWCR_bws_mul.RData")
betas0 <- betas <- as.matrix(GWCR_betas_fix[[r]][,1:var.n])
bws0 <- rep(GWCR_bws_fix[r,2], var.n)

# Construct initial adjusted dependent variable Z
eta.i <- rowSums(betas*x1)
lpl <- cox_derivatives_cpp(eta.i, time, status, num_threads = 5)
dd <- -lpl$second_derivative
d.i <- diag(dd)
u.i <- lpl$first_derivative
z.i <- eta.i + solve(d.i)%*%u.i 
gc()

# Initialize additive terms f_j
f.i <- betas*x1 
f.i2 <-rowSums(f.i)

# Hatmatrix initialization
if(hatmatrix){
  Shat <- matrix(nrow=dp.n, ncol=dp.n)
  S.arrays <- array(dim=c(var.n, dp.n, dp.n))
  C <- array(dim=c(dp.n,var.n,dp.n))
  Beta_SE <- matrix(nrow=dp.n, ncol=var.n)
}


#######################################################
## Backfitting algorithm 
#######################################################
iteration2 <- 0 
criteria.outer <- 10000000

while((iteration2 < max.iterations) && criteria.outer > threshold2){ 
  print(paste0("iteration2: ", iteration2," / ","criteria.outer: ",criteria.outer))
  eta.i.se <- eta.i
  f.i_old <- f.i
  f.i <- betas0 * x1
  f.i2 <-rowSums(f.i)
  socf_nu <- socf_de <- 0
  
  # Preliminary GWR of Z on X
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
  
  # Inner loop: Iterate over each covariate
  while((iteration1 < max.iterations) && criteria.inner > threshold1){
    for(i in 1:var.n){
      dMat <- dMats[[var.dMat.indx[i]]]  
      z.res <- z.i - f.i2 + f.i[,i] 
      # Bandwidth selection
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
      
      # Update coefficients and hat matrices
      res <- gwr.q2(matrix(x1[,i], ncol=1), z.res, d.i, dp.locat, adaptive=adaptive, hatmatrix = hatmatrix, bw=bw.i, kernel=kernel.id,dMat=dMat)
      betai <- res[[1]]
      betas[,i] <- betai
      f.i[,i] <- betas[,i]*x1[,i]
      f.i2 <- rowSums(f.i)
      Si <- res[[2]] # Aj
      S.arrayi <- S.arrays[i,,] # Rj_old
      S.arrays[i,,] <- Si%*%S.arrayi + Si - Si%*%Shat # Rj
      Shat <- Shat- S.arrayi + S.arrays[i,,] # S 
    }
    bws.vars <- rbind(bws.vars, bws) 
    
    # Check convergence for inner loop via log-likelihood
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
  # Recalculate Z based on updated linear predictor
  eta.i <- rowSums(betas*x1)
  lpl <- cox_derivatives_cpp(eta.i, time, status, num_threads = 5)
  dd <- -lpl$second_derivative
  d.i <- diag(dd)
  u.i <- lpl$first_derivative
  z.i <- eta.i + solve(d.i)%*%u.i
  gc()
  
  # Outer loop convergence check
  socf_nu <- sum(colMeans((f.i-f.i_old)^2))
  socf_de <- sum(rowSums((f.i)^2))
  socf <- sqrt(socf_nu/socf_de) 
  criteria.outer <- socf
  
  iteration2 <- iteration2+1 
}

# Store iteration results
last_row <- bws.vars[nrow(bws.vars), ]
results <- rbind(results, c(r, last_row))
betas_df <- as.data.frame(betas) %>% mutate(loc=data1$loc)
betahat[r] <- list(betas_df)
truebetas[r] <- list(data1[,c("beta1","beta2","beta3")])

## Final Standard Error estimation
lpl2 <- cox_derivatives_cpp(eta.i.se, time, status, num_threads = 5)
dd <- -lpl2$second_derivative
S_arrays_Cpp <- aperm(S.arrays, c(2, 3, 1)) 
beta.se <- calculateBetaSE_Cpp(as.matrix(x1), dd, S_arrays_Cpp)
res.se[r] <- list(cbind(beta.se,data1$loc)) 

}  
colnames(results) <- c("r", paste0("mgwcr", 1:(ncol(results) - 1)))

MGWCR_bws_fix <- results
MGWCR_betas_fix <- betahat
MGWCR_trues <- truebetas
MGWCR_se <- res.se




###################################################################
### Post-simulation Performance Analysis
###################################################################
# standard error
betahat_with_r <- lapply(1:length(MGWCR_betas_fix), function(i) {
  df <- as.data.frame(MGWCR_betas_fix[[i]])
  df$r <- i
  return(df)
})
betahat_df.m <- do.call(rbind, betahat_with_r)%>% group_by(loc, r) %>% slice(1)
betas_sd.m <- betahat_df.m %>% group_by(loc) %>%
  summarise(V1 = sd(X1),
            V2 = sd(X2),
            V3 = sd(X3))
MESE <- colMeans(betas_sd.m);MESE

res.se_df <- do.call(rbind, MGWCR_se)
se_avg <- as.data.frame(res.se_df)%>% group_by(V4) %>%
  summarise(V1 = mean(V1),
            V2 = mean(V2),
            V3 = mean(V3))
MASE <- colMeans(se_avg);MASE

# print(betas_sd.m, n = Inf)
# print(se_avg, n = Inf)


# performance
process <- function(parMat, trueBetas, standarderror){
  MAB <- mean(abs(parMat - trueBetas))
  MESE <-sd(parMat)
  MSE <- mean((parMat - trueBetas)^2)
  MASE <- mean(standarderror)
  return(c(MAB = MAB,  MSE = MSE, MESE=MESE, MASE=MASE))
}

results.m <- lapply(1:3, function(i) {
  sapply(1:iter, function(r) {
    process(MGWCR_betas_fix[[r]][,i], MGWCR_trues[[r]][,i], MGWCR_se[[r]][,i])
  })
})
mean_results <- lapply(results.m, rowMeans)
for (i in 1:3) print(mean_results[[i]])

# estimates(betas)
betas_mgwcr <- do.call(rbind, MGWCR_betas_fix)
betas_mgwcr2 <- as.data.frame(betas_mgwcr)%>% group_by(loc) %>%
  summarise(V1 = mean(X1),
            V2 = mean(X2),
            V3 = mean(X3))
print(betas_mgwcr2, n = Inf)

# mab, mse
mse_mgwcr2<-cbind(t(results.m[[1]][1:2,]),t(results.m[[2]][1:2,]),t(results.m[[3]][1:2,]))
