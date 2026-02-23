###################################################################
## GW Cox: Geographically Weighted Cox Regression
## Simulation Scenario 1 - Same Variation Scale 
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
### Generates survival data with two covariates 
###################################################################

data.generator2 <- function(N, lam, probC, tau){
  getdata.f <- function(id, tau, lam, X1, X2, coor) {
    
    # Define true beta surfaces with uniform spatial scale
    tbeta1 <- 1/24 * (coor[1]+coor[2])
    tbeta2 <- 1/24 * (coor[1]-coor[2])
    
    # Hazard function following Cox Proportional Hazards model
    lam <- lam * exp(tbeta1 * X1 + tbeta2 * X2)
    
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
                      tau = tau, X1 = X1, X2 = X2, beta1=tbeta1, beta2=tbeta2)
    return(tmp)
  }
  
  # Censoring time logic
  if (probC == 0) {
    CC <- rep(tau, N)
  } else {
    CC <- rexp(N, rate = ((-1) * log(1 - probC)))
    CC <- ifelse(CC > tau, tau, CC)
  }
  
  # Multivariate normal covariates
  mean_vec <- c(0, 0)  
  cov_matrix <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  # Generate multivariate normal covariates
  covariates <- mvrnorm(n = N, mu = mean_vec, Sigma = cov_matrix)
  X1 <- covariates[, 1]
  X2 <- covariates[, 2]
  
  # Coordinates representing the locations of a dataset
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

# Main simulation results containers
results <- data.frame()
betahat <- list()
truebetas <- list()
res.se <-  list()

# Run simulation loop (1000 iterations)
iter=1000
for (r in 1:1000){
  print(r)
  set.seed(r)
  data1 <- data.generator2(N = 1250,
                           lam = 1,
                           probC = 0,
                           tau = 10)
kk
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
xx1 <- x 

if (is.null(nrow(xx1))) xx1 <- matrix(xx1, nrow=length(xx1))
if(n.cent>1){
  predictors.centered.means <- colMeans(xx1[,predictor.centered])
} else {
  predictors.centered.means <- mean(xx1[,predictor.centered])
}

if(n.cent>1){
  xx1[,predictor.centered] <- scale(xx1[,predictor.centered], scale=FALSE)
}

colnames(xx1) <-  colnames(x) # x : original, xx1 : centered x


# Distance matrix calculation
dMat <- gw.dist(dp.locat,longlat=longlat)

###################################################################### 
## Initialization 
######################################################################

# Construct initial adjusted dependent variable Z
betas <- matrix(0, nrow=nrow(data1), ncol=var.n)
eta.i <- rowSums(betas*xx1)
lpl <- cox_derivatives_cpp(eta.i, time, status, num_threads = 5)
dd <- -lpl$second_derivative
d.i <- diag(dd)
u.i <- lpl$first_derivative
z.i <- eta.i + solve(d.i)%*%u.i 
gc()

# Initialize additive terms f_j
f.i <- betas*xx1


#######################################################
## Backfitting algorithm 
#######################################################

iteration = 0 
criterion.val <- 10000000  
while ((iteration < max.iterations) && criterion.val > threshold){
  # Update the linear predictor (eta) and linearized dependent variable (Z)
  eta.i <- rowSums(betas*xx1)
  lpl <- cox_derivatives_cpp(eta.i, time, status, num_threads = 5)
  dd <- -lpl$second_derivative
  d.i <- diag(dd)
  u.i <- lpl$first_derivative
  z.i <- eta.i + solve(d.i)%*%u.i 
  gc()
  
  f.i_old <- f.i
  diffi <- 0
  f.i2 <- 0
  #Global spatial bandwidth selection
  bw1<-bw.gwr2(as.matrix(xx1), z.i,  dp.locat, dd, approach = approach, kernel=kernel.id, adaptive = adaptive, dMat, verbose = verbose, nlower = nlower,parallel.method=parallel.method,parallel.arg=parallel.arg)
  print(paste0("iteration: ",iteration,"////", criterion.val," ///  bw:",bw1))
  # Weighted regression to update local coefficients
  res <- gwr.q2(as.matrix(xx1), z.i, d.i, dp.locat, adaptive=adaptive, hatmatrix = hatmatrix, bw=bw1, kernel=kernel.id,dMat=dMat)
  betas <- res[[1]]
  bws.vars <- c(bws.vars, bw1)
  
  f.i <- betas * xx1
  diffi <- sum(colMeans((f.i - f.i_old)^2))
  f.i2 <- rowSums(f.i)
  socf <- sqrt(diffi/sum(f.i2^2)) 
  criterion.val <- socf
  iteration <- iteration+1 
  
}
# Optimal single bandwidth for the GW Cox model
last_row <- bws.vars[length(bws.vars)]
results <- rbind(results, c(r, last_row))
betas_df <- as.data.frame(betas) %>% mutate(loc=data1$loc)
betahat[r] <- list(betas_df)
truebetas[r] <- list(data1[,c("beta1","beta2")])

## Standard Error Estimation
S.arrays <- array(dim=c(var.n, dp.n, dp.n))
for (i in 1:var.n){
  res <- gwr.q2(matrix(xx1[,i], ncol=1), z.i, d.i, dp.locat, adaptive=adaptive, hatmatrix = hatmatrix, bw=last_row, kernel=kernel.id,dMat=dMat)
  S.arrays[i,,] <- res[[2]] # Rj
  
}
S_arrays_Cpp <- aperm(S.arrays, c(2, 3, 1)) 
beta.se <- calculateBetaSE_Cpp(as.matrix(xx1), dd, S_arrays_Cpp)
res.se[r] <- list(cbind(beta.se,data1$loc)) 
}

GWCR_bws_fix <- results
GWCR_betas_fix <- betahat
GWCR_trues <- truebetas
GWCR_se <- res.se
colnames(GWCR_bws_fix) <- c("r", paste0("gwcr", 1:(ncol(GWCR_bws_fix) - 1)))

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
    process(GWCR_betas_fix[[r]][,i], GWCR_trues[[r]][,i], GWCR_se[[r]][,i])
  })
})
mean_results <- lapply(results.m, rowMeans)
for (i in 1:2) print(mean_results[[i]])

# estimates(betas)
betas_gwcr <- do.call(rbind, GWCR_betas_fix)
betas_gwcr2 <- as.data.frame(betas_gwcr)%>% group_by(loc) %>%
  summarise(V1 = mean(X1),
            V2 = mean(X2))
print(betas_gwcr2, n = Inf)

# mab, mse
mse_gwcr1<-cbind(t(results.m[[1]][1:2,]),t(results.m[[2]][1:2,]))
