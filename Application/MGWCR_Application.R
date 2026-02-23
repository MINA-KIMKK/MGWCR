###################################################################
## MGWCR: Multiscale Geographically Weighted Cox Regression
## Empirical Study on Depression Risk (New York State)
## Nurses' Health Study (NHS) data
###################################################################

library(MASS)
library(foreach)
library(doParallel)
library(cvTools)
library(dplyr)
library(survival)
library(sp)
library(parallel)
library(Rcpp)
library(RcppParallel)
library(lubridate)
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


#############################################################
### Data Loading and Feature Engineering
#############################################################
# Load final dataset containing coordinates and survival variables
dat.final <- read.csv("~/final_data.csv")
data.ori <- dat.final
coordinates <- data.ori[, c("loc1", "loc2")]

# Define a SpatialPointsDataFrame for geographic modeling
data <- SpatialPointsDataFrame(coords = coordinates, data = data.ori)

# Survival time, status(depression), and five risk factor covariates.
varnames = c("time","status","agemo","marriage","income2","pm25","NDVI270")
formula = Surv(time,status==1) ~ agemo + marriage + income2 + pm25 + NDVI270

# Initial Settings
longlat=TRUE # Great Circle Distance calculation
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
max.iterations1=50
max.iterations2=20
threshold1=10^-6
threshold2=10^-5
num_threads=16
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
predictor.centered <- rep(TRUE, length.out=var.n) 
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

# Construct geographic distance matrix using Great Circle Distance
dMat <- gw.dist(dp.locat,longlat=longlat)


###################################################################### 
## Initialization with GW Cox estimates 
######################################################################

load("~/GWCR_betas_app_ny.RData")
load("~/GWCR_bws_app_ny.RData")
betas0 <- betas <- as.matrix(GWCR_betas_app[,1:var.n])
bws0 <- rep(GWCR_bws_app, var.n)

# Construct initial adjusted dependent variable Z
eta.i <- rowSums(betas*x1)
lpl <- cox_derivatives_cpp(eta.i, time, status, num_threads = num_threads) 
dd <- -lpl$second_derivative
dd[dd == 0] <- 1e-10
d.i <- diag(dd)
u.i <- lpl$first_derivative
z.i <- eta.i + diag(1/dd)%*%u.i
gc()

# Initialize the additive term fj
f.i <- betas*x1
f.i2 <-rowSums(f.i)

# Hatmatrix initialization
if(hatmatrix){
  Shat <- matrix(nrow=dp.n, ncol=dp.n)
  S.arrays <- array(dim=c(var.n, dp.n, dp.n))
  C <- array(dim=c(dp.n,var.n,dp.n))
}



#######################################################
## Backfitting algorithm 
#######################################################

iteration2 <- 0 
criteria.outer <- 10000000

while((iteration2 < max.iterations2) && criteria.outer > threshold2){ 
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
  while((iteration1 < max.iterations1) && criteria.inner > threshold1){
    print(paste0("iter1: ", iteration1," / ","criteria.inner: ",criteria.inner))
    for(i in 1:var.n){
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
      
      result <- update_all_rows_cpp(Si, S.arrayi, Shat, num_threads = num_threads) 
      
      S.arrays[i,,] <- result$temp1 + Si - result$temp2
      Shat <- Shat - S.arrayi + S.arrays[i,,]
      
    }
    print(bws)
    bws.vars <- rbind(bws.vars, bws) 

    # Check convergence for inner loop via log-likelihood
    betas.df <- as.data.frame(betas) 
    
    
    colnames(betas.df) <- c("V1","V2","V3","V4","V5")
    tmp <- cbind(betas.df,x1,time,status) %>%
      mutate(X1beta1 = agemo*V1,
             X2beta2 = marriage*V2,
             X3beta3 = income2*V3,
             X4beta4 = pm25*V4,
             X5beta5 = NDVI270*V5)
    m1 <- coxph(Surv(time, status)~X1beta1 + X2beta2 + X3beta3 + X4beta4 + X5beta5,
                data=tmp, init=c(1,1,1,1,1), control=coxph.control(iter.max=0))
    
    loglik.new <- m1$loglik[2]
    criteria.inner <- abs(1-loglik.old/loglik.new)
    loglik.old <- loglik.new
    iteration1 <- iteration1 + 1
  }
  # Recalculate Z based on updated linear predictor
  eta.i <- rowSums(betas*x1)
  lpl <- cox_derivatives_cpp(eta.i, time, status, num_threads = num_threads)
  dd <- -lpl$second_derivative
  d.i <- diag(dd)
  dd[dd == 0] <- 1e-10
  u.i <- lpl$first_derivative
  z.i <- eta.i + diag(1/dd)%*%u.i
  gc()
  
  # Outer loop convergence check
  socf_nu <- sum(colMeans((f.i-f.i_old)^2))
  socf_de <- sum(rowSums((f.i)^2))
  socf <- sqrt(socf_nu/socf_de) 
  criteria.outer <- socf
  
  iteration2 <- iteration2+1 
  
}

# Store iteration results
results <- bws.vars[nrow(bws.vars), ] # Optimal bandwidths
betas_df <- as.data.frame(betas) %>% mutate(loc=dat.final$loc)

## Final Standard Error estimation
lpl2 <- cox_derivatives_cpp(eta.i.se, time, status, num_threads = num_threads)
dd <- -lpl2$second_derivative
S_arrays_Cpp <- aperm(S.arrays, c(2, 3, 1))
dd[dd == 0] <- 1e-10
beta.se <- calculateBetaSE_Cpp(as.matrix(x1), dd, S_arrays_Cpp) 
res.se <- cbind(beta.se,dat.final$fips)

MGWCR_bws_app <- results
MGWCR_betas_app <- betas_df
MGWCR_se_app <- res.se



#############################################################
### Inference - Standard Error and P-value Calculation
#############################################################
# coefficients
colnames(MGWCR_betas_app) <- c("agemo","white0","marriage","income2","NDVI270")
MGWCR_betas_app <- MGWCR_betas_app %>% mutate(loc=as.data.frame(MGWCR_se_app)$V6)
MGWCR_betas_app$NDVI270 <- MGWCR_betas_app$NDVI270*0.1
as.data.frame(MGWCR_betas_app %>% group_by(loc) %>% slice(1))

# standard error
MGWCR_se_app<- as.data.frame(MGWCR_se_app)
sorted_data <- MGWCR_se_app[order(MGWCR_se_app$V6, rowSums(is.na(MGWCR_se_app)), decreasing = FALSE), ]
# as.data.frame(sorted_data %>% group_by(V6) %>% slice(1))

# p-value
loc_counts <- as.data.frame(table(MGWCR_betas_app$loc)) %>% mutate(loc=as.numeric(as.character(Var1))) %>% rename(count=Freq)
betas_res <- as.data.frame(MGWCR_betas_app %>% group_by(loc) %>% slice(1))
se_res <- as.data.frame(sorted_data %>% group_by(V6) %>% slice(1)) %>% rename(loc = V6)
MGWCR_res <- betas_res %>%
  left_join(se_res,by="loc") %>%
  left_join(loc_counts, by="loc") %>%
  na.omit()

k <- 5
n_tests <- 5
MGWCR_res$p_beta1 <- 2 * (1-pt(abs(MGWCR_res$agemo/MGWCR_res$V1),df=sum(MGWCR_res$count)-k-1))
MGWCR_res$p_beta2 <- 2 * (1-pt(abs(MGWCR_res$white0/MGWCR_res$V2),df=sum(MGWCR_res$count)-k-1))
MGWCR_res$p_beta3 <- 2 * (1-pt(abs(MGWCR_res$marriage/MGWCR_res$V3),df=sum(MGWCR_res$count)-k-1))
MGWCR_res$p_beta4 <- 2 * (1-pt(abs(MGWCR_res$income2/MGWCR_res$V4),df=sum(MGWCR_res$count)-k-1))
MGWCR_res$p_beta5 <- 2 * (1-pt(abs(MGWCR_res$NDVI90/MGWCR_res$V5),df=sum(MGWCR_res$count)-k-1))
# MGWCR_res$p_beta6 <- 2 * (1-pt(abs(MGWCR_res$NDVI270/MGWCR_res$V6),df=sum(MGWCR_res$count)-k-1))
# print(MGWCR_res[, c("p_beta1", "p_beta2", "p_beta3", "p_beta4", "p_beta5","loc")])

# Bonferroni correction
MGWCR_res$p_beta1_bonf <- pmin(MGWCR_res$p_beta1 * n_tests, 1)
MGWCR_res$p_beta2_bonf <- pmin(MGWCR_res$p_beta2 * n_tests, 1)
MGWCR_res$p_beta3_bonf <- pmin(MGWCR_res$p_beta3 * n_tests, 1)
MGWCR_res$p_beta4_bonf <- pmin(MGWCR_res$p_beta4 * n_tests, 1)
MGWCR_res$p_beta5_bonf <- pmin(MGWCR_res$p_beta5 * n_tests, 1)
# MGWCR_res$p_beta6_bonf <- pmin(MGWCR_res$p_beta6 * n_tests, 1)

print(MGWCR_res[, c("p_beta1_bonf", "p_beta2_bonf", "p_beta3_bonf", "p_beta4_bonf", "p_beta5_bonf", "loc")])
print(MGWCR_res[, c("p_beta1_bonf", "p_beta2_bonf", "p_beta3_bonf", "p_beta4_bonf", "p_beta5_bonf", "p_beta6_bonf", "loc")])

