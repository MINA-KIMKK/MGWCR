###################################################################
## GW Cox: Geographically Weighted Cox Regression
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
max.iterations=100
threshold=10^-5
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
## Initialization with Cox estimates(global model)
######################################################################

cox.glb <- coxph(formula,data=data.ori)
beta0 <- cox.glb$coefficients
betas <- matrix(rep(beta0, nrow(dp.locat)), nrow = nrow(dp.locat), byrow = TRUE)
bws0 <- rep(max(dMat), var.n)

# Initialize the additive term fj
f.i <- betas*x1


#######################################################
## Backfitting algorithm 
#######################################################

iteration = 0 
criterion.val <- 10000000  

while ((iteration < max.iterations) && criterion.val > threshold){
  print(paste0("iteration: ", iteration+1," / ","criterion.val: ",criterion.val))
  # Update the linear predictor (eta) and linearized dependent variable (Z)
  eta.i <- rowSums(betas*x1)
  lpl <- cox_derivatives_cpp(eta.i, time, status, num_threads = num_threads) 
  dd <- -lpl$second_derivative
  dd[dd == 0] <- 1e-10
  d.i <- diag(dd)
  u.i <- lpl$first_derivative
  z.i <- eta.i + diag(1/dd)%*%u.i
  gc()
  
  f.i_old <- f.i
  diffi <- 0
  f.i2 <- 0
  #Global spatial bandwidth selection
  bw1<-bw.gwr2(as.matrix(x1), z.i,  dp.locat, dd, approach = approach, kernel=kernel.id, adaptive = adaptive, dMat, verbose = verbose, nlower = nlower,parallel.method=parallel.method,parallel.arg=parallel.arg)
  # Weighted regression to update local coefficients
  res <- gwr.q2(as.matrix(x1), z.i, d.i, dp.locat, adaptive=adaptive, hatmatrix = hatmatrix, bw=bw1, kernel=kernel.id,dMat=dMat)
  betas <- res[[1]]
  bws.vars <- c(bws.vars, bw1)
  print(bw1)
  
  f.i <- betas * x1
  diffi <- sum(colMeans((f.i - f.i_old)^2))
  f.i2 <- rowSums(f.i)
  
  socf <- sqrt(diffi/sum(f.i2^2)) 
  criterion.val <- socf

  iteration <- iteration+1 
}
# Optimal single bandwidth for the GW Cox model
results <-  bws.vars[length(bws.vars)] 
betas_df <- as.data.frame(betas) %>% mutate(loc=dat.final$fips)

## Standard Error Estimation
S.arrays <- array(dim=c(var.n, dp.n, dp.n))
for (i in 1:var.n){
  res <- gwr.q2(matrix(x1[,i], ncol=1), z.i, d.i, dp.locat, adaptive=adaptive, hatmatrix = hatmatrix, bw=results, kernel=kernel.id,dMat=dMat)
  S.arrays[i,,] <- res[[2]] # Rj
  
}
S_arrays_Cpp <- aperm(S.arrays, c(2, 3, 1)) 
beta.se <- calculateBetaSE_Cpp(as.matrix(x1), dd, S_arrays_Cpp)
res.se <- cbind(beta.se,dat.final$fips)

GWCR_bws_app <- results
GWCR_betas_app <- betas_df
GWCR_betas_se_app <- res.se


#############################################################
### Inference - Standard Error and P-value Calculation
#############################################################

# coefficients
GWCR_betas_app$NDVI270 <- GWCR_betas_app$NDVI270*0.1
as.data.frame(GWCR_betas_app %>% group_by(loc) %>% slice(1))

# standard error
GWCR_se_app<- as.data.frame(GWCR_betas_se_app)
sorted_data <- GWCR_se_app[order(GWCR_se_app$V6, rowSums(is.na(GWCR_se_app)), decreasing = FALSE), ]
as.data.frame(sorted_data %>% group_by(V7) %>% slice(1))

# p-value
loc_counts <- as.data.frame(table(GWCR_betas_app$loc)) %>% mutate(loc=as.numeric(as.character(Var1))) %>% rename(count=Freq)
betas_res <- as.data.frame(GWCR_betas_app %>% group_by(loc) %>% slice(1))
se_res <- as.data.frame(sorted_data %>% group_by(V6) %>% slice(1)) %>% rename(loc = V6)
colnames(betas_res) <- c("agemo","marriage","income2","pm25","NDVI270","loc")
GWCR_res <- betas_res %>%
  left_join(se_res,by="loc") %>%
  rename(beta1=agemo,beta2=marriage,beta3=income2,beta4=pm25,beta5=NDVI270,
         se1=V1,se2=V2,se3=V3,se4=V4,se5=V5) %>%
  left_join(loc_counts, by="loc") %>%
  na.omit()

k <- 5
n_tests <- 5

GWCR_res$p_beta1 <- 2 * (1-pt(abs(GWCR_res$beta1/GWCR_res$se1),df=sum(GWCR_res$count)-k-1))
GWCR_res$p_beta2 <- 2 * (1-pt(abs(GWCR_res$beta2/GWCR_res$se2),df=sum(GWCR_res$count)-k-1))
GWCR_res$p_beta3 <- 2 * (1-pt(abs(GWCR_res$beta3/GWCR_res$se3),df=sum(GWCR_res$count)-k-1))
GWCR_res$p_beta4 <- 2 * (1-pt(abs(GWCR_res$beta4/GWCR_res$se4),df=sum(GWCR_res$count)-k-1))
GWCR_res$p_beta5 <- 2 * (1-pt(abs(GWCR_res$beta5/GWCR_res$se5),df=sum(GWCR_res$count)-k-1))
GWCR_res$p_beta6 <- 2 * (1-pt(abs(GWCR_res$beta5/GWCR_res$se6),df=sum(GWCR_res$count)-k-1))

# Bonferroni correction
GWCR_res$p_beta1_bonf <- pmin(GWCR_res$p_beta1 * n_tests, 1)
GWCR_res$p_beta2_bonf <- pmin(GWCR_res$p_beta2 * n_tests, 1)
GWCR_res$p_beta3_bonf <- pmin(GWCR_res$p_beta3 * n_tests, 1)
GWCR_res$p_beta4_bonf <- pmin(GWCR_res$p_beta4 * n_tests, 1)
GWCR_res$p_beta5_bonf <- pmin(GWCR_res$p_beta5 * n_tests, 1)
GWCR_res$p_beta6_bonf <- pmin(GWCR_res$p_beta6 * n_tests, 1)

print(GWCR_res[, c("p_beta1_bonf", "p_beta2_bonf", "p_beta3_bonf", "p_beta4_bonf", "p_beta5_bonf", "loc")])
print(GWCR_res[, c("p_beta1_bonf", "p_beta2_bonf", "p_beta3_bonf", "p_beta4_bonf", "p_beta5_bonf", "p_beta6_bonf", "loc")])

