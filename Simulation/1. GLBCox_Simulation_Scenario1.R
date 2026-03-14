###################################################################
## Global Cox regression model
## Simulation Scenario 1 - Same Variation Scale 
###################################################################
library(MASS)
library(foreach)
library(doParallel)
library(dplyr)
library(survival)
select <- dplyr::select
mutate <- dplyr::mutate
rename <- dplyr::rename 
summarise <- dplyr::summarise 


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


# Run simulation loop (1000 iterations)
registerDoParallel(cores = 8)

data.all <- foreach(rr = 1:1000,.combine=rbind,.packages=c("doParallel","dplyr","survival")) %dopar% {
  
  set.seed(rr)
  data1 <- data.generator2(N = 1250,
                           lam = 1,
                           probC = 0,
                           tau = 10)
  
  true.beta <- data1 %>% select(beta1,beta2)
  
  glb.cox <- coxph(Surv(estop, estatus==1) ~ X1 + X2, data1)
  glb.coef1 <- glb.cox$coefficients[1]
  glb.coef2 <- glb.cox$coefficients[2]
  
  return(data.frame(rr,data1$loc,glb.coef1,glb.coef2,true.beta,t(sqrt(diag(vcov(glb.cox))))))
}

names(data.all) <- c("r","loc","glb.coef1","glb.coef2","beta1","beta2","se.glb.coef1","se.glb.coef2")
data.cp <- data.all %>%
  mutate(bias.glb.beta1 = abs(glb.coef1 - beta1),
         bias.glb.beta2 = abs(glb.coef2 - beta2))

## estimates
betas_glb1 <- as.data.frame(data.cp %>% group_by(loc) %>% 
                summarise(mean(glb.coef1),mean(glb.coef2)))
names(betas_glb1) <- c("r","V1","V2")

## true beta
betas_true1 <- as.data.frame(data.cp %>% group_by(loc) %>% 
                summarise(mean(beta1),mean(beta2)))


## MSE
mse_glb1 <- as.data.frame(data.cp %>% group_by(r) %>% 
                summarise(mean(bias.glb.beta1^2), mean(bias.glb.beta2^2)))
names(mse_glb1) <- c("r","glb_mse1","glb_mse2")
mse_glb1 %>% select(-r)

## standard error
# empirical sd
sdbetas_glb1 <- as.data.frame(data.cp %>% group_by(loc) %>% 
                summarise(sd(glb.coef1), sd(glb.coef2)))

# average se
seavg_glb1 <- as.data.frame(data.cp %>% group_by(loc) %>% 
                summarise(mean(se.glb.coef1), mean(se.glb.coef2)))
