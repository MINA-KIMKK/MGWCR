###################################################################
## Global Cox regression model
## Simulation Scenario 2 - Varied Spatial Heterogeneity (Constant, Medium, High)
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



registerDoParallel(cores = 8)

data.all <- foreach(rr = 1:1000,.combine=rbind,.packages=c("doParallel","dplyr","survival")) %dopar% {
  
  set.seed(rr)
  data1 <- data.generator1(N = 1250,
                           lam = 1,
                           probC = 0,
                           tau = 10)
  
  true.beta <- data1 %>% select(beta1,beta2,beta3)
  
  glb.cox <- coxph(Surv(estop, estatus==1) ~ X1 + X2 + X3, data1)
  glb.coef1 <- glb.cox$coefficients[1]
  glb.coef2 <- glb.cox$coefficients[2]
  glb.coef3 <- glb.cox$coefficients[3]
  
  
  return(data.frame(rr,data1$loc,glb.coef1,glb.coef2,glb.coef3,true.beta,t(sqrt(diag(vcov(glb.cox))))))
}

names(data.all) <- c("r","loc","glb.coef1","glb.coef2","glb.coef3","beta1","beta2","beta3","se.glb.coef1","se.glb.coef2","se.glb.coef3")
data.cp <- data.all %>%
  mutate(bias.glb.beta1 = abs(glb.coef1 - beta1),
         bias.glb.beta2 = abs(glb.coef2 - beta2),
         bias.glb.beta3 = abs(glb.coef3 - beta3))


## estimates
betas_glb2 <- as.data.frame(data.cp %>% group_by(loc) %>% 
                summarise(mean(glb.coef1),mean(glb.coef2),mean(glb.coef3)))
names(betas_glb2) <- c("r","V1","V2","V3")

## true beta
betas_true2 <- as.data.frame(data.cp %>% group_by(loc) %>% 
                summarise(mean(beta1),mean(beta2),mean(beta3)))

## mab,mse
mse_glb2 <- as.data.frame(data.cp %>% group_by(rr) %>% 
                summarise(mean(bias.glb.beta1^2),mean(bias.glb.beta2^2),mean(bias.glb.beta3^2)))
names(mse_glb2) <- c("r","glb_mse1","glb_mse2","glb_mse3")
mse_glb2 %>% select(-r)

## standard error
# empirical sd
sdbetas_glb2 <- as.data.frame(data.cp %>% group_by(loc) %>% 
                summarise(sd(glb.coef1), sd(glb.coef2), sd(glb.coef3)))

# average se
seavg_glb2 <- as.data.frame(data.cp %>% group_by(loc) %>% 
                summarise(mean(se.glb.coef1), mean(se.glb.coef2), mean(se.glb.coef3)))
