# HOIF-Car

The *HOIFCar* package facilitates the estimation of treatment effects through a range of covariate adjustment methods. It supports the calculation of point estimates and variance estimates for treatment effects. Additionally, the package provides oracle bias, oracle variance, and the corresponding approximated variance for the adjusted estimators derived from higher-order influence functions (HOIF).
   

Complete derivation for the exact variance of $\hat{\tau}_{\mathsf{adj}, 2}^{\dagger}$ under the randomization-based framework can be found [here](https://github.com/Cinbo-Wang/HOIF-Car/blob/main/var-db.pdf).

## Installation
```R
devtools::install_github("Cinbo-Wang/HOIF-Car")
```

## Example

### adj2 and adj2c estimators
We first consider the Completely Randomized Experiment (CRE) with moderately high-dimensional covariates under the randomization-based (design-based) framework. 
```R
# Oracle Estimatation 
## Linear setting
set.seed(100)
n <- 500
p <- 50
beta <- rt(p,3)

X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 3)
Y1 <- as.numeric(X %*% beta)
pi1 <- 0.5
n1 <- ceiling(n*pi1)

require(HOIFCar)
result_adj_db <- get_oracle_bias_var_adj_db(X = X,Y1 = Y1,n1 = n1) # from LYW paper
result_adj2c <- get_oracle_bias_var_adj2c(X = X,Y1 = Y1,n1 = n1)
result_adj2_3 <- get_oracle_bias_var_adj_2_3(X = X, Y1 = Y1,n1 = n1)
unlist(result_adj_db)
unlist(result_adj2c)
unlist(result_adj2_3)



## Nonlinear setting
n <- 500;
alpha <- 0.2;
set.seed(1000)
p <- ceiling(n * alpha)
Sigma_true <- matrix(0, nrow = p, ncol = p)
for(i in 1:p){
  for(j in 1:p){
    Sigma_true[i, j] <- 0.1**(abs(i - j))
  }
}

X <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)
beta <- rt(p,3)
or_baseline <- sign(X %*% beta) * abs(X %*% beta)^(1/2) + sin(X %*% beta)
epsilon1 <- epsilon0 <- rt(n, 3)
Y1 <- 1 + as.numeric(or_baseline) + epsilon1


pi1 <- 2/3
n1 <- ceiling(n * pi1)

result_adj_db <- get_oracle_bias_var_adj_db(X = X, Y1 = Y1, n1 = n1) # from LYW paper
result_adj2c <- get_oracle_bias_var_adj2c(X = X, Y1 = Y1, n1 = n1)
result_adj2_3 <- get_oracle_bias_var_adj_2_3(X = X, Y1 = Y1, n1 = n1)
unlist(result_adj_db)
unlist(result_adj2c)
unlist(result_adj2_3)



# Realistic estimation
set.seed(100)
n <- 500
p <- n * 0.3
beta <- runif(p, -1 / sqrt(p), 1 / sqrt(p))

X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 3)
Y1 <- as.numeric(X %*% beta)
Y0 <- rep(0, n)

pi1 <- 2/3
n1 <- ceiling(n * pi1)
ind <- sample(n, size = n1)
A <- rep(0, n)
A[ind] <- 1
Y <- Y1 * A + Y0 * (1 - A)

Xc <- scale(X, scale = FALSE)
H <- Xc %*% MASS::ginv(t(Xc) %*% Xc) %*% t(Xc)

## Estimate mean treat on the treatment arm
result_ls <- esti_mean_treat(X, Y, A, H)
point_est <- result_ls$point_est
var_est <- result_ls$var_est
print(paste0('True mean treat:', round(mean(Y1), digits = 3), '.'))
print('Absolute bias:')
print(abs(point_est - mean(Y1)))
print('Estimated variance:')
print(var_est)

## Estimate ATE using adj2 and adj2c estimators
Xc <- cbind(1, scale(X, scale = FALSE))
result.adj2.adj2c.random.ls <-  fit.adj2.adj2c.Random(Y, Xc, A)
result.adj2.adj2c.random.ls

```
We next consider the simple randomization experiments with moderately high-dimensional covariates under the Super-Population based framework. 
```R
set.seed(120)
alpha0 <- 0.1; 
n <- 400; 

p0 <- ceiling(n * alpha0)
beta0_full <- 1 / (1:p0) ^ (1 / 2) * (-1) ^ c(1:p0)
beta <- beta0_full / norm(beta0_full,type='2')

Sigma_true <- matrix(0, nrow = p0, ncol = p0)
for (i in 1:p0) {
  for (j in 1:p0) {
    Sigma_true[i, j] <- 0.1 ** (abs(i - j))
  }
}

X <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)

lp0 <- X %*% beta
delta_X <- 1  -  1/4 * X[, 2] -  1/8 * X[, 3]
lp1 <- lp0 + delta_X

Y0 <- lp0 + rnorm(n)
Y1 <- lp1 + rnorm(n)


pi1 <- 1 / 2
A <- rbinom(n, size = 1, prob = pi1)
Y <- A * Y1 + (1 - A) * Y0


# Estimate ATE, EY1 and EY0 using adj2 and adj2c estimators
Xc <- cbind(1, scale(X, scale = FALSE))
result.adj2.adj2c.sp.ate.ls <- fit.adj2.adj2c.Super(Y, Xc, A, intercept = TRUE,
                                                    target = 'ATE', lc = TRUE)
result.adj2.adj2c.sp.ate.ls
result.adj2.adj2c.sp.treat.ls <- fit.adj2.adj2c.Super(Y, Xc, A, intercept = TRUE,
                                                      target = 'EY1', lc = TRUE)
result.adj2.adj2c.sp.treat.ls
result.adj2.adj2c.sp.control.ls <- fit.adj2.adj2c.Super(Y, Xc, A, intercept = TRUE,
                                                        target = 'EY0', lc = TRUE)
result.adj2.adj2c.sp.control.ls
```
We next consider the Covariate-adaptive randomization with moderately high-dimensional covariates under the Super-Population based framework. 
```R
  set.seed(120)
  alpha0 <- 0.1;
  n <- 400;
  S <- as.factor(sample(c("0-30","31-50",">50"),n,replace = T,prob=c(0.2,0.4,0.4)))
  ns_min = min(table(S))
  
  p0 <- ceiling(ns_min * alpha0)
  beta0_full <- 1 / (1:p0) ^ (1 / 2) * (-1) ^ c(1:p0)
  beta <- beta0_full / norm(beta0_full,type='2')
  
  Sigma_true <- matrix(0, nrow = p0, ncol = p0)
  for (i in 1:p0) {
    for (j in 1:p0) {
      Sigma_true[i, j] <- 0.1 ** (abs(i - j))
    }
  }
  
  X <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)
  
  lp0 <- X %*% beta
  delta_X <- 1  -  1/4 * X[, 2] -  1/8 * X[, 3]
  lp1 <- lp0 + delta_X
  
  Y0 <- lp0 + rnorm(n)
  Y1 <- lp1 + rnorm(n)
  
  
  pi1 <- 1 / 2
  
  # We use stratified block randomization as an example, simple randomization
  # is also valid by setting S = rep(1,n) and A = rbinom(n,1,pi1)
  
  sbr <- function(S,nA,p,block_size=10){
    N <- length(S)
    B <- block_size
    A <- rep(0,N)
    nS <- length(unique(S))
    for(s in 1:nS){
      ind_s <- which(S==s)
      n_s <- length(ind_s)
      A_s <- rep(0,n_s)
      numB <- floor(n_s/B)
      rem <- n_s - numB*B
      size_A <- B*p[s]
      if(numB==0){
        size_rem = floor(rem*p[s])
        size_rem[1] = rem - sum(size_rem[-1])
        A_s[(B*numB+1):n_s] <- sample(rep(0:(nA-1),size_rem),size=rem,replace = F)
      }else{
        for(i in 1:numB){
          A_s[(B*(i-1)+1):(B*i)] <- sample(rep(0:(nA-1),size_A),size=B,replace = F)
        }
        if(rem>0){
          size_rem = floor(rem*p[s])
          size_rem[1] = rem - sum(size_rem[-1])
          A_s[(B*numB+1):n_s] <- sample(rep(0:(nA-1),size_rem),size=rem,replace = F)
        }
      }
      A[ind_s] <- A_s
    }
    return(A)
  }
  
  
  A <- sbr(as.numeric(S),2,rep(pi1,3),block_size = 4)
  
  Y <- A * Y1 + (1 - A) * Y0
  
  Xc <- cbind(1, scale(X, scale = FALSE))
  result.adj2.adj2c.car.ate.ls <- fit.adj2.adj2c.CAR(Y, Xc,S, A, intercept = TRUE,
                                                     target = 'ATE')
  result.adj2.adj2c.car.ate.ls
  result.adj2.adj2c.car.treat.ls <- fit.adj2.adj2c.CAR(Y, Xc,S, A, intercept = TRUE,
                                                       target = 'EY1')
  result.adj2.adj2c.car.treat.ls
  result.adj2.adj2c.car.control.ls <- fit.adj2.adj2c.CAR(Y, Xc,S, A, intercept = TRUE,
                                                         target = 'EY0')
  result.adj2.adj2c.car.control.ls

```

### JASA and JASA-cal estimators
```R
generate_data_SR <- function(n, family, pi1, p_n_ratio = 0.05, seed = 123){
  set.seed(seed)
  alpha0 <- 0.15
  p0 <- ceiling(round(n * alpha0))
  beta0_full <- 1/(1:p0)^(1/4)*(-1)^c(1:p0)
  Sigma_true <- matrix(0, nrow = p0, ncol = p0)
  for (i in 1:p0) {
    for (j in 1:p0) {
      Sigma_true[i, j] <- 0.1 ** (abs(i - j))
    }
  }

  if(family != 'poisson'){
    X <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)
  }else{
    X0 <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)
    X <- pmin(pmax(X0, -3), 3)
    rm(X0)
  }

  beta <- beta0_full / norm(beta0_full,type='2')
  lp0 <- X %*% beta

  delta_X <- 1 - 1/2 * pmin(X[, 1]^2, 5) + 1/4 * X[, 1:10] %*% beta[1:10]
  lp1 <- lp0 + delta_X


  if (family == 'binomial') {
    r0 <- plogis(2 * lp0)
    r1 <- plogis(2 * lp1)
    Y1 <- rbinom(n, size=1, prob=r1)
    Y0 <- rbinom(n, size=1, prob=r0)
  }else if(family == 'poisson'){
    # quantile(lp1);quantile(lp0)
    lp1_tran <- pmin(lp1, 4)
    lp0_tran <- pmin(lp0, 4)
    r1 <- exp(lp1_tran)
    r0 <- exp(lp0_tran)

    Y1 <- rpois(n,r1)
    Y0 <- rpois(n,r0)
  }else if(family == 'gaussian'){
    r1 <- lp1;
    r0 <- lp0
    Y1 <- r1 + rnorm(n)
    Y0 <- r0 + rnorm(n)
  }

  A <- rbinom(n, size=1, prob=pi1)
  Y <- A * Y1 + (1 - A) * Y0

  p <- ceiling(round(n * p_n_ratio))
  if(p > ncol(X)){
    if(family != 'poisson'){
      X_noise <- rmvt(n, sigma = diag(p - ncol(X)), df = 3)
    }else{
      X0_noise <- rmvt(n, sigma = diag(p - ncol(X)), df = 3)
      X_noise <- pmin(pmax(X0_noise, -3), 3)
    }
    X_obs <- cbind(X, X_noise)
  }else{
    X_obs <- X[, 1:p, drop = FALSE]
  }

  data_ls <- list(
    X = X_obs, Y = Y, A = A,
    Y1 = Y1, Y0 = Y0,
    r1 = r1, r0 = r0
  )
  return(data_ls)
}


n <- 400; pi1 <- 1/3

## Continuous outcomes
family <- 'gaussian'; p_n_ratio <- 0.05
data_ls <- generate_data_SR(n, family, pi1, p_n_ratio)
X <- data_ls$X;
A <- data_ls$A
Y <- data_ls$Y

Xc <- scale(X, scale = FALSE)
Xc_aug <- cbind(1, Xc)
result.jasa.ls <- fit.JASA(Y, Xc_aug, A, family, opt_obj = 'beta')
result.jasa.ls$tau_vec
result.jasa.ls$var_tau_vec

## Count outcomes
family <- 'poisson'; p_n_ratio <- 0.05
data_ls <- generate_data_SR(n, family, pi1, p_n_ratio)
X <- data_ls$X;
A <- data_ls$A
Y <- data_ls$Y

Xc <- scale(X, scale = FALSE)
Xc_aug <- cbind(1, Xc)
result.jasa.ls <- fit.JASA(Y, Xc_aug, A, family, opt_obj = 'mu')
result.jasa.ls$tau_vec
result.jasa.ls$var_tau_vec

## Binary outcomes
family <- 'binomial'; p_n_ratio <- 0.05
data_ls <- generate_data_SR(n, family, pi1, p_n_ratio)
X <- data_ls$X;
A <- data_ls$A
Y <- data_ls$Y

Xc <- scale(X, scale = FALSE)
Xc_aug <- cbind(1, Xc)
result.jasa.ls <- fit.JASA(Y, Xc_aug, A, family, opt_obj = 'beta')
result.jasa.ls$tau_vec
result.jasa.ls$var_tau_vec
```


## References
Zhao, S., Wang, X., Liu, L., & Zhang, X. (2024). Covariate adjustment in randomized experiments motivated by higher-order influence functions. arXiv preprint. https://arxiv.org/abs/2411.08491

Gu, Y., Liu, L., & Ma, W. (2025). Assumption-lean covariate adjustment under covariate adaptive randomization when $ p= o (n) $. arXiv preprint arXiv:2512.20046.

Xin Lu, Fan Yang, and Yuhao Wang. Debiased regression adjustment in completely randomized experiments with moderately high-dimensional covariates. arXiv preprint arXiv:2309.02073, 2023.

Marlena S. Bannick, Jun Shao, Jingyi Liu, Yu Du, Yanyao Yi, Ting Ye (2025). A General Form of Covariate Adjustment in Randomized Clinical Trials. Biometrika.

