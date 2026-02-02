
#' Covariate Adjustment via JAckknife Score-based Adjustment (JASA) for Generalized Linear Working Models
#'
#' @description
#' Implements the Jackknife Score-Based Adjustment (JASA) method and its calibration for covariate adjustment in simple randomized experiments where covariate dimension p may be large relatively to sample size n. Handles Continuous, Binary, and Poisson outcomes.
#'
#' @param Y Numeric vector of outcome values.
#' @param X Matrix of centered baseline covariates (may include intercept). Dimensions: n x p.
#' @param A Binary treatment vector (0 = control, 1 = treatment). Assignment is assumed to follow simple randomization.
#' @param family  GLM family specification: \code{"gaussian"}, \code{"binomial"}, or \code{"poisson"}. Default: \code{"gaussian"}.
#' @param pi1 The assignment probability for the simple randomization. If NULL (default), the empirical assignment probability is used.
#' @param is.parallel Boolean for parallelization. Default: FALSE.
#' @param core_num Number of cores for parallel computing (default is 4 if is.parallel = TRUE).
#' @param opt_obj  Ways to optimization: 'beta' (GLM coefficients) or 'mu' (expected outcomes). Default: 'beta'.
#'
#' @return List with components:
#' \item{tau_vec}{Named vector of average treatment effect estimates (\code{JASA}, \code{JASA-cal})}
#' \item{tau0_vec}{Treatment effect estimates on the control group}
#' \item{tau1_vec}{Treatment effect estimates on the treatment group}
#' \item{var_tau_vec}{Variance estimates for average treatment effect estimates }
#' \item{var_tau0_vec}{Variance estimates for treatment effect estimates on the control group}
#' \item{var_tau1_vec}{Variance estimates for treatment effect estimates on the treatment group}
#' \item{y_hat_mat}{Matrix of predicted outcomes (columns: control, treatment) using varied Leave-One-Out strategy.}
#' \item{obj_value_mat}{Matrix of objective values from optimization (columns: control, treatment) }
#'
#' @export
#'
#' @examples
#' generate_data_SR <- function(n, family, pi1, p_n_ratio = 0.05, seed = 123){
#'   set.seed(seed)
#'   alpha0 <- 0.15
#'   p0 <- ceiling(round(n * alpha0))
#'   beta0_full <- 1/(1:p0)^(1/4)*(-1)^c(1:p0)
#'   Sigma_true <- matrix(0, nrow = p0, ncol = p0)
#'   for (i in 1:p0) {
#'     for (j in 1:p0) {
#'       Sigma_true[i, j] <- 0.1 ** (abs(i - j))
#'     }
#'   }
#'
#'   if(family != 'poisson'){
#'     X <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)
#'   }else{
#'     X0 <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)
#'     X <- pmin(pmax(X0, -3), 3)
#'     rm(X0)
#'   }
#'
#'   beta <- beta0_full / norm(beta0_full,type='2')
#'   lp0 <- X %*% beta
#'
#'   delta_X <- 1 - 1/2 * pmin(X[, 1]^2, 5) + 1/4 * X[, 1:10] %*% beta[1:10]
#'   lp1 <- lp0 + delta_X
#'
#'
#'   if (family == 'binomial') {
#'     r0 <- plogis(2 * lp0)
#'     r1 <- plogis(2 * lp1)
#'     Y1 <- rbinom(n, size=1, prob=r1)
#'     Y0 <- rbinom(n, size=1, prob=r0)
#'   }else if(family == 'poisson'){
#'     # quantile(lp1);quantile(lp0)
#'     lp1_tran <- pmin(lp1, 4)
#'     lp0_tran <- pmin(lp0, 4)
#'     r1 <- exp(lp1_tran)
#'     r0 <- exp(lp0_tran)
#'
#'     Y1 <- rpois(n,r1)
#'     Y0 <- rpois(n,r0)
#'   }else if(family == 'gaussian'){
#'     r1 <- lp1;
#'     r0 <- lp0
#'     Y1 <- r1 + rnorm(n)
#'     Y0 <- r0 + rnorm(n)
#'   }
#'
#'   A <- rbinom(n, size=1, prob=pi1)
#'   Y <- A * Y1 + (1 - A) * Y0
#'
#'   p <- ceiling(round(n * p_n_ratio))
#'   if(p > ncol(X)){
#'     if(family != 'poisson'){
#'       X_noise <- rmvt(n, sigma = diag(p - ncol(X)), df = 3)
#'     }else{
#'       X0_noise <- rmvt(n, sigma = diag(p - ncol(X)), df = 3)
#'       X_noise <- pmin(pmax(X0_noise, -3), 3)
#'     }
#'     X_obs <- cbind(X, X_noise)
#'   }else{
#'     X_obs <- X[, 1:p, drop = FALSE]
#'   }
#'
#'   data_ls <- list(
#'     X = X_obs, Y = Y, A = A,
#'     Y1 = Y1, Y0 = Y0,
#'     r1 = r1, r0 = r0
#'   )
#'   return(data_ls)
#' }
#'
#'
#' n <- 400; pi1 <- 1/3
#'
#' family <- 'gaussian'; p_n_ratio <- 0.05
#' data_ls <- generate_data_SR(n, family, pi1, p_n_ratio)
#' X <- data_ls$X;
#' A <- data_ls$A
#' Y <- data_ls$Y
#'
#' \dontrun{
#' Xc <- scale(X, scale = FALSE)
#' Xc_aug <- cbind(1, Xc)
#' result.jasa.ls <- fit.JASA(Y, Xc_aug, A, family, opt_obj = 'beta')
#' result.jasa.ls$tau_vec
#' result.jasa.ls$var_tau_vec
#'
#'
#' family <- 'poisson'; p_n_ratio <- 0.05
#' data_ls <- generate_data_SR(n, family, pi1, p_n_ratio)
#' X <- data_ls$X;
#' A <- data_ls$A
#' Y <- data_ls$Y
#'
#' Xc <- scale(X, scale = FALSE)
#' Xc_aug <- cbind(1, Xc)
#' result.jasa.ls <- fit.JASA(Y, Xc_aug, A, family, opt_obj = 'mu')
#' result.jasa.ls$tau_vec
#' result.jasa.ls$var_tau_vec
#'
#'
#' family <- 'binomial'; p_n_ratio <- 0.05
#' data_ls <- generate_data_SR(n, family, pi1, p_n_ratio)
#' X <- data_ls$X;
#' A <- data_ls$A
#' Y <- data_ls$Y
#'
#' Xc <- scale(X, scale = FALSE)
#' Xc_aug <- cbind(1, Xc)
#' result.jasa.ls <- fit.JASA(Y, Xc_aug, A, family, opt_obj = 'beta')
#' result.jasa.ls$tau_vec
#' result.jasa.ls$var_tau_vec
#' }
#' @import parallel
#' @import glmnet
#' @import foreach
#' @import doParallel
#' @import stats
#' @import brglm2
fit.JASA <- function(Y, X, A, family = 'gaussian', pi1 = NULL,
                     is.parallel = FALSE, core_num = 4, opt_obj = c('beta', 'mu')[1]){
  # require(doParallel);require(foreach);require(stats); require(brglm2)

  if(is.null(pi1)){
    pi1_hat <- mean(A)
  }else{
    pi1_hat <- pi1
  }

  N <- length(Y)
  data_df <- data.frame(cbind(Y, X))
  colnames(data_df) <- c('Y', paste0('X', 1:ncol(X)))

  fit.glm <- stats::glm(
    Y ~ -1 + .,
    data = data_df,
    family = family,
    method = brglm2::brglmFit,
    control = brglm2::brglmControl(
      type = 'AS_mixed',
      maxit = 1000,
      slowit = 0.1,
      epsilon = 1e-4
    )
  )
  beta1_hat <- beta0_hat <- as.numeric(stats::coef(fit.glm))

  ## Update  initial estimates
  fit1.firth <- stats::glm(
    Y ~ -1 + .,
    data = data_df,
    family = family,
    method = brglm2::brglmFit,
    subset = (A == 1),
    control = brglm2::brglmControl(type = 'AS_mixed', maxit = 500, slowit = 0.1)
  )
  beta1_firth <- as.numeric(stats::coef(fit1.firth))
  if(sum(is.na(beta1_firth)) == 0){
    if(fit1.firth$converged & (max(abs(beta1_firth)) < 1e5)){
      beta1_hat <- beta1_firth
    }
  }


  fit0.firth <- stats::glm(
    Y ~ -1 + .,
    data = data_df,
    family = family,
    method = brglm2::brglmFit,
    subset = (A == 0),
    control = brglm2::brglmControl(type = 'AS_mixed', maxit = 500, slowit = 0.1)
  )
  beta0_firth <- as.numeric(stats::coef(fit0.firth))
  if(sum(is.na(beta0_firth)) == 0){
    if(fit0.firth$converged & (max(abs(beta0_firth)) < 1e5)){
      beta0_hat <- beta0_firth
    }
  }



  if(is.parallel){
    type <- ifelse(.Platform$OS.type=='windows','PSOCK','FORK')
    cl <- parallel::makeCluster(core_num,type);
    doParallel::registerDoParallel(cl)

    result.jasa.ls <- foreach::foreach(i = 1:N, .export = c("fit.JASA.i", "fit.JASA.i.poisson",
                                                   "fit.JASA.i.binomial", "fit.JASA.i.gaussian"))%dopar%{
      # source('./fun_JAckknifingScoreAdjustment.R')
      result.jasa.tmp <- fit.JASA.i(Y, X, A, family, pi1_hat, beta1_hat, beta0_hat, i, opt_obj)
      result.jasa.tmp
    }
    doParallel::stopImplicitCluster()
    parallel::stopCluster(cl)
  }else{
    result.jasa.ls <- list()
    for(i in 1:N){
      result.jasa.ls[[i]] <- fit.JASA.i(Y, X, A, family, pi1_hat, beta1_hat, beta0_hat, i, opt_obj)
    }
  }

  y_hat_mat <- matrix(0, nrow = N, ncol = 2)
  obj_value_mat <- matrix(0, nrow = N, ncol = 2)
  for(i in 1:N){
    y_hat_mat[i, ] <- result.jasa.ls[[i]]$yhat_i_negi
    obj_value_mat[i, ] <- result.jasa.ls[[i]]$l0_l1_value
  }

  yc_hat <- y_hat_mat[, 1]
  yt_hat <- y_hat_mat[, 2]

  # calibration
  Xtilde <- data.frame(y = Y, mu_c = yc_hat, mu_t = yt_hat)
  fit0 <- lm(y ~ ., data = Xtilde, subset = (A == 0))
  fit1 <- lm(y ~ ., data = Xtilde, subset = (A == 1))

  yc_hat_cal <- predict(fit0, newdata = Xtilde)
  yt_hat_cal <- predict(fit1, newdata = Xtilde)

  # quantile(yc_hat_cal)
  # quantile(yt_hat_cal)
  p <- ncol(X)
  if(p < 10){
    width <- max(Y) - min(Y)
    if(family == 'poisson'){
      lb <- - 0.1 * width; ub <- max(Y) + 0.1 * width
    }else if(family == 'binomial'){
      lb <- -0.1; ub <- 1.1
    }else if(family == 'gaussian'){
      lb <- min(Y) - 0.1 * width; ub <- max(Y) + 0.1 * width
    }
  }else{
    lb <- -Inf; ub <- Inf
  }


  cond1 <- (min(yc_hat_cal) > lb & max(yc_hat_cal) < ub)
  cond2 <- (min(yt_hat_cal) > lb & max(yt_hat_cal) < ub)
  if(!(cond1 & cond2)){
    yc_hat_cal <- yc_hat
    yt_hat_cal <- yt_hat
  }


  y_hat_cal_mat <- cbind(yc_hat_cal, yt_hat_cal)

  psi1_jasa_vec <- A * Y / pi1_hat + (1 - A / pi1_hat) * yt_hat
  psi0_jasa_vec <- (1 - A) * Y / (1 - pi1_hat) + (1 - (1 - A) / (1 - pi1_hat)) * yc_hat

  psi1_jasa_cal_vec <- A * Y / pi1_hat + (1 - A / pi1_hat) * yt_hat_cal
  psi0_jasa_cal_vec <- (1 - A) * Y / (1 - pi1_hat) + (1 - (1 - A) / (1 - pi1_hat)) * yc_hat_cal


  tau1_jasa <- mean(psi1_jasa_vec, na.rm = TRUE)
  tau1_jasa_cal <- mean(psi1_jasa_cal_vec, na.rm = TRUE)

  tau0_jasa <- mean(psi0_jasa_vec, na.rm = TRUE)
  tau0_jasa_cal <- mean(psi0_jasa_cal_vec, na.rm = TRUE)


  tau0_jasa_vec <- c(tau0_jasa, tau0_jasa_cal)
  tau1_jasa_vec <- c(tau1_jasa, tau1_jasa_cal)
  tau_jasa_vec <- tau1_jasa_vec - tau0_jasa_vec

  names(tau0_jasa_vec) <- names(tau1_jasa_vec) <- names(tau_jasa_vec) <- c('JASA','JASA-cal')


  ### Variance Estimate by Influence function
  var_infl_jasa <- 1 / N * var(psi1_jasa_vec - psi0_jasa_vec, na.rm = TRUE)
  var_infl_jasa_cal <- 1 / N * var(psi1_jasa_cal_vec - psi0_jasa_cal_vec, na.rm = TRUE)
  var_infl_tau_vec <- c(var_infl_jasa, var_infl_jasa_cal)
  names(var_infl_tau_vec) <- c('JASA','JASA-cal')


  var1_infl_jasa <- 1 / N * var(psi1_jasa_vec, na.rm = TRUE)
  var1_infl_jasa_cal <- 1 / N * var(psi1_jasa_cal_vec, na.rm = TRUE)
  var_infl_tau1_vec <- c(var1_infl_jasa, var1_infl_jasa_cal)
  names(var_infl_tau1_vec) <- c('JASA','JASA-cal')

  var0_infl_jasa <- 1 / N * var(psi0_jasa_vec, na.rm = TRUE)
  var0_infl_jasa_cal <- 1 / N * var(psi0_jasa_cal_vec, na.rm = TRUE)
  var_infl_tau0_vec <- c(var0_infl_jasa, var0_infl_jasa_cal)
  names(var_infl_tau0_vec) <- c('JASA','JASA-cal')


  return(list(
    tau_vec = tau_jasa_vec,
    tau0_vec = tau0_jasa_vec,
    tau1_vec = tau1_jasa_vec,
    var_tau_vec = var_infl_tau_vec,
    var_tau0_vec = var_infl_tau0_vec,
    var_tau1_vec = var_infl_tau1_vec,
    y_hat_mat = y_hat_mat,
    obj_value_mat = obj_value_mat
  ))

}

#' @keywords internal
fit.JASA.i <- function(Y, X, A, family, pi1_hat, beta1_hat, beta0_hat, i, opt_obj){

  if(family=='poisson'){
    result.jasa.i <- fit.JASA.i.poisson(Y, X, A, pi1_hat, beta1_hat, beta0_hat, i, opt_obj)
  }else if(family=='binomial'){
    result.jasa.i <- fit.JASA.i.binomial(Y, X, A, pi1_hat, beta1_hat, beta0_hat, i, opt_obj)
  }else if(family=='gaussian'){
    result.jasa.i <- fit.JASA.i.gaussian(Y, X, A, pi1_hat, beta1_hat, beta0_hat, i, opt_obj)
  }

  return(list(
    yhat_i_negi = result.jasa.i$yhat_i_negi,
    l0_l1_value = result.jasa.i$l0_l1_value
  ))
}



#' @keywords internal
#' @import stats
#' @import BB
fit.JASA.i.binomial <- function(Y, X, A, pi1_hat, beta1_hat, beta0_hat, i, opt_obj){
  # require(stats); require(BB)
  n <- length(A)
  Y_negi = Y[-i]
  X_negi = X[-i, , drop=F]
  A_negi = A[-i]

  Y1tilde_negi <- as.numeric(A_negi * Y_negi / pi1_hat)
  Y0tilde_negi <- as.numeric((1 - A_negi) * Y_negi/ (1 - pi1_hat))

  ## Treatment arm
  # Step 1: check the existence of solution.
  l1_mu_fun <- function(mu){
    term1 <- t(X_negi) %*% Y1tilde_negi / (n - 1)
    term2 <- t(X_sub) %*% mu / nrow(X_sub)
    score_diff <- term1 - term2
    fval <- sum(score_diff**2 * 10^6) / (2 * length(score_diff))
    return(fval)
  }
  grad_l1_mu_fun <- function(mu){
    grad <- - 1 / nrow(X_sub) * X_sub %*% (t(X_negi) %*% Y1tilde_negi / (n - 1) - t(X_sub) %*% mu / nrow(X_sub))
    grad <- grad / length(grad) * 10^6
    return(grad)
  }

  # X_sub <- X_negi[A_negi == 1, ] ## (almost) standard leave one out：
  lb <- 0
  ub <- 1
  loc0 <- which(A == 0)
  n0 <- length(loc0)
  loc0_sub <- setdiff(loc0, i)
  X_sub <- rbind(X[i, ], X_negi[A_negi == 1, ], X[loc0_sub, ])

  mu_init <- as.numeric(1 / (1 + exp(-X_sub %*% beta1_hat)))
  mu1.sol <- optim(
    par = mu_init,
    fn = l1_mu_fun,
    gr = grad_l1_mu_fun,
    method = "L-BFGS-B",
    lower = lb,
    upper = ub,
    control = list(
      maxit = 5000,
      factr = 1e-10 / .Machine$double.eps  # 设置相对精度容差
    )
  )
  l1_value <- abs(mu1.sol$value)
  y1_hat_i <- mu1.sol$par[1]


  if(opt_obj == 'beta'){

    # Step 2: Use BBsolve if there exists possible solutions
    is_converge <- FALSE
    if(l1_value < 1e-6){
      est_fun1 <- function (beta) {
        term1 <- as.numeric(t(X_negi) %*% Y1tilde_negi / (n - 1))
        u <- as.numeric(1 / (1 + exp(-X_sub %*% beta)))
        term2 <- as.numeric(t(X_sub) %*% u / nrow(X_sub))
        return(term1 - term2)
      }
      beta1.sol <- BBsolve(par = rep(0, ncol(X_negi)), fn = est_fun1, quiet = TRUE, control = list(tol=1e-7))
      beta1_hat_negi <- beta1.sol$par
      y1_hat_sub <- as.numeric(1 / (1 + exp(-X_sub %*% beta1_hat_negi)))
      l1_value <- l1_mu_fun(y1_hat_sub)
      y1_hat_i <- y1_hat_sub[1]
      # est_sub_fun1(beta1.sol$par)
      is_converge <- ifelse(beta1.sol$convergence == 0, TRUE, FALSE)
    }

    # Step 3: If no solution exist or BBsolve does not converge, use optim
    if(!is_converge){
      l1_fun <- function(beta){
        term1 <- as.numeric(t(X_negi) %*% Y1tilde_negi / (n - 1))
        u <- as.numeric(1 / (1 + exp(-X_sub %*% beta)))
        term2 <- as.numeric(t(X_sub) %*% u / nrow(X_sub))
        score_diff <- term1 - term2
        fval <- sum(score_diff**2 * 10^6) / (2 * length(score_diff))
        return(fval)
      }
      grad_l1_fun <- function(beta){
        u <- as.numeric(1 / (1 + exp(-X_sub %*% beta)))
        grad <- - 1 / nrow(X_sub) * t(X_sub) %*% diag(u * (1 - u)) %*% X_sub %*% (t(X_negi) %*% Y1tilde_negi / (n - 1) - t(X_sub) %*% u / nrow(X_sub))
        grad <- grad / length(grad) * 10^6
        return(grad)
      }

      lb <- min(beta1_hat)
      ub <- max(beta1_hat)

      beta1.sol <- optim(
        par = beta1_hat,
        fn = l1_fun,
        gr = grad_l1_fun,
        method = "L-BFGS-B",
        lower = lb,
        upper = ub,
        control = list(
          maxit = 5000,
          factr = 1e-10 / .Machine$double.eps  # 设置相对精度容差
        )
      )
      beta1_hat_negi <- beta1.sol$par
      l1_value <- beta1.sol$value
      y1_hat_i <- as.numeric(1 / (1 + exp(- X[i, ] %*% beta1_hat_negi)))
    }

  }


  ## Control arm
  ## reduced form：
  l0_mu_fun <- function(mu){
    term1 <- t(X_negi) %*% Y0tilde_negi / (n - 1)
    term2 <- t(X_sub) %*% mu / nrow(X_sub)
    score_diff <- term1 - term2
    fval <- sum(score_diff**2 * 10^6) / (2 * length(score_diff))
    return(fval)
  }
  grad_l0_mu_fun <- function(mu){
    grad <- - 1 / nrow(X_sub) * X_sub %*% (t(X_negi) %*% Y0tilde_negi / (n - 1) - t(X_sub) %*% mu / nrow(X_sub))
    grad <- grad / length(grad) * 10^6
    return(grad)
  }
  lb <- 0
  ub <- 1
  loc1 <- which(A == 1)
  n1 <- length(loc1)

  loc1_sub <- setdiff(loc1, i)
  X_sub <- rbind(X[i, ], X_negi[A_negi == 0, ], X[loc1_sub, ])

  mu_init <- as.numeric(1 / (1 + exp(-X_sub %*% beta0_hat)))
  mu0.sol <- optim(
    par = mu_init,
    fn = l0_mu_fun,
    gr = grad_l0_mu_fun,
    method = "L-BFGS-B",
    lower = lb,
    upper = ub,
    control = list(
      maxit = 5000,
      factr = 1e-10 / .Machine$double.eps
    )
  )
  l0_value <- abs(mu0.sol$value)
  y0_hat_i <- mu0.sol$par[1]

  if(opt_obj == 'beta'){
    is_converge <- FALSE
    if(l0_value < 1e-6){
      est_fun0 <- function (beta) {
        term1 <- as.numeric(t(X_negi) %*% Y0tilde_negi / (n - 1))
        u <- as.numeric(1 / (1 + exp(-X_sub %*% beta)))
        term2 <- as.numeric(t(X_sub) %*% u / nrow(X_sub))
        return(term1 - term2)
      }
      beta0.sol <- BBsolve(par = rep(0, ncol(X_negi)), fn = est_fun0, quiet = TRUE, control = list(tol=1e-7))
      beta0_hat_negi <- beta0.sol$par
      y0_hat_sub <- as.numeric(1 / (1 + exp(-X_sub %*% beta0_hat_negi)))
      l0_value <- l0_mu_fun(y0_hat_sub)
      y0_hat_i <- y0_hat_sub[1]
      is_converge <- ifelse(beta0.sol$convergence == 0, TRUE, FALSE)
    }

    if(!is_converge){
      l0_fun <- function(beta){
        term1 <- as.numeric(t(X_negi) %*% Y0tilde_negi / (n - 1))
        u <- as.numeric(1 / (1 + exp(-X_sub %*% beta)))
        term2 <- as.numeric(t(X_sub) %*% u / nrow(X_sub))
        score_diff <- term1 - term2
        fval <- sum(score_diff**2 * 10^6) / (2 * length(score_diff))
        return(fval)
      }
      grad_l0_fun <- function(beta){
        u <- as.numeric(1 / (1 + exp(-X_sub %*% beta)))
        grad <- - 1 / nrow(X_sub) * t(X_sub) %*% diag(u * (1 - u)) %*% X_sub %*% (t(X_negi) %*% Y0tilde_negi / (n - 1) - t(X_sub) %*% u / nrow(X_sub))
        grad <- grad / length(grad) * 10^6
        return(grad)
      }

      lb <- min(beta0_hat)
      ub <- max(beta0_hat)

      beta0.sol <- optim(
        par = beta0_hat,
        fn = l0_fun,
        gr = grad_l0_fun,
        method = "L-BFGS-B",
        lower = lb,
        upper = ub,
        control = list(
          maxit = 5000,
          factr = 1e-10 / .Machine$double.eps  # 设置相对精度容差
        )
      )
      beta0_hat_negi <- beta0.sol$par
      y0_hat_sub <- as.numeric(1 / (1 + exp(-X_sub %*% beta0_hat_negi)))
      l0_value <- l0_mu_fun(y0_hat_sub)
      # l0_fun(beta0_hat_negi)
      y0_hat_i <- y0_hat_sub[1]
    }
  }

  yhat_i_negi <- c(y0_hat_i, y1_hat_i)
  l0_l1_value = c(l0_value, l1_value)
  names(l0_l1_value) <- c('l0', 'l1')

  return(list(
    yhat_i_negi = yhat_i_negi,
    l0_l1_value = l0_l1_value
  ))
}

#' @keywords internal
#' @import stats
#' @import BB
fit.JASA.i.gaussian <- function(Y, X, A, pi1_hat, beta1_hat, beta0_hat, i, opt_obj){
  # require(stats); require(BB)
  n <- length(A)
  Y_negi = Y[-i]
  X_negi = X[-i, , drop=F]
  A_negi = A[-i]
  X2 <- t(X) %*% X

  Y1tilde_negi <- as.numeric(A_negi * Y_negi / pi1_hat)
  Y0tilde_negi <- as.numeric((1 - A_negi) * Y_negi/ (1 - pi1_hat))

  ## Treatment arm
  ## Step 1:
  l1_mu_fun <- function(mu){
    term1 <- t(X_negi) %*% Y1tilde_negi / (n - 1)
    term2 <- t(X) %*% mu / n
    score_diff <- term1 - term2
    fval <- sum(score_diff**2 * 10^6) / (2 * length(score_diff))
    return(fval)
  }
  grad_l1_mu_fun <- function(mu){
    grad <- - 1 / n * X %*% (t(X_negi) %*% Y1tilde_negi / (n - 1) - t(X) %*% mu / n)
    grad <- grad / length(grad) * 10^6
    return(grad)
  }
  mu_init <- as.numeric(X %*% beta1_hat)
  lb <- rep(-Inf, ncol(X))
  ub <- rep(Inf, ncol(X))
  mu1.sol <- optim(
    par = mu_init,
    fn = l1_mu_fun,
    gr = grad_l1_mu_fun,
    method = "L-BFGS-B",
    lower = lb,
    upper = ub,
    control = list(
      maxit = 5000,
      factr = 1e-10 / .Machine$double.eps
    )
  )
  y1_hat_i <- mu1.sol$par[i]
  l1_value <- mu1.sol$value


  if(opt_obj == 'beta'){
    is_converge <- FALSE
    # Step 2:
    if(l1_value < 1e-6){
      est_fun1 <- function (beta) {
        term1 <- as.numeric(t(X_negi) %*% Y1tilde_negi / (n - 1))
        term2 <- as.numeric(X2 %*% beta / n)
        return(term1 - term2)
      }
      beta1.sol <- BBsolve(par = rep(0, ncol(X_negi)), fn = est_fun1, quiet = TRUE, control = list(tol=1e-7))
      beta1_hat_negi <- beta1.sol$par
      y1_hat <- as.numeric(X %*% beta1_hat_negi)
      l1_value <- l1_mu_fun(y1_hat)
      y1_hat_i <- y1_hat[i]
      is_converge <- ifelse(beta1.sol$convergence == 0, TRUE, FALSE)
    }

    # Step 3
    if(!is_converge){
      l1_fun <- function(beta){
        term1 <- as.numeric(t(X_negi) %*% Y1tilde_negi / (n - 1))
        term2 <- as.numeric(X2 %*% beta / n)
        score_diff <- term1 - term2
        fval <- sum(score_diff**2 * 10^6) / (2 * length(score_diff))
        return(fval)
      }
      grad_l1_fun <- function(beta){
        u <- as.numeric(X %*% beta)
        grad <- - 1 / n * X2 %*% (t(X_negi) %*% Y1tilde_negi / (n - 1) - t(X) %*% u / n)
        grad <- grad / length(grad) * 10^6
        return(grad)
      }

      lb <- min(beta1_hat)
      ub <- max(beta1_hat)

      beta1.sol <- optim(
        par = beta1_hat,
        fn = l1_fun,
        gr = grad_l1_fun,
        method = "L-BFGS-B",
        lower = lb,
        upper = ub,
        control = list(
          maxit = 5000,
          factr = 1e-10 / .Machine$double.eps  # 设置相对精度容差
        )
      )
      beta1_hat_negi <- beta1.sol$par
      l1_value <- beta1.sol$value
      y1_hat_i <- as.numeric(X[i, ] %*% beta1_hat_negi)
    }
  }


  ## Control arm
  # Step 1
  l0_mu_fun <- function(mu){
    term1 <- t(X_negi) %*% Y0tilde_negi / (n - 1)
    term2 <- t(X) %*% mu / n
    score_diff <- term1 - term2
    fval <- sum(score_diff**2 * 10^6) / (2 * length(score_diff))
    return(fval)
  }
  grad_l0_mu_fun <- function(mu){
    grad <- - 1 / n * X %*% (t(X_negi) %*% Y0tilde_negi / (n - 1) - t(X) %*% mu / n)
    grad <- grad / length(grad) * 10^6
    return(grad)
  }
  mu_init <- as.numeric(X %*% beta0_hat)
  lb <- rep(-Inf, ncol(X))
  ub <- rep(Inf, ncol(X))
  mu0.sol <- optim(
    par = mu_init,
    fn = l0_mu_fun,
    gr = grad_l0_mu_fun,
    method = "L-BFGS-B",
    lower = lb,
    upper = ub,
    control = list(
      maxit = 5000,
      factr = 1e-10 / .Machine$double.eps  # 设置相对精度容差
    )
  )

  l0_value <- mu0.sol$value
  y0_hat_i <- mu0.sol$par[i]

  if(opt_obj == 'beta'){
    is_converge <- FALSE
    # Step 2
    if(l0_value < 1e-6){
      est_fun0 <- function (beta) {
        term1 <- as.numeric(t(X_negi) %*% Y0tilde_negi / (n - 1))
        term2 <- as.numeric(X2 %*% beta / n)
        return(term1 - term2)
      }
      beta0.sol <- BBsolve(par = rep(0, ncol(X_negi)), fn = est_fun0, quiet = TRUE, control = list(tol=1e-7))
      beta0_hat_negi <- beta0.sol$par
      y0_hat <- as.numeric(X %*% beta0_hat_negi)
      l0_value <- l0_mu_fun(y0_hat)
      y0_hat_i <- y0_hat[i]
      is_converge <- ifelse(beta0.sol$convergence == 0, TRUE, FALSE)
    }

    # Step 3:
    if(!is_converge){
      l0_fun <- function(beta){
        term1 <- as.numeric(t(X_negi) %*% Y0tilde_negi / (n - 1))
        term2 <- as.numeric(X2 %*% beta / n)
        score_diff <- term1 - term2
        fval <- sum(score_diff**2 * 10^6) / (2 * length(score_diff))
        return(fval)
      }
      grad_l0_fun <- function(beta){
        u <- as.numeric(X %*% beta)
        grad <- - 1 / n * X2 %*% (t(X_negi) %*% Y0tilde_negi / (n - 1) - t(X) %*% u / n)
        grad <- grad / length(grad) * 10^6
        return(grad)
      }

      lb <- min(beta0_hat)
      ub <- max(beta0_hat)

      beta0.sol <- optim(
        par = beta0_hat,
        fn = l0_fun,
        gr = grad_l0_fun,
        method = "L-BFGS-B",
        lower = lb,
        upper = ub,
        control = list(
          maxit = 5000,
          factr = 1e-10 / .Machine$double.eps  # 设置相对精度容差
        )
      )
      beta0_hat_negi <- beta0.sol$par
      l0_value <- beta0.sol$value
      y0_hat_i <- as.numeric(X[i, ] %*% beta0_hat_negi)
    }
  }

  yhat_i_negi <- c(y0_hat_i, y1_hat_i)
  l0_l1_value = c(l0_value, l1_value)
  names(l0_l1_value) <- c('l0', 'l1')

  return(list(
    yhat_i_negi = yhat_i_negi,
    l0_l1_value = l0_l1_value
  ))

}

#' @keywords internal
#' @import stats
#' @import BB
fit.JASA.i.poisson <- function(Y, X, A, pi1_hat, beta1_hat, beta0_hat, i, opt_obj){
  # require(stats); require(BB)
  n <- length(A)
  Y_negi = Y[-i]
  X_negi = X[-i, , drop=F]
  A_negi = A[-i]

  Y1tilde_negi <- as.numeric(A_negi * Y_negi / pi1_hat)
  Y0tilde_negi <- as.numeric((1 - A_negi) * Y_negi/ (1 - pi1_hat))

  ## Treatment arm
  # Step 1
  l1_mu_fun <- function(mu){
    term1 <- t(X_negi) %*% Y1tilde_negi / (n - 1)
    term2 <- t(X_sub) %*% mu / nrow(X_sub)
    score_diff <- term1 - term2
    fval <- sum(score_diff**2) / (2 * length(score_diff))
    return(fval)
  }
  grad_l1_mu_fun <- function(mu){
    grad <- - 1 / nrow(X_sub) * X_sub %*% (t(X_negi) %*% Y1tilde_negi / (n - 1) - t(X_sub) %*% mu / nrow(X_sub))
    grad <- grad / length(grad)
    return(grad)
  }

  lb <- 0
  ub <- Inf
  loc0 <- which(A == 0)
  n0 <- length(loc0)

  loc0_sub <- setdiff(loc0, i)
  X_sub <- rbind(X[i, ], X_negi[A_negi == 1, ], X[loc0_sub, ])

  mu_init <- rep(1, nrow(X_sub))
  mu1.sol <- optim(
    par = mu_init,
    fn = l1_mu_fun,
    gr = grad_l1_mu_fun,
    method = "L-BFGS-B",
    lower = lb,
    upper = ub,
    control = list(
      maxit = 5000,
      factr = 1e-10 / .Machine$double.eps
    )
  )
  l1_value <- abs(mu1.sol$value)
  y1_hat_i <- mu1.sol$par[1]


  if(opt_obj == 'beta'){
    # Step 2:
    is_converge <- FALSE
    if(l1_value < 1e-6){
      est_fun1 <- function (beta) {
        term1 <- as.numeric(t(X_negi) %*% Y1tilde_negi / (n - 1))
        lp <- pmin(as.numeric(X_sub %*% beta), 700)
        term2 <- as.numeric(t(X_sub) %*% exp(lp) / nrow(X_sub))
        return(term1 - term2)
      }
      beta1.sol <- BBsolve(par = rep(0,ncol(X_negi)), fn = est_fun1, quiet = TRUE, control = list(tol=1e-7))
      beta1_hat_negi <- beta1.sol$par
      is_converge <- ifelse(beta1.sol$convergence == 0, TRUE, FALSE)
      y1_hat_sub <- as.numeric(exp(X_sub %*% beta1_hat_negi))
      l1_value <- l1_mu_fun(y1_hat_sub)
      y1_hat_i <- y1_hat_sub[1]
    }

    # Step 3:
    if(!is_converge){
      l1_fun <- function(beta){
        term1 <- as.numeric(t(X_negi) %*% Y1tilde_negi / (n - 1))
        lp <- pmin(as.numeric(X_sub %*% beta), 700)
        term2 <- as.numeric(t(X_sub) %*% exp(lp) / nrow(X_sub))
        score_diff <- term1 - term2
        fval <- sum(score_diff**2) / (2 * length(score_diff))
        return(fval)
      }
      grad_l1_fun <- function(beta){
        lp <- as.numeric(X_sub %*% beta)
        lp <- pmin(lp, 700)
        u <- as.numeric(exp(lp))
        grad <- - 1 / nrow(X_sub) * t(X_sub) %*% diag(u) %*% X_sub %*% (t(X_negi) %*% Y1tilde_negi / (n - 1) - t(X_sub) %*% u / nrow(X_sub))
        grad <- grad / length(grad)
        return(grad)
      }


      lb <- min(beta1_hat)
      ub <- max(beta1_hat)

      beta1.sol <- optim(
        par = beta1_hat,
        fn = l1_fun,
        gr = grad_l1_fun,
        method = "L-BFGS-B",
        lower = lb,
        upper = ub,
        control = list(
          maxit = 5000,
          factr = 1e-10 / .Machine$double.eps  # 设置相对精度容差
        )
      )
      beta1_hat_negi <- beta1.sol$par
      l1_value <- beta1.sol$value
      y1_hat_i <- as.numeric(exp(X[i, ] %*% beta1_hat_negi))
    }
  }


  ## Control arm
  # Step 1:
  l0_mu_fun <- function(mu){
    term1 <- t(X_negi) %*% Y0tilde_negi / (n - 1)
    term2 <- t(X_sub) %*% mu / nrow(X_sub)
    score_diff <- term1 - term2
    fval <- sum(score_diff**2) / (2 * length(score_diff))
    return(fval)
  }
  grad_l0_mu_fun <- function(mu){
    grad <- - 1 / nrow(X_sub) * X_sub %*% (t(X_negi) %*% Y0tilde_negi / (n - 1) - t(X_sub) %*% mu / nrow(X_sub))
    grad <- grad / length(grad)
    return(grad)
  }

  lb <- 0
  ub <- Inf
  loc1 <- which(A == 1)
  n1 <- length(loc1)
  loc1_sub <- setdiff(loc1, i)
  X_sub <- rbind(X[i, ], X_negi[A_negi == 0, ], X[loc1_sub, ])

  mu_init <- rep(1, nrow(X))
  mu0.sol <- optim(
    par = mu_init,
    fn = l0_mu_fun,
    gr = grad_l0_mu_fun,
    method = "L-BFGS-B",
    lower = lb,
    upper = ub,
    control = list(
      maxit = 5000,
      factr = 1e-10 / .Machine$double.eps  # 设置相对精度容差
    )
  )
  l0_value <- abs(mu0.sol$value)
  y0_hat_i <- mu0.sol$par[1]


  if(opt_obj == 'beta'){
    is_converge <- FALSE
    if(l0_value < 1e-6){
      est_fun0 <- function (beta) {
        term1 <- as.numeric(t(X_negi) %*% Y0tilde_negi / (n - 1))
        lp <- pmin(as.numeric(X_sub %*% beta), 700)
        term2 <- as.numeric(t(X_sub) %*% exp(lp) / nrow(X_sub))
        return(term1 - term2)
      }
      beta0.sol <- BBsolve(par = rep(0, ncol(X_negi)), fn = est_fun0, quiet = TRUE, control = list(tol=1e-7))
      beta0_hat_negi <- beta0.sol$par
      is_converge <- ifelse(beta0.sol$convergence == 0, TRUE, FALSE)
      y0_hat_sub <- as.numeric(exp(X_sub %*% beta0_hat_negi))
      l0_value <- l0_mu_fun(y0_hat_sub)
      y0_hat_i <- y0_hat_sub[1]
    }

    # Step 3
    if(!is_converge){
      l0_fun <- function(beta){
        term1 <- as.numeric(t(X_negi) %*% Y0tilde_negi / (n - 1))
        lp <- pmin(as.numeric(X_sub %*% beta), 700)
        term2 <- as.numeric(t(X_sub) %*% exp(lp) / nrow(X_sub))
        score_diff <- term1 - term2
        fval <- sum(score_diff**2) / (2 * length(score_diff))
        return(fval)
      }
      grad_l0_fun <- function(beta){
        lp <- as.numeric(X_sub %*% beta)
        lp <- pmin(lp, 700)
        u <- as.numeric(exp(lp))
        grad <- - 1 / nrow(X_sub) * t(X_sub) %*% diag(u) %*% X_sub %*% (t(X_negi) %*% Y0tilde_negi / (n - 1) - t(X_sub) %*% u / nrow(X_sub))
        grad <- grad / length(grad)
        return(grad)
      }

      lb <- min(beta0_hat)
      ub <- max(beta0_hat)

      beta0.sol <- optim(
        par = beta0_hat,
        fn = l0_fun,
        gr = grad_l0_fun,
        method = "L-BFGS-B",
        lower = lb,
        upper = ub,
        control = list(
          maxit = 5000,
          factr = 1e-10 / .Machine$double.eps
        )
      )
      beta0_hat_negi <- beta0.sol$par
      y0_hat_sub <- as.numeric(exp(X_sub %*% beta0_hat_negi))
      l0_value <- l0_fun(beta0_hat_negi)
      y0_hat_i <- y0_hat_sub[1]
    }

  }

  yhat_i_negi <- c(y0_hat_i, y1_hat_i)
  l0_l1_value = c(l0_value, l1_value)
  names(l0_l1_value) <- c('l0', 'l1')


  return(list(
    yhat_i_negi = yhat_i_negi,
    l0_l1_value = l0_l1_value
  ))
}

