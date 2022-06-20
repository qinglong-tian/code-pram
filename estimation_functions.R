# Fit X*~Y: for X is binary outcome
fit_x_star_y_logistic <- function(data)
{
  fit <- glm(X_star ~ Y, family = "binomial", data = data)
  
  return(fit)
}

# Method 1: \tilde{E}(X|Y) is the same as \tilde{E}(X^\ast|Y)
Tilde_E_X_Y_Naive <- function(data, p00, p11)
{
  glmFit <- fit_x_star_y_logistic(data)
  glmFit$fitted.values
}

# Method 2: Solving the linear equations
Tilde_E_X_Y_Solve <- function(data, p00, p11)
{
  glmFit <- fit_x_star_y_logistic(data)
  p10 <- 1 - p00
  fittedValNaive <- glmFit$fitted.values
  fittedValSolve <- (fittedValNaive - p10) / (p11 - p10)
  fittedValSolve[fittedValSolve < 0] <- 0
  fittedValSolve[fittedValSolve > 1] <- 1
  
  return(fittedValSolve)
}

# Compute EXP(-Odd)
Compute_Exp_M_Odd <- function(beta_hat, yVec)
{
  exp(-c(cbind(1, yVec) %*% matrix(beta_hat, ncol = 1)))
}

# Compute Gamma
Compute_Gamma <- function(beta_hat, yVec)
{
  exp_m_odd <- Compute_Exp_M_Odd(beta_hat, yVec)
  exp_m_odd / (1 + exp_m_odd) ^ 2
}

# Compute Omega
Compute_Omega <- function(beta_hat, yVec)
{
  gammaVec <- Compute_Gamma(beta_hat, yVec)
  e_y2_gamma <- mean(yVec ^ 2 * gammaVec)
  e_y_gamma <- mean(yVec * gammaVec)
  e_gamma <- mean(gammaVec)
  
  e_dU_dbeta <- matrix(nrow = 2, ncol = 2)
  e_dU_dbeta[1, 1] <- e_gamma
  e_dU_dbeta[1, 2] <- e_y_gamma
  e_dU_dbeta[2, 1] <- e_y_gamma
  e_dU_dbeta[2, 2] <- e_y2_gamma
  e_dU_dbeta <- -e_dU_dbeta
  solve(e_dU_dbeta)
}

Compute_ABCD <- function(beta_hat, yVec, tilde_e)
{
  omega <- Compute_Omega(beta_hat, yVec)
  bb1 <- cbind(1 - tilde_e, yVec * (1 - tilde_e))
  AC <- t(omega %*% t(bb1))
  
  bb2 <- cbind(-tilde_e,-yVec * tilde_e)
  BD <- t(omega %*% t(bb2))
  
  return(list(AC = AC, BD = BD))
}

Compute_C_Func <- function(ABCD_list, p00, p11)
{
  AC <- ABCD_list$AC
  BD <- ABCD_list$BD
  
  p01 <- 1 - p11
  p10 <- 1 - p00
  
  c11 <- (p00 * AC[, 1] - p01 * BD[, 1]) / (p11 * p00 - p01 * p10)
  c1m1 <- (p11 * BD[, 1] - p10 * AC[, 1]) / (p11 * p00 - p01 * p10)
  
  c21 <- (p00 * AC[, 2] - p01 * BD[, 2]) / (p11 * p00 - p01 * p10)
  c2m1 <- (p11 * BD[, 2] - p10 * AC[, 2]) / (p11 * p00 - p01 * p10)
  
  c_x_ast_1 <- cbind(c11, c21)
  c_x_ast_m1 <- cbind(c1m1, c2m1)
  
  return(list(c_func_x_ast_1 = c_x_ast_1, c_func_x_ast_m1 = c_x_ast_m1))
}

# Compute Efficient IF
Compute_Efficient_IF_Logistic <- function(beta_hat, dat, p00, p11, tilde_Func)
{
  yVec <- dat$Y
  e_tilde <- tilde_Func(dat, p00, p11)
  ABCD_List <- Compute_ABCD(beta_hat, yVec, e_tilde)
  c_func <- Compute_C_Func(ABCD_List, p00, p11)
  
  C_1 <- c_func$c_func_x_ast_1
  C_m1 <- c_func$c_func_x_ast_m1
  
  X_star <- dat$X_star
  
  IFVec <- matrix(nrow = length(yVec), ncol = 2)
  for (i in 1:length(yVec))
  {
    if (X_star[i] == 0)
    {
      IFVec[i,] <- C_m1[i,]
    }
    else
    {
      IFVec[i,] <- C_1[i,]
    }
  }
  
  tilde_e_u_y_1 <- e_tilde-1/(1+Compute_Exp_M_Odd(beta_hat, yVec))
  tilde_e_u_y_2 <- yVec*e_tilde-yVec/(1+Compute_Exp_M_Odd(beta_hat, yVec))
  omegaMat <- Compute_Omega(beta_hat, yVec)
  
  added_term <- omegaMat %*% t(cbind(tilde_e_u_y_1, tilde_e_u_y_2))
  added_term <- t(added_term)
  
  IFVec <- IFVec+added_term
  
  return(IFVec)
}

Compute_Efficient_IF_Logistic_Sum <- function(beta_hat,
                                              dat,
                                              p00,
                                              p11,
                                              tilde_Func)
{
  Compute_Efficient_IF_Logistic(beta_hat,
                                dat,
                                p00,
                                p11,
                                tilde_Func) -> mat
  sum(colMeans(mat) ^ 2)
}

# Compute Direct-Solve Method
Compute_Beta_Direct_Solve <- function(data, p00, p11)
{
  Tilde_E_X_Y_Solve(data, p00, p11) -> tilde_e_x_y
  to_remove <- which(tilde_e_x_y %in% c(0, 1))
  if (length(to_remove) == 0)
  {
    tilde_e_x_y_clean <- tilde_e_x_y
    y_clean <- data$Y
  }
  else
  {
    tilde_e_x_y_clean <- tilde_e_x_y[-to_remove]
    y_clean <- data$Y[-to_remove]
  }
  
  odd_clean <- log(tilde_e_x_y_clean / (1 - tilde_e_x_y_clean))
  newDat <- data.frame(resp = odd_clean, expl = y_clean)
  lm(resp ~ expl, data = newDat)$coefficients
}
