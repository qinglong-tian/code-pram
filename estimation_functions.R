# Fit X*~Y: for X is binary outcome
fit_x_star_y_logistic <- function(data)
{
  fit <- glm(X_star ~ Y, family = "binomial", data = data)
  
  return(fit)
}

# Predict Binary $X$
E_X_Y <- function(yVec, betaVecY)
{
  oddVec <- cbind(1, yVec) %*% matrix(betaVecY, ncol = 1)
  probVec <- 1 / (1 + exp(-oddVec))
  
  as.numeric(probVec > 0.5)
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

# Compute A - Logistic
Compute_A <- function(beta_vec, tilde_e_x_y, yVec)
{
  oddVec <- cbind(1, yVec) %*% matrix(beta_vec, ncol = 1)
  gammaVec <- 1 / (1 + exp(-oddVec))
  (tilde_e_x_y - 1) / (gammaVec ^ 2 - gammaVec) / var(yVec) * (yVec * mean(yVec) -
                                                                 mean(yVec ^ 2))
}

# Compute C - Logistic
Compute_C <- function(beta_vec, tilde_e_x_y, yVec)
{
  oddVec <- cbind(1, yVec) %*% matrix(beta_vec, ncol = 1)
  gammaVec <- 1 / (1 + exp(-oddVec))
  (tilde_e_x_y - 1) / (gammaVec ^ 2 - gammaVec) / var(yVec) * (mean(yVec) -
                                                                 yVec)
}

# Compute B - Logistic
Compute_B <- function(beta_vec, tilde_e_x_y, yVec)
{
  oddVec <- cbind(1, yVec) %*% matrix(beta_vec, ncol = 1)
  gammaVec <- 1 / (1 + exp(-oddVec))
  tilde_e_x_y / (gammaVec ^ 2 - gammaVec) / var(yVec) * (yVec * mean(yVec) -
                                                           mean(yVec ^ 2))
}

# Compute D - Logistic
Compute_D <- function(beta_vec, tilde_e_x_y, yVec)
{
  oddVec <- cbind(1, yVec) %*% matrix(beta_vec, ncol = 1)
  gammaVec <- 1 / (1 + exp(-oddVec))
  tilde_e_x_y / (gammaVec ^ 2 - gammaVec) / var(yVec) * (mean(yVec) - yVec)
}

# Compute \tilde{c}: [c(x^\ast=1, y), c(x^\ast=0, y)]
Compute_Tilde_C <- function(beta_vec, tilde_e_x_y, yVec, p00, p11)
{
  AVec <- Compute_A(beta_vec, tilde_e_x_y, yVec)
  BVec <- Compute_B(beta_vec, tilde_e_x_y, yVec)
  CVec <- Compute_C(beta_vec, tilde_e_x_y, yVec)
  DVec <- Compute_D(beta_vec, tilde_e_x_y, yVec)
  
  p01 <- 1 - p11
  p10 <- 1 - p00
  
  der <- (p11 * p00 - p01 * p10)
  
  c_x1_1 <- (p00 * AVec - p01 * BVec) / der
  c_x1_2 <- (p00 * CVec - p01 * DVec) / der
  
  c_x0_1 <- (p11 * BVec - p10 * AVec) / der
  c_x0_2 <- (p11 * DVec - p10 * CVec) / der
  
  cbind(c_x1_1, c_x1_2, c_x0_1, c_x0_2)
}

# Compute \Omega*\tilde{E}(U|Y)
Compute_Added_Term <- function(beta_vec, yVec, tilde_e_x_y)
{
  oddVec <- cbind(1, yVec) %*% matrix(beta_vec, ncol = 1)
  gammaVec <- 1 / (1 + exp(-oddVec))
  
  (tilde_e_x_y - gammaVec) / (gammaVec ^ 2 - gammaVec) / var(yVec) -> firstTerm
  
  secondTerm1 <- mean(yVec ^ 2) - yVec * mean(yVec)
  secondTerm2 <- yVec - mean(yVec)
  
  return(cbind(firstTerm * secondTerm1, firstTerm * secondTerm2))
}

# Compute Efficient IF
Compute_Efficient_IF_Logistic <-
  function(beta_hat,
           dat_logistic,
           p00,
           p11,
           tilde_Func)
  {
    X_star <- dat_logistic$X_star
    yVec <- dat_logistic$Y
    # Flexible
    flex_if <- matrix(0, ncol = 2, nrow = length(yVec))
    
    tilde_e_x_y_naive <- tilde_Func(dat_logistic, p00, p11)
    Compute_Tilde_C(beta_hat, tilde_e_x_y_naive, yVec, p00, p11) -> tilde_c_flex
    Compute_Added_Term(beta_hat, yVec, tilde_e_x_y_naive) -> added_flex
    
    for (i in 1:length(yVec))
    {
      flex_if[i,] <-
        ifelse(X_star[i] == 1, tilde_c_flex[i, c(1, 2)], tilde_c_flex[i, c(3, 4)])
    }
    
    flex_if - added_flex
  }

Compute_Efficient_IF_Logistic_Sum <- function(beta_hat,
                                              dat_logistic,
                                              p00,
                                              p11,
                                              tilde_Func)
{
  Compute_Efficient_IF_Logistic(beta_hat,
                                dat_logistic,
                                p00,
                                p11,
                                tilde_Func) -> mat
  sum(colMeans(mat)^2)
}
