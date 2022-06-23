# Fit X*~Y: for X is binary outcome
fit_x_star_y_logistic <- function(data)
{
  fit <- glm(X_star ~ Y, family = "binomial", data = data)
  
  return(fit)
}

# Make Mat Q (inverse of P)
Compute_Q <- function(p00, p11)
{
  PMat <- matrix(ncol = 2, nrow = 2)
  PMat[1, 1] <- p00
  PMat[1, 2] <- 1 - p00
  PMat[2, 1] <- 1 - p11
  PMat[2, 2] <- p11
  
  return(solve(PMat))
}

# Compute Exp(-beta_0-beta_1*y)
Compute_Exp <- function(beta_hat, yVec)
{
  oddVec <- cbind(1, yVec) %*% matrix(beta_hat, ncol = 1)
  return(exp(-oddVec))
}

# Compute Prob in Logistic
Compute_Expit <- function(beta_hat, yVec)
{
  exp_ <- Compute_Exp(beta_hat, yVec)
  return(1 / (1 + exp_))
}

# Compute Omega Logistic

Compute_Gamma <- function(beta_hat, yVec)
{
  exp_ <- Compute_Exp(beta_hat, yVec)
  return(exp_ / (1 + exp_) ^ 2)
}

Compute_Omega_Logistic <- function(beta_hat, yVec)
{
  gammaVec <- Compute_Gamma(beta_hat, yVec)
  Omega <- matrix(nrow = 2, ncol = 2)
  
  # Assign Values
  Omega[1, 1] <- mean(gammaVec)
  Omega[1, 2] <- mean(yVec * gammaVec)
  Omega[2, 1] <- Omega[1, 2]
  Omega[2, 2] <- mean(yVec ^ 2 * gammaVec)
  
  return(-solve(Omega))
}

# Compute U Function Logistic
Compute_U_Logistic <- function(beta_hat, X, yVal)
{
  expit <- Compute_Expit(beta_hat, yVal)
  return(c(X - expit, (X - expit) * yVal))
}

# Compute Mat B in Logistic
Compute_B_Logistic <- function(beta_hat, yVal, Omg_Log)
{
  B <- matrix(ncol = 2, nrow = 2)
  
  Ux0 <- Compute_U_Logistic(beta_hat, 0, yVal)
  Ux1 <- Compute_U_Logistic(beta_hat, 1, yVal)
  
  B[, 1] <- Omg_Log %*% matrix(Ux0, ncol = 1)
  B[, 2] <- Omg_Log %*% matrix(Ux1, ncol = 1)
  
  return(B)
}

Compute_C_Logistic <- function(beta_hat, yVal, Omg_Log, p00, p11)
{
  BMat <- Compute_B_Logistic(beta_hat, yVal, Omg_Log)
  QMat <- Compute_Q(p00, p11)
  
  return(QMat %*% t(BMat))
}

Compute_IF_Logistic <- function(beta_hat, yVec, xStarVec, p00, p11)
{
  Omg_Log <- Compute_Omega_Logistic(beta_hat, yVec)
  IF_Mat <- matrix(nrow = length(yVec), ncol = 2)
  for (i in 1:length(yVec))
  {
    yVal <- yVec[i]
    xStar <- xStarVec[i]
    CMat <- Compute_C_Logistic(beta_hat, yVal, Omg_Log, p00, p11)
    IF_Mat[i, ] <- CMat[xStar + 1, ]
  }
  return(IF_Mat)
}

Compute_IF_Logistic_ColSum <- function(beta_hat, yVec, xStarVec, p00, p11)
{
  IFMat <- Compute_IF_Logistic(beta_hat, yVec, xStarVec, p00, p11)
  sum(colMeans(IFMat)^2)
}
