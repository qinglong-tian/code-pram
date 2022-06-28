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

Compute_C_Logistic <- function(beta_hat, yVal, Omg_Log, p00, p11, xStarVec, use_q)
{
  BMat <- Compute_B_Logistic(beta_hat, yVal, Omg_Log)
  if(use_q)
  {
    QMat <- Compute_Q(p00, p11)
  }
  else
  {
    PX_ast <- c(1-mean(xStarVec), mean(xStarVec))
    QMat <- t(Compute_X_Given_X_ast(PX_ast, p00, p11))
  }
  
  return(QMat %*% t(BMat))
}

Compute_IF_Logistic <- function(beta_hat, data, p00, p11, no_omega, use_q)
{
  yVec <- data$Y
  xStarVec <- data$X_star
  if(no_omega)
  {
    Omg_Log <- diag(2)
  }
  else
  {
    Omg_Log <- Compute_Omega_Logistic(beta_hat, yVec)
  }
  IF_Mat <- matrix(nrow = length(yVec), ncol = 2)
  for (i in 1:length(yVec))
  {
    yVal <- yVec[i]
    xStar <- xStarVec[i]
    CMat <- Compute_C_Logistic(beta_hat, yVal, Omg_Log, p00, p11, xStarVec, use_q)
    IF_Mat[i,] <- CMat[xStar + 1,]
  }
  return(IF_Mat)
}

Compute_IF_Logistic_ColSum <-
  function(beta_hat, data, p00, p11, no_omega = F, use_q = T)
  {
    IFMat <- Compute_IF_Logistic(beta_hat, data, p00, p11, no_omega, use_q)
    sum(colMeans(IFMat) ^ 2)
  }

# Solving P(X|Y): no adjustment
Compute_P_Inv <- function(p00, p11)
{
  P <- matrix(ncol = 2, nrow = 2)
  P[1, 1] <- p00
  P[1, 2] <- 1 - p11
  P[2, 1] <- 1 - p00
  P[2, 2] <- p11
  
  return(solve(P))
}

Compute_P_X_star_Given_Y <-
  function(data,
           probit = F,
           ommit_intercept = F)
  {
    if (probit)
    {
      link_func <- "probit"
    }
    else
    {
      link_func <- "logit"
    }
    
    if (ommit_intercept)
    {
      fit_glm <-
        glm(X_star ~ 0 + Y,
            family = binomial(link = link_func),
            data = data)
    }
    else
    {
      fit_glm <-
        glm(X_star ~ Y,
            family = binomial(link = link_func),
            data = data)
    }
    return(fit_glm$fitted.values)
  }

Compute_IF_By_Solve_Method2 <- function(beta_hat, data, p00, p11, probit, ommit_intercept)
{
  prob_x_star_1_y <- Compute_P_X_star_Given_Y(data, probit, ommit_intercept)
  pinv <- Compute_P_Inv(p00, p11)
  yVec <- data$Y
  xStar <- data$X_star
  IFMat <- matrix(nrow = length(yVec), ncol = 2)
  for (i in 1:length(yVec))
  {
    pxs1y <- prob_x_star_1_y[i]
    pxs0y <- 1-pxs1y
    pxy <- c(pinv %*% matrix(c(pxs0y, pxs1y), ncol = 1))
    px0y <- pxy[1]
    px1y <- pxy[2]
    if(xStar[i] == 1)
    {
      px1yxs <- p11/pxs1y*px1y
      px0yxs <- (1-p00)/pxs1y*px0y
    }
    else
    {
      px1yxs <- (1-p11)/pxs0y*px1y
      px0yxs <- p00/pxs0y*px0y
    }
    u0 <- Compute_U_Logistic(beta_hat, 0, yVec[i])
    u1 <- Compute_U_Logistic(beta_hat, 1, yVec[i])
    c(cbind(u0, u1) %*% matrix(c(px0yxs, px1yxs), ncol = 1)) -> IFMat[i,]
  }
  return(IFMat)
}

Compute_IF_By_Solve_Method2_ColSum <-
  function(beta_hat,
           data,
           p00,
           p11,
           probit = F,
           ommit_intercept = F)
  {
    Compute_IF_By_Solve_Method2(beta_hat,
                                data,
                                p00,
                                p11,
                                probit = probit,
                                ommit_intercept = ommit_intercept) -> ifMat
    sum(colMeans(ifMat) ^ 2)
  }

# Oracle



# Naive


