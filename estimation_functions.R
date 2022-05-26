# Compute E(X)
find_true_EX_1 <- function(yMu, ySigma, betaXY, B=10000)
{
  yVec <- rnorm(B, yMu, ySigma)
  oddVec <- cbind(1, yVec) %*% matrix(betaXY, ncol = 1)
  probVec <- 1/(1+exp(-oddVec))
  xVec <- as.numeric(runif(B) < probVec)
  
  mean(xVec)
}

# Fit X*~Y
fit_x_star_y <- function(data)
{
  fit <- glm(X_star~Y, family = "binomial", data = data)
  
  return(fit)
}

# Compute \tilde{E}(X|Y)
# If fit_x_y is computed using fit_x_star_y(), then it is the flexible method
# If fit_x_y is computed with the X~Y, then this is the oracle
compute_tilde_EXY <- function(yVec, fit_x_y)
{
  coef_x_y <- fit_x_y$coefficients
  oddVec <- cbind(1, yVec) %*% matrix(coef_x_y, ncol = 1)
  probVec <- 1/(1+exp(-oddVec))
  
  return(probVec)
}

# Compute \Pr(X|Y) by solving linear systems
make_p_matrix <- function(p)
{
  pMat <- matrix(ncol = 2, nrow = 2)
  pMat[1,1] <- p
  pMat[2,2] <- p
  pMat[1,2] <- 1-p
  pMat[2,1] <- 1-p
  
  return(pMat)
}

# Compute [\Pr(X^\ast=0|Y), \Pr(X^\ast=1|Y)]^T first
# Then output [\Pr(X=0|Y), \Pr(X=1|Y)]^T
compute_pr_x_y <- function(y, fit_x_star_y, pMatInv)
{
  coef_x_star_y <- fit_x_star_y$coefficients
  odd <- cbind(1, y) %*% matrix(coef_x_star_y, ncol = 1)
  prob <- 1/(1+exp(-odd))
  
  v <- matrix(c(1-prob, prob), ncol = 1)
  u <- pMatInv%*%v
  return(u[2,1])
}

# Compute \tilde{E} by solving the linear systems
compute_tilde_EXY_solve <- function(yVec, fit_x_star_y, p)
{
  n <- length(yVec)
  pMatInv <- solve(make_p_matrix(p))
  out <- numeric(n)
  for(i in 1:n)
  {
    out[i] <- compute_pr_x_y(yVec[i], fit_x_star_y, pMatInv)
  }
  return(out)
}
