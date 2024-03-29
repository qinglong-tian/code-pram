#####################################
# The PRAM is performed on X
# Consider two settings:
# (1) X is outcome, Y is covariate
# (2) X is covariate, Y is outcome
#####################################

x2xstar <- function(xVec, p00, p11)
# For binary x & x^\ast
{
  n <- length(xVec)
  p <- p11*xVec+p00*(1-xVec)
  xStarVec <- ifelse(runif(n)<p, xVec, 1-xVec)
  return(xStarVec)
}

generate_dat_logistic <- function(n, yMu, ySigma, betaXY, p00, p11, verbose = F, probit = F)
# X is binary outcome, Y is continuous covariate
{
  yVec <- rnorm(n, yMu, ySigma)
  oddVec <- cbind(1, yVec) %*% matrix(betaXY, ncol = 1)
  if (!probit)
  {
    probVec <- 1/(1+exp(-oddVec))
  }
  else
  {
    probVec <- pnorm(oddVec)
  }
  xVec <- as.numeric(runif(n) < probVec)
  xStarVec <- x2xstar(xVec, p00, p11)
  if(verbose)
  {
    cat(paste("The proportion of X=1 is ", 100*mean(xVec), "%.", sep = ""))
    cat("\n")
    cat(paste("The proportion of X*=1 is ", 100*mean(xStarVec), "%.", sep = ""))
    cat("\n")
  }
  
  outDF <- cbind(xVec, xStarVec, yVec)
  colnames(outDF) <- c("X_original", "X_star", "Y")
  outDF <- as.data.frame(outDF)
  return(outDF)
}

generate_dat_linear <- function(n, px, betaYX, sdYX, p00, p11, verbose = F)
# X is binary as covariate, Y is response using linear model
{
  xVec <- as.numeric(runif(n)<px)
  yVec <- c(cbind(1, xVec) %*% matrix(betaYX, ncol = 1))+rnorm(n, 0, sdYX)
  xStarVec <- x2xstar(xVec, p00, p11)
  
  if(verbose)
  {
    cat(paste("The proportion of X=1 is ", 100*mean(xVec), "%.", sep = ""))
    cat("\n")
    cat(paste("The proportion of X*=1 is ", 100*mean(xStarVec), "%.", sep = ""))
    cat("\n")
  }
  
  outDF <- cbind(xVec, xStarVec, yVec)
  colnames(outDF) <- c("X_original", "X_star", "Y")
  outDF <- as.data.frame(outDF)
  return(outDF)
}
