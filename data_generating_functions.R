#####################################
# The PRAM is performed on X
# Consider two settings:
# (1) X is outcome, Y is covariate
# (2) Y is covariate, Y is outcome
#####################################

x2xstar <- function(xVec, p)
{
  n <- length(xVec)
  xStarVec <- ifelse(runif(n)<p, xVec, 1-xVec)
  return(xStarVec)
}

generate_dat_1 <- function(n, yMu, ySigma, betaXY, p, verbose = F)
{
  yVec <- rnorm(n, yMu, ySigma)
  oddVec <- cbind(1, yVec) %*% matrix(betaXY, ncol = 1)
  probVec <- 1/(1+exp(-oddVec))
  xVec <- as.numeric(runif(n) < probVec)
  xStarVec <- x2xstar(xVec, p)
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
