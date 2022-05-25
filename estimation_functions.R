find_true_EX_1 <- function(yMu, ySigma, betaXY, B=10000)
{
  yVec <- rnorm(B, yMu, ySigma)
  oddVec <- cbind(1, yVec) %*% matrix(betaXY, ncol = 1)
  probVec <- 1/(1+exp(-oddVec))
  xVec <- as.numeric(runif(B) < probVec)
  
  mean(xVec)
}
