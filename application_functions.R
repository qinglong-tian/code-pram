# Prepare for the transformation functions
make_p_matrix <- function(p, n)
{
  pOut <- matrix(ncol = n, nrow = n)
  for (i in 1:n)
  {
    for (j in 1:n)
    {
      if (i == j)
      {
        pOut[i,j] <- p
      }
      else
      {
        pOut[i,j] <- (1-p)/(n-1)
      }
    }
  }
  return(pOut)
}

x_2_x_ast <- function(x, pMat)
{
  n <- nrow(pMat)
  pdfMat <- pMat[,c(x)]
  apply(pdfMat, 2, function(pdfVec)
  {
    sample(1:n, size = 1, replace = T, pdfVec)
  }) -> x_ast
  x_ast <- factor(x_ast)
  
  return(x_ast)
}

# Estimation
fit_naive_estimator <- function(dat)
{
  lm(income2005 ~ age + edu_level + gender, data = dat) -> lmfit
  return(list(fit = lmfit,
              coef = lmfit$coefficients))
}

fit_oracle_estimator <- function(dat)
{
  lm(income2005 ~ age + edu_star + gender, data = dat) -> lmfit
  return(list(fit = lmfit,
              coef = lmfit$coefficients))
}

