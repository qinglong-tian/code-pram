# Prepare for the transformation functions
make_p_matrix <- function(p, K)
{
  pOut <- matrix(ncol = K, nrow = K)
  for (i in 1:K)
  {
    for (j in 1:K)
    {
      if (i == j)
      {
        pOut[i, j] <- p
      }
      else
      {
        pOut[i, j] <- (1 - p) / (K - 1)
      }
    }
  }
  return(pOut)
}

x_2_x_ast <- function(x, pMat)
{
  n <- nrow(pMat)
  pdfMat <- pMat[, c(x)]
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
  lm(income2005 ~ age + edu_star + gender, data = dat) -> lmfit
  return(list(fit = lmfit,
              coef = lmfit$coefficients))
}

fit_oracle_estimator <- function(dat)
{
  lm(income2005 ~ age + edu_new_level + gender, data = dat) -> lmfit
  return(list(fit = lmfit,
              coef = lmfit$coefficients))
}

# The alternative method

Compute_X_ast_Given_YZ <- function(dat, linkFunc = "logit")
{
  glmfit <-
    glm(
      edu_star ~ income2005 + gender + age,
      family = binomial(link = linkFunc),
      data = dat
    )
  glmfit$fitted.values
}

# Compute SD
extract_info <- function(results, name)
{
  sapply(results, function(x)
  {
    x[[name]]
  }) -> datMat
  apply(datMat, 1, sd)
}
