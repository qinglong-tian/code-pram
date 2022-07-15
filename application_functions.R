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

# Compute U Matrix Given age And gender
Compute_U_Mat_App <- function(betaVal, yVal, age, gender, K)
{
  gender2 <- ifelse(gender == 2, 1, 0)
  eduMat <- rbind(0, diag(K - 1))
  xMat <- cbind(1, age, eduMat, gender2)
  sumVec <- xMat %*% matrix(betaVal, ncol = 1)
  resid <- matrix(yVal - sumVec,
                  nrow = K,
                  ncol = K + 2,
                  byrow = F)
  
  resid * xMat
}

# Compute Efficient Score
Compute_Efficient_Score_App <-
  function(betaVal,
           pMatInv,
           yVec,
           ageVec,
           genderVec,
           eduVec,
           K)
  {
    num <- length(ageVec)
    EffMat <- matrix(nrow = num, ncol = K + 2)
    for (i in 1:num)
    {
      age <- ageVec[i]
      gender <- genderVec[i]
      yVal <- yVec[i]
      edu <- eduVec[i]
      pVec <- pMatInv[, edu]
      uMat <- Compute_U_Mat_App(betaVal, yVal, age, gender, K)
      EffMat[i,] <- c(t(uMat) %*% pVec)
    }
    
    return(EffMat)
  }

Compute_Efficient_Sum_App <-
  function(betaVal,
           pMatInv,
           yVec,
           ageVec,
           genderVec,
           eduVec,
           K)
  {
    Compute_Efficient_Score_App(betaVal, pMatInv, yVec, ageVec, genderVec, eduVec, K) -> mat
    sum(colMeans(mat) ^ 2)
  }

Compute_U_j_X_k_Z_i_Y_i <-
  function(betaVal, y, edu, age, gender, j, K)
  {
    eduMat <- rbind(0, diag(K - 1))
    gender2 <- ifelse(gender == 2, 1, 0)
    xVec <- c(1, age, eduMat[edu,], gender2)
    resid <- y - sum(betaVal * xVec)
    resid * xVec[j]
  }

Compute_d_U_j_X_k_Z_i_Y_i <-
  function(betaVal, y, edu, age, gender, j, m, K)
  {
    eduMat <- rbind(0, diag(K - 1))
    gender2 <- ifelse(gender == 2, 1, 0)
    xVec <- c(1, age, eduMat[edu,], gender2)

    - xVec[m] * xVec[j]
  }

Compute_Efficient_Sum_Deriv_App <-
  function(betaVal,
           pMatInv,
           yVec,
           ageVec,
           genderVec,
           eduVec,
           K)
  {
    firstTerm <- 0
    n <- length(ageVec)
    for (i in 1:n)
    {
      y <- yVec[i]
      age <- ageVec[i]
      gender <- genderVec[i]
      edu <- genderVec[i]
      for (k in 1:K)
      {
        for (j in 1:(K + 2))
        {
          w <- pMatInv[k, edu]
          u <-
            Compute_U_j_X_k_Z_i_Y_i(betaVal, y, k, age, gender, j, K)
          firstTerm <- firstTerm + w * u
        }
      }
    }
    firstTerm <- 2 * firstTerm

    derivVec <- numeric(K + 2)
    for (m in 1:(K + 2))
    {
      for (i in 1:n)
      {
        y <- yVec[i]
        age <- ageVec[i]
        gender <- genderVec[i]
        edu <- genderVec[i]
        for (k in 1:K)
        {
          for (j in 1:(K + 2))
          {
            w <- pMatInv[k, edu]
            du <-
              Compute_d_U_j_X_k_Z_i_Y_i(betaVal, y, k, age, gender, j, m, K)
            derivVec[m] <- derivVec[m] + w * du
          }
        }
      }
    }

    derivVec * firstTerm
  }
