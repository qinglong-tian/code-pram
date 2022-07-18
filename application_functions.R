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
      EffMat[i, ] <- c(t(uMat) %*% pVec)
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

Compute_X_Given_YZ <- function(pMatInv, prob_X_ast_YZ)
{
  t(pMatInv %*% rbind(1 - prob_X_ast_YZ, prob_X_ast_YZ))
}

Compute_X_Given_X_star_YZ <-
  function(prob_X_ast_YZ, prob_x_yz, pMat, dat)
  {
    num_of_obs <- nrow(prob_x_yz)
    outMat <- matrix(nrow = num_of_obs, ncol = 2)
    for (i in 1:num_of_obs)
    {
      eduStar <- dat$edu_star[i]
      p_xast_y <-
        ifelse(eduStar == 1, 1 - prob_X_ast_YZ[i], prob_X_ast_YZ[i])
      
      p_xast_x <- pMat[eduStar, 1]
      p_x_y <- prob_x_yz[i, 1]
      rawProb <- p_xast_x / p_xast_y * p_x_y
      outMat[i, 1] <-
        ifelse(rawProb < 0,
               0,
               ifelse(rawProb > 1, 1, rawProb))
      outMat[i, 2] <- 1 - outMat[i, 1]
    }
    
    return(outMat)
  }

Compute_Conditional_U <-
  function(betaVal, dat, pMat, pMatInv, K, linkFunc = "logit")
  {
    prob_X_ast_YZ <- Compute_X_ast_Given_YZ(dat, linkFunc)
    prob_x_yz <- Compute_X_Given_YZ(pMatInv, prob_X_ast_YZ)
    condProb <-
      Compute_X_Given_X_star_YZ(prob_X_ast_YZ, prob_x_yz, pMat, dat)
    
    num_of_obs <- nrow(dat)
    yVec <- dat$income2005
    ageVec <- dat$age
    genderVec <- dat$gender
    
    outU <- matrix(nrow = num_of_obs, ncol = K + 2)
    for (i in 1:num_of_obs)
    {
      age <- ageVec[i]
      gender <- genderVec[i]
      yVal <- yVec[i]
      uMat <- Compute_U_Mat_App(betaVal, yVal, age, gender, K)
      outU[i, ] <- c(matrix(condProb[i, ], ncol = 2) %*% uMat)
    }
    
    sum(colMeans(outU) ^ 2)
  }
