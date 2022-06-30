#######################################
# Useful Snippets
colSds <- function(mat)
{
  apply(mat, MARGIN = 2, sd)
}

remove_outliers <- function(mat)
{
  IQRs <- apply(mat, 2, IQR)
  Bounds <- apply(mat, 2, quantile, prob = c(0.25, 0.75))
  Bounds[1, ] <- Bounds[1, ] - 2 * IQRs
  Bounds[2, ] <- Bounds[2, ] + 2 * IQRs
  
  mat_clean <- NULL
  removed <- NULL
  for (i in 1:nrow(mat))
  {
    bool_b1 <- (mat[i, 1] > Bounds[1, 1] & mat[i, 1] < Bounds[2, 1])
    bool_b2 <- (mat[i, 2] > Bounds[1, 2] & mat[i, 2] < Bounds[2, 2])
    
    if (bool_b1 & bool_b2)
    {
      mat_clean <- rbind(mat_clean, mat[i, ])
    }
    else
    {
      removed <- c(removed, i)
    }
  }
  cat(paste(100*length(removed)/nrow(mat),"% have been removed!", sep = ""))
  return(list(mat_clean = mat_clean, removed = removed))
}

Compute_Efficiency <- function(mat1_d, mat2_n, trueVal)
{
  mat1_d <- mat1_d$mat_clean
  mat2_n <- mat2_n$mat_clean
  
  num1 <- nrow(mat1_d)
  num2 <- nrow(mat2_n)
  
  true1 <- matrix(trueVal, nrow = num1, ncol = length(trueVal), byrow = T)
  true2 <- matrix(trueVal, nrow = num2, ncol = length(trueVal), byrow = T)

  colMeans((mat1_d-true1)^2) -> m1
  colMeans((mat2_n-true2)^2) -> m2
  
  return(m2/m1)
}

#######################################
A1 <- t(sapply(results, function(x) {x$BetaEff}))
A1 <- remove_outliers(A1)

C1 <- t(sapply(results, function(x) {x$BetaEU}))
C1 <- remove_outliers(C1)

C2 <- t(sapply(results, function(x) {x$BetaEU0}))
C2 <- remove_outliers(C2)

Oracle <- t(sapply(results, function(x) {x$Oracle}))
Oracle <- remove_outliers(Oracle)
Naive <- t(sapply(results, function(x) {x$Naive}))
Naive <- remove_outliers(Naive)

colMeans(A1$mat_clean)
colMeans(C1$mat_clean)
colMeans(C2$mat_clean)
colMeans(Oracle$mat_clean)
colMeans(Naive$mat_clean)

colSds(A1$mat_clean)
colSds(C1$mat_clean)
colSds(C2$mat_clean)
colSds(Naive$mat_clean)
colSds(Oracle$mat_clean)

Compute_Efficiency(C1, A1, trueVal = betaXY)
