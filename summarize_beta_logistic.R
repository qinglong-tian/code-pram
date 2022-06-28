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
  
  return(list(mat_clean = mat_clean, removed = removed))
}

Compute_Efficiency <- function(mat1_d, mat2_n, trueVal)
{
  num1 <- nrow(mat1_d)
  num2 <- nrow(mat2_n)
  
  true1 <- matrix(nrow = num1, ncol = length(trueVal))
  true2 <- matrix(nrow = num2, ncol = length(trueVal))

  colMeans(mat1_d-true1)^2 -> m1
  colMeans(mat2_n-true2)^2 -> m2
  
  return(m2/m1)
}

#######################################
A1 <- t(sapply(results, function(x) {x$A1}))
A1 <- remove_outliers(A1)$mat_clean

C1 <- t(sapply(results, function(x) {x$C1}))
C1 <- remove_outliers(C1)$mat_clean

C2 <- t(sapply(results, function(x) {x$C2}))
C2 <- remove_outliers(C2)$mat_clean

colMeans(A1)
colMeans(C1)
colMeans(C2)

colSds(A1)
colSds(C1)
colSds(C2)
