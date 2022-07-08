#######################################
library(stringr)
# Useful Snippets
colSds <- function(mat)
{
  apply(mat, MARGIN = 2, sd)
}

remove_outliers <- function(mat)
{
  IQRs <- apply(mat, 2, IQR)
  Bounds <- apply(mat, 2, quantile, prob = c(0.25, 0.75))
  Bounds[1,] <- Bounds[1,] - 2.5 * IQRs
  Bounds[2,] <- Bounds[2,] + 2.5 * IQRs
  
  mat_clean <- NULL
  removed <- NULL
  for (i in 1:nrow(mat))
  {
    bool_b1 <- (mat[i, 1] > Bounds[1, 1] & mat[i, 1] < Bounds[2, 1])
    bool_b2 <- (mat[i, 2] > Bounds[1, 2] & mat[i, 2] < Bounds[2, 2])
    
    if (bool_b1 & bool_b2)
    {
      mat_clean <- rbind(mat_clean, mat[i,])
    }
    else
    {
      removed <- c(removed, i)
    }
  }
  cat(paste(100 * length(removed) / nrow(mat), "% have been removed!\n", sep = ""))
  return(list(mat_clean = mat_clean, removed = removed))
}

Compute_Efficiency <- function(mat1_d, mat2_n, trueVal)
{
  mat1_d <- mat1_d$mat_clean
  mat2_n <- mat2_n$mat_clean
  
  num1 <- nrow(mat1_d)
  num2 <- nrow(mat2_n)
  
  true1 <-
    matrix(trueVal,
           nrow = num1,
           ncol = length(trueVal),
           byrow = T)
  true2 <-
    matrix(trueVal,
           nrow = num2,
           ncol = length(trueVal),
           byrow = T)
  
  colMeans((mat1_d - true1) ^ 2) -> m1
  colMeans((mat2_n - true2) ^ 2) -> m2
  
  return(m2 / m1)
}
#######################################

Read_in_Data <- function(dir, name, trueVal)
{
  nVec <- str_match(name, "results_n_(.*?)_")[2] %>% as.numeric()
  p00Vec <- str_match(name, "_p00_(.*?)_p11")[2] %>% as.numeric()
  p11Vec <- str_match(name, "_p11_(.*?)_.RDS")[2] %>% as.numeric()
  
  n <- rep(nVec, 26)
  p00 <- rep(p00Vec, 26)
  p11 <- rep(p11Vec, 26)

  rst <- readRDS(paste(dir, name, sep = ""))
  rst1 <- rst$Estimates
  rst2 <- rst$Pert
  
  sapply(rst2, function(x) {
    xR <- remove_outliers(x)
    x <- xR$mat_clean
    apply(x, MARGIN = 2, sd)
  }) -> sdMat
  sdMat <- t(sdMat)
  sdMatR <- remove_outliers(sdMat)
  r1 <- sdMatR$removed
  
  sapply(rst1, function(x) {
    x$BetaEff
  }) -> effMat
  effMat <- t(effMat)
  effMatR <- remove_outliers(effMat)
  r2 <- effMatR$removed
  
  sapply(rst1, function(x) {
    x$BetaEU
  }) -> EUMat
  EUMat <- t(EUMat)
  EUMatR <- remove_outliers(EUMat)
  r3 <- EUMatR$removed
  
  sapply(rst1, function(x) {
    x$BetaEU0
  }) -> EU0Mat
  EU0Mat <- t(EU0Mat)
  EU0MatR <- remove_outliers(EU0Mat)
  r4 <- EU0MatR$removed
  
  sapply(rst1, function(x) {
    x$Oracle
  }) -> OracleMat
  OracleMat <- t(OracleMat)
  OracleMatR <- remove_outliers(OracleMat)
  r5 <- OracleMatR$removed

  sapply(rst1, function(x) {
    x$Naive
  }) -> NaiveMat
  NaiveMat <- t(NaiveMat)
  NaiveMatR <- remove_outliers(NaiveMat)
  r6 <- NaiveMatR$removed
  
  rmv <- unique(c(r2, r3, r4))

  if (length(union(rmv, r1)) == 0)
  {
    effMat2 <- effMat
  }
  else
  {
    sdMat <- sdMat[-union(rmv, r1), ]
    effMat2 <- effMat[-union(rmv, r1), ]
  }
  
  
  CP1 <- numeric(nrow(effMat2))
  CP1 -> CP2
  for(i in 1:nrow(effMat2))
  {
    lwb1 <- effMat2[i,1]-1.96*sdMat[i,1]
    upb1 <- effMat2[i,1]+1.96*sdMat[i,1]
    
    CP1[i] <- as.numeric(trueVal[1] < upb1 & trueVal[1] > lwb1)
    
    lwb2 <- effMat2[i,2]-1.96*sdMat[i,2]
    upb2 <- effMat2[i,2]+1.96*sdMat[i,2]
    
    CP2[i] <- as.numeric(trueVal[2] < upb2 & trueVal[2] > lwb2)
  }
  
  if (length(rmv) != 0)
  {
    NaiveMat <- NaiveMat[-rmv, ]
    OracleMat <- OracleMat[-rmv, ]
    EU0Mat <- EU0Mat[-rmv, ]
    EUMat <- EUMat[-rmv, ]
    effMat <- effMat[-rmv, ]
  }
  
  EffMean <- colMeans(effMat)
  EffSd <- colSds(effMat)
  
  EUMean <- colMeans(EUMat)
  EUSd <- colSds(EUMat)
  
  EU0Mean <- colMeans(EU0Mat)
  EU0Sd <- colSds(EU0Mat)
  
  OracleMean <- colMeans(OracleMat)
  OracleSd <- colSds(OracleMat)
  
  NaiveMean <- colMeans(NaiveMat)
  NaiveSd <- colSds(NaiveMat)
  
  PertSd <- colMeans(sdMat)
  
  
  
  Efficiency_Eff_EU <- Compute_Efficiency(EUMatR, effMatR, trueVal)
  
  Parameters <- rep(c("beta[0]", "beta[1]"), 13)
  Values <- c(
    EffMean,
    EUMean,
    EU0Mean,
    OracleMean,
    NaiveMean,
    EffSd,
    EUSd,
    EU0Sd,
    OracleSd,
    NaiveSd,
    PertSd,
    Efficiency_Eff_EU,
    mean(CP1),
    mean(CP2)
  )
  Type <- c(rep("Estimate", 10),
            rep("SE", 10),
            rep("SD", 2),
            rep("Efficiency", 2),
            rep("Coverage", 2))
  Method <- c(rep(rep(
    c("Efficient", "U1", "U2", "Oracle", "Naive"), each = 2
  ), 2),
  rep("Est. Efficient", 2),
  rep("EffNU1", 2),
  rep("Perturbation",2))
  
  
  out <- data.frame(Parameters, Values, Type, Method, n, p00, p11)
  rownames(out) <- NULL
  
  return(out)
}

#######################################
