library(parallel)
library(Rcpp)
sourceCpp("estimation_functions_fast.cpp")
source("data_generating_functions.R")
source("estimation_functions.R")
#######################################
yMu <- 0.5
ySigma <- 1
betaXY <- c(-1, 1.5)
B1 <- 10
B2 <- 10
#######################################
n <- 1000
p00 <- .75
p11 <- .75
#######################################
data_mc_list <- mclapply(1:B1, function(x) {
  generate_dat_logistic(n,
                        yMu,
                        ySigma,
                        betaXY,
                        p00,
                        p11,
                        verbose = F,
                        probit = F)
},
mc.cores = detectCores())

rexp_vec_list <- mclapply(1:B2, function(x) {
  rexp(n)
})

QMat <- Compute_Q_CPP(p00, p11)

mclapply(data_mc_list, function(data)
{
  Compute_P_Inv(p00, p11) -> pinv
  
  glm(X_star~Y, family = binomial(link = "logit"), data = data) -> glmf1
  glmf1$fitted.values -> prob_x_star_1_y_1
  
  glm(X_star~0+Y, family = binomial(link = "logit"), data = data) -> glmf2
  glmf1$fitted.values -> prob_x_star_1_y_2
  
  A_Beta <-
    optim(
      betaXY,
      Compute_Eff_IF_Logisitic_Sum_CPP,
      yVec = data$Y,
      xVec = data$X_star,
      QMat = QMat
    )$par
  C1_Beta <-
    optim(
      betaXY,
      Compute_IF_Solving_Sum_CPP,
      p00 = p00,
      p11 = p11,
      pinv = pinv,
      prob_x_star_1_y = prob_x_star_1_y_1,
      yVec = data$Y,
      xStar = data$X_star
    )$par
  C2_Beta <-
    optim(
      betaXY,
      Compute_IF_Solving_Sum_CPP,
      p00 = p00,
      p11 = p11,
      pinv = pinv,
      prob_x_star_1_y = prob_x_star_1_y_2,
      yVec = data$Y,
      xStar = data$X_star
    )$par
  Oracle <- Compute_Oracle_Beta_Logistic(data)
  Naive <- Compute_Naive_Beta_Logistic(data)
  return(
    list(
      BetaEff = A_Beta,
      BetaEU = C1_Beta,
      BetaEU0 = C2_Beta,
      Oracle = Oracle,
      Naive = Naive
    )
  )
},
mc.cores = detectCores()) -> results

mclapply(data_mc_list, function(data) {
  betaHatMat <- matrix(nrow = B2, ncol = length(betaXY))
  for (i in 1:B2)
  {
    rexpVec <- rexp_vec_list[[i]]
    betaHat <-
      optim(
        betaXY,
        Compute_Eff_IF_Logisitic_Weighted_Sum_CPP,
        yVec = data$Y,
        xVec = data$X_star,
        QMat = QMat,
        rexpVec = rexpVec
      )$par
    betaHatMat[i, ] <- betaHat
  }
  return(betaHatMat)
},
mc.cores = detectCores()) -> results2

saveRDS(
  list(Estimates = results, Pert = results2),
  file = paste("dat1/results_n_", n, "_p00_", p00, "_p11_", p11, "_.RDS")
)