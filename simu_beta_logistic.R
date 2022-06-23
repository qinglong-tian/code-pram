library(parallel)
source("data_generating_functions.R")
source("estimation_functions.R")

yMu <- 0
ySigma <- 1
betaXY <- c(-1, 1)
n <- 200

B1 <- 300

p00 <- 0.9
p11 <- 0.85

data_mc_list <- mclapply(1:B1, function(x) {
  generate_dat_logistic(n, yMu, ySigma, betaXY, p00, p11, verbose = F)
},
mc.cores = detectCores())

mclapply(data_mc_list, function(data)
{
  optim(
    betaXY,
    Compute_Efficient_IF_Logistic_Sum,
    dat = data,
    p00 = p00,
    p11 = p11,
    tilde_Func = Tilde_E_X_Y_Naive
  )$par -> betaNaive
  optim(
    betaXY,
    Compute_Efficient_IF_Logistic_Sum,
    dat = data,
    p00 = p00,
    p11 = p11,
    tilde_Func = Tilde_E_X_Y_Solve
  )$par -> betaSolve
  Compute_Beta_Direct_Solve(data = data, p00, p11) -> betaDirect
  glm(X_original ~ Y, data = data, family = binomial)$coefficients -> betaOracle
  
  return(
    list(
      BetaNaive = betaNaive,
      BetaSolve = betaSolve,
      BetaDirect = betaDirect,
      BetaOracle = betaOracle
    )
  )
},
mc.cores = detectCores()) -> results

#########################

library(tidyverse)
sapply(results, function (x) {
  x$BetaNaive
}) %>% rowMeans() - betaXY
sapply(results, function (x) {
  x$BetaNaive
}) %>% apply(MARGIN = 1, FUN = sd)

sapply(results, function (x) {
  x$BetaSolve
}) %>% rowMeans() - betaXY
sapply(results, function (x) {
  x$BetaSolve
}) %>% apply(MARGIN = 1, FUN = sd)

sapply(results, function (x) {
  x$BetaDirect
}) %>% rowMeans() - betaXY
sapply(results, function (x) {
  x$BetaDirect
}) %>% apply(MARGIN = 1, FUN = sd)

sapply(results, function (x) {
  x$BetaOracle
}) %>% rowMeans() - betaXY
sapply(results, function (x) {
  x$BetaOracle
}) %>% apply(MARGIN = 1, FUN = sd)
