library(parallel)
source("data_generating_functions.R")
source("estimation_functions.R")
#######################################
yMu <- 0
ySigma <- 0.8
betaXY <- c(1, -1)
n <- 150

B1 <- 1000

p00 <- 0.75
p11 <- 0.75
#######################################
data_mc_list <- mclapply(1:B1, function(x) {
  generate_dat_logistic(n, yMu, ySigma, betaXY, p00, p11, verbose = F, probit = F)
},
mc.cores = detectCores())

dat1 <- data_mc_list[[1]]
PX_ast <- c(1-mean(dat1$X_star), mean(dat1$X_star))
PX <- c(1-mean(dat1$X_original), mean(dat1$X_original))
table(dat1$X_original)/n

mclapply(data_mc_list, function(data)
  {
  A_Beta <- optim(betaXY, Compute_IF_Logistic_ColSum, data = data, p00 = p00, p11 = p11, no_omega = T, use_q = T)$par
  C1_Beta <- optim(betaXY, Compute_IF_By_Solve_Method2_ColSum, data = data, p00 = p00, p11 = p11, probit = F)$par
  C2_Beta <- optim(betaXY, Compute_IF_By_Solve_Method2_ColSum, data = data, p00 = p00, p11 = p11, probit = F, ommit_intercept = T)$par
  return(list(
    A1 = A_Beta,
    C1 = C1_Beta,
    C2 = C2_Beta
  ))
},
mc.cores = detectCores()-2) -> results
