library(parallel)
library(Rcpp)
sourceCpp("estimation_functions_fast.cpp")
source("data_generating_functions.R")
source("estimation_functions.R")
source("summarize_beta_logistic.R")
####################################
index <- commandArgs(TRUE)
index <- as.numeric(index[1])
pValVec <- seq(from = 0.75, to = 0.95, by = 0.01)
ppList <- expand.grid(pValVec, pValVec)
####################################
n <- 1000
index <- 29
p00 <- ppList[index,1]
p11 <- ppList[index,2]
px <- 0.5
betaYX <- c(-1,1)
sdYX <- 1
B1 <- 100
B2 <- 500

t0 <- Sys.time()
####################################
data_mc_list <- mclapply(1:B1, function(x) {
  generate_dat_linear(n, px, betaYX, sdYX, p00, p11)
},
mc.cores = detectCores())

rexp_vec_list <- mclapply(1:B2, function(x) {
  rexp(n)
})

QMat <- Compute_Q_CPP(p00, p11)
Compute_P_Inv(p00, p11) -> pinv

mclapply(data_mc_list, function(data)
{
  glm(X_star~Y, family = binomial(link = "logit"), data = data) -> glmf1
  glmf1$fitted.values -> prob_x_star_1_y_1
  
  glm(X_star~0+Y, family = binomial(link = "logit"), data = data) -> glmf2
  glmf2$fitted.values -> prob_x_star_1_y_2
  
  A_Beta <-
    optim(
      betaYX,
      Compute_Eff_IF_Linear_Sum_CPP,
      yVec = data$Y,
      xVec = data$X_star,
      QMat = QMat
    )$par
  
  C1_Beta <-
    optim(
      betaYX,
      Compute_IF_Solving_Linear_Sum_CPP,
      p00 = p00,
      p11 = p11,
      pinv = pinv,
      prob_x_star_1_y = prob_x_star_1_y_1,
      yVec = data$Y,
      xStar = data$X_star
    )$par
  
  C2_Beta <-
    optim(
      betaYX,
      Compute_IF_Solving_Linear_Sum_CPP,
      p00 = p00,
      p11 = p11,
      pinv = pinv,
      prob_x_star_1_y = prob_x_star_1_y_2,
      yVec = data$Y,
      xStar = data$X_star
    )$par
  
  Oracle <- Compute_Oracle_Beta_Linear(data)
  names(Oracle) <- NULL
  
  Naive <- Compute_Naive_Beta_Linear(data)
  names(Naive) <- NULL
  
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
  betaHatMat <- matrix(nrow = B2, ncol = length(betaYX))
  for (i in 1:B2)
  {
    rexpVec <- rexp_vec_list[[i]]
    betaHat <-
      optim(
        betaYX,
        Compute_Eff_IF_Linear_Weighted_Sum_CPP,
        yVec = data$Y,
        xVec = data$X_star,
        QMat = QMat,
        rexpVec = rexpVec
      )$par
    betaHatMat[i,] <- betaHat
  }
  return(betaHatMat)
},
mc.cores = detectCores()) -> results2

####################################
Sys.time()-t0

dir2dat <- "dat4/"
filename <- paste("results_n_", n, "_p00_", p00, "_p11_", p11, "_.RDS", sep = "")

saveRDS(
  list(Estimates = results, Pert = results2),
  file = paste(dir2dat, filename, sep = "")
)

savedir <- "dat4/"
savename <- paste("df_", index, ".RDS", sep = "")
saveRDS(
  Read_in_Data(dir2dat, filename, betaYX),
  file = paste(savedir, savename, sep = "")
)

file.remove(paste(dir2dat, filename, sep = ""))