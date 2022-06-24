library(parallel)
source("data_generating_functions.R")
source("estimation_functions.R")

colSds <- function(mat)
{
  apply(mat, MARGIN = 2, sd)
}

yMu <- 0
ySigma <- 1
betaXY <- c(-1, 1)
n <- 1000

B1 <- 200

p00 <- 0.9
p11 <- 0.9

data_mc_list <- mclapply(1:B1, function(x) {
  generate_dat_logistic(n, yMu, ySigma, betaXY, p00, p11, verbose = F)
},
mc.cores = detectCores())

mclapply(data_mc_list, function(data)
  {
  A_Beta <- optim(betaXY, Compute_IF_Logistic_ColSum, data = data, p00 = p00, p11 = p11)$par
  B1_Beta <- optim(betaXY, Compute_IF_By_Solving_Colsum, data = data, p00 = p00, p11 = p11)$par
  B2_Beta <- optim(betaXY, Compute_IF_By_Solving_Colsum, data = data, p00 = p00, p11 = p11, probit = T)$par
  C1_Beta <- optim(betaXY, Compute_IF_By_Solve_Method2_ColSum, data = data, p00 = p00, p11 = p11, probit = F)$par
  C2_Beta <- optim(betaXY, Compute_IF_By_Solve_Method2_ColSum, data = data, p00 = p00, p11 = p11, probit = T)$par
  D_Beta <- Estimate_Beta_EM(data, p00, p11, probit = F, tol = 1e-7)
  E_Beta <- glm(X_original~Y, family = binomial(link = "logit"), data = data)
  return(list(
    A_Beta,
    B1_Beta,
    B2_Beta,
    C1_Beta,
    C2_Beta,
    D_Beta,
    E_Beta
  ))
},
mc.cores = detectCores()-2) -> results

A <- t(sapply(results, function(x) {x[[1]]}))
B1 <- t(sapply(results, function(x) {x[[2]]}))
B2 <- t(sapply(results, function(x) {x[[3]]}))
C1 <- t(sapply(results, function(x) {x[[4]]}))
C2 <- t(sapply(results, function(x) {x[[5]]}))
D <- t(sapply(results, function(x) {x[[6]]}))
E <- t(sapply(results, function(x) {x[[7]]$coefficients}))

colSds(A)
colSds(B1)
colSds(B2)
colSds(C1)
colSds(C2)
colSds(D)
colSds(E)

colMeans(A)
colMeans(B1)
colMeans(B2)
colMeans(C1)
colMeans(C2)
colMeans(D)
colMeans(E)
