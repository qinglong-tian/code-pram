library(parallel)
source("data_generating_functions.R")
source("estimation_functions.R")

colSds <- function(mat)
{
  apply(mat, MARGIN = 2, sd)
}

yMu <- 0
ySigma <- 1
betaXY <- c(0, 0)
n <- 400

B1 <- 200

p00 <- 0.8
p11 <- 0.8

data_mc_list <- mclapply(1:B1, function(x) {
  generate_dat_logistic(n, yMu, ySigma, betaXY, p00, p11, verbose = F, probit = F)
},
mc.cores = detectCores())

dat1 <- data_mc_list[[1]]
PX_ast <- c(1-mean(dat1$X_star), mean(dat1$X_star))
PX <- c(1-mean(dat1$X_original), mean(dat1$X_original))
table(dat1$X_original)/n

Compute_X_Given_X_ast(PX_ast, p00, p11) -> m1
t(m1)

Compute_X_Given_X_ast_Oracle(PX_ast, PX, p00, p11) -> m1_oracle
t(m1_oracle)

Compute_Q(p00, p11) -> m2
m2

mclapply(data_mc_list, function(data)
  {
  A_Beta <- optim(betaXY, Compute_IF_Logistic_ColSum, data = data, p00 = p00, p11 = p11, no_omega = T, use_q = T)$par
  A_Beta_No_Q <- optim(betaXY, Compute_IF_Logistic_ColSum, data = data, p00 = p00, p11 = p11, no_omega = T, use_q = F)$par
  B1_Beta <- optim(betaXY, Compute_IF_By_Solving_Colsum, data = data, p00 = p00, p11 = p11)$par
  B2_Beta <- optim(betaXY, Compute_IF_By_Solving_Colsum, data = data, p00 = p00, p11 = p11, probit = T)$par
  C1_Beta <- optim(betaXY, Compute_IF_By_Solve_Method2_ColSum, data = data, p00 = p00, p11 = p11, probit = F)$par
  C2_Beta <- optim(betaXY, Compute_IF_By_Solve_Method2_ColSum, data = data, p00 = p00, p11 = p11, probit = T)$par
  # D_Beta <- Estimate_Beta_EM(data, p00, p11, probit = F, tol = 1e-7)
  # E_Beta <- glm(X_original~Y, family = binomial(link = "logit"), data = data)
  return(list(
    A1 = A_Beta,
    A2 = A_Beta_No_Q,
    B1 = B1_Beta,
    B2 = B2_Beta,
    C1 = C1_Beta,
    C2 = C2_Beta
    # D_Beta,
    # E_Beta
  ))
},
mc.cores = detectCores()-2) -> results

A1 <- t(sapply(results, function(x) {x$A1}))
A2 <- t(sapply(results, function(x) {x$A2}))

B1 <- t(sapply(results, function(x) {x$B1}))
B2 <- t(sapply(results, function(x) {x$B2}))


# B1 <- t(sapply(results, function(x) {x[[2]]}))
# B2 <- t(sapply(results, function(x) {x[[3]]}))
C1 <- t(sapply(results, function(x) {x$C1}))
C2 <- t(sapply(results, function(x) {x$C2}))
# D <- t(sapply(results, function(x) {x[[6]]}))
# E <- t(sapply(results, function(x) {x[[7]]$coefficients}))

# colSds(A)
# # colSds(B1)
# # colSds(B2)
# # colSds(D)
# # colSds(E)

colMeans(A1)
colMeans(A2)



# colMeans(B1)
# colMeans(B2)
colMeans(C1)
colMeans(C2)

colSds(A1)
colSds(A2)
colSds(C1)
colSds(C2)

# colMeans(D)
# colMeans(E)

colMeans(A1)
colMeans(A2)
colMeans(C1)
colMeans(C2)
