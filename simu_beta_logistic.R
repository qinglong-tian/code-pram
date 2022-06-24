library(parallel)
source("data_generating_functions.R")
source("estimation_functions.R")

yMu <- 0
ySigma <- 1
betaXY <- c(-1, 1)
n <- 500

B1 <- 200

p00 <- 0.9
p11 <- 0.9

data_mc_list <- mclapply(1:B1, function(x) {
  generate_dat_logistic(n, yMu, ySigma, betaXY, p00, p11, verbose = F)
},
mc.cores = detectCores())

mclapply(data_mc_list, function(data)
  {
  beta1 <- optim(betaXY, Compute_IF_By_Solving_Colsum, data = data, p00 = p00, p11 = p11)$par
  beta2 <- optim(betaXY, Compute_IF_By_Solving_Colsum, data = data, p00 = p00, p11 = p11, probit = T)$par
  
  beta3 <- optim(betaXY, Compute_IF_Logistic_ColSum, data = data, p00 = p00, p11 = p11)$par
  
  beta4 <- optim(betaXY, Compute_IF_By_Solve_Method2_ColSum, data = data, p00 = p00, p11 = p11, probit = F)$par
  beta41 <- optim(betaXY, Compute_IF_By_Solve_Method2_ColSum, data = data, p00 = p00, p11 = p11, probit = T)$par
  
  beta5 <- Estimate_Beta_EM(data, p00, p11, probit = F, tol = 1e-7)
  beta6 <- Estimate_Beta_EM(data, p00, p11, probit = T, tol = 1e-7)
  return(list(
    Beta1 = beta1,
    Beta2 = beta2,
    Beta3 = beta3,
    Beta4 = beta4,
    Beta41 = beta41,
    Beta5 = beta5,
    Beta6 = beta6
  ))
},
mc.cores = detectCores()-2) -> results
# 
# optim(betaXY, Compute_IF_By_Solving_Colsum, data = data, p00 = p00, p11 = p11)
# optim(betaXY, Compute_IF_By_Solving_Colsum, data = data, p00 = p00, p11 = p11, probit = T)
# 
# optim(betaXY, Compute_IF_Logistic_ColSum, data = data, p00 = p00, p11 = p11)

# Method 1: Logit

apply(sapply(results, function(dat) {
  dat$Beta1
}),
MARGIN = 1,
sd)

apply(sapply(results, function(dat) {
  dat$Beta1
}),
MARGIN = 1,
mean
)

# Method 1: Probit

apply(sapply(results, function(dat) {
  dat$Beta2
}),
MARGIN = 1,
sd)

apply(sapply(results, function(dat) {
  dat$Beta2
}),
MARGIN = 1,
mean
)

# Proposed Efficient Method

apply(sapply(results, function(dat) {
  dat$Beta3
}),
MARGIN = 1,
sd)

apply(sapply(results, function(dat) {
  dat$Beta3
}),
MARGIN = 1,
mean)

# Method 2-approach1-Logit

apply(sapply(results, function(dat) {
  dat$Beta4
}),
MARGIN = 1,
sd)

apply(sapply(results, function(dat) {
  dat$Beta4
}),
MARGIN = 1,
mean)

# Method 2-approach1-probit

apply(sapply(results, function(dat) {
  dat$Beta41
}),
MARGIN = 1,
sd)

apply(sapply(results, function(dat) {
  dat$Beta41
}),
MARGIN = 1,
mean)

# Method 2-EM- Logit

apply(sapply(results, function(dat) {
  dat$Beta5
}),
MARGIN = 1,
sd)

apply(sapply(results, function(dat) {
  dat$Beta5
}),
MARGIN = 1,
mean)


# Method 2-EM- Probit

apply(sapply(results, function(dat) {
  dat$Beta6
}),
MARGIN = 1,
sd)

apply(sapply(results, function(dat) {
  dat$Beta6
}),
MARGIN = 1,
mean)

