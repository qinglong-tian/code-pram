library(tidyverse)
library(readxl)
library(parallel)
source("application_functions.R")
Rcpp::sourceCpp("application_functions_fast.cpp")
#####################################
p <- 0.8
B <- 500
#####################################
## Level of Eduction: 1, 2, 3, 4, 5, 6, 7, 8
results <- mclapply(1:B, function(x)
{
  read_excel("real_dat/dataforuse.xlsx") %>%
    mutate(income2005 = aux, edu_level = `level of education` - 1) %>%
    filter(edu_level != 98) %>% select(age, edu_level, gender, income2005) %>%
    na.omit %>%
    mutate(edu_new_level = ifelse(edu_level <= 4, 1, 2)) -> dat
  K <- length(unique(dat$edu_new_level))
  pMat <- make_p_matrix(p, K)
  dat %>%
    mutate(edu_star = x_2_x_ast(edu_new_level, pMat)) %>%
    mutate(
      gender = factor(gender),
      edu_new_level = factor(edu_new_level),
      edu_star = factor(edu_star)
    ) -> dat
  
  # Proposed Method
  betaVal1 <- fit_oracle_estimator(dat)$coef
  betaVal2 <- fit_naive_estimator(dat)$coef
  pMatInv <- solve(pMat)
  
  optim(
    par = betaVal1,
    fn = Compute_Efficient_Sum_App_Cpp,
    method = "BFGS",
    pMatInv = pMatInv,
    yVec = dat$income2005,
    ageVec = dat$age,
    genderVec = dat$gender,
    eduVec = dat$edu_star,
    K = K
  ) -> optimEff
  
  probVec <- Compute_X_ast_Given_YZ(dat, linkFunc = "logit")
  optim(
    par = betaVal1,
    fn = Compute_Conditional_U_Cpp,
    method = "BFGS",
    pMat = pMat,
    pMatInv = pMatInv,
    prob_X_ast_YZ = probVec,
    eduVec = dat$edu_star,
    ageVec = dat$age,
    yVec = dat$income2005,
    genderVec = dat$gender,
  ) -> optimM
  
  return(list(
    betaEff = optimEff$par,
    betaM = optimM$par,
    betaB = betaVal2,
    betaO = betaVal1
  ))
},
mc.cores = detectCores())
saveRDS(results, file = paste("application_results_p_", p, ".RDS", sep = ""))

#########
collect_results <- function(tmp, name)
{
  sapply(tmp, function(x)
  {
    x[[name]]
  }) -> outMat
  t(outMat)
}
#########

tmp <- readRDS("application_results_p_0.8.RDS")
effMat <- collect_results(tmp, "betaEff")
colMeans(effMat)
mMat <- collect_results(tmp, "betaM")
colMeans(mMat)
nMat <- collect_results(tmp, "betaNaive")
colMeans(nMat)

