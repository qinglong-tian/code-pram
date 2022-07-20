library(tidyverse)
library(readxl)
library(parallel)
source("application_functions.R")
Rcpp::sourceCpp("application_functions_fast.cpp")
#####################################
extract_info <- function(results, name)
{
  sapply(results, function(x)
  {
    x[[name]]
  }) -> datMat
  apply(datMat, 1, sd)
}
#####################################
p <- 0.75
B <- 500
#####################################
## Level of Eduction: 1, 2, 3, 4, 5, 6, 7, 8
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
  ) %>% mutate(age = scale(age, center = F)[, 1],
               income2005 = scale(income2005, center = F)[, 1]) -> dat

oracle_glm <- fit_oracle_estimator(dat)
naive_glm <- fit_naive_estimator(dat)

betaVal1 <- oracle_glm$coef
betaVal2 <- naive_glm$coef
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

probVec1 <- Compute_X_ast_Given_YZ(dat, linkFunc = "probit")
optim(
  par = betaVal1,
  fn = Compute_Conditional_U_Cpp,
  method = "BFGS",
  pMat = pMat,
  pMatInv = pMatInv,
  prob_X_ast_YZ = probVec1,
  eduVec = dat$edu_star,
  ageVec = dat$age,
  yVec = dat$income2005,
  genderVec = dat$gender,
) -> optimM1

probVec2 <- Compute_X_ast_Given_YZ(dat, linkFunc = "logit")
optim(
  par = betaVal1,
  fn = Compute_Conditional_U_Cpp,
  method = "BFGS",
  pMat = pMat,
  pMatInv = pMatInv,
  prob_X_ast_YZ = probVec2,
  eduVec = dat$edu_star,
  ageVec = dat$age,
  yVec = dat$income2005,
  genderVec = dat$gender,
) -> optimM2

##
B1 <- 500
n <- nrow(dat)
r_exp_list <- lapply(1:B1, function(x) {
  rexp(n)
})

mclapply(r_exp_list, function(rexpVec)
{
  optim(
    par = betaVal1,
    fn = Compute_Efficient_Weighted_Sum_App_Cpp,
    method = "BFGS",
    pMatInv = pMatInv,
    yVec = dat$income2005,
    ageVec = dat$age,
    genderVec = dat$gender,
    eduVec = dat$edu_star,
    K = K,
    rexpVec = rexpVec
  ) -> optimEff
  
  probVec1 <- Compute_X_ast_Given_YZ(dat, linkFunc = "probit")
  optim(
    par = betaVal1,
    fn = Compute_Conditional_Weigted_U_Cpp,
    method = "BFGS",
    pMat = pMat,
    pMatInv = pMatInv,
    prob_X_ast_YZ = probVec1,
    eduVec = dat$edu_star,
    ageVec = dat$age,
    yVec = dat$income2005,
    genderVec = dat$gender,
    rexpVec = rexpVec
  ) -> optimM1
  
  probVec2 <- Compute_X_ast_Given_YZ(dat, linkFunc = "logit")
  optim(
    par = betaVal1,
    fn = Compute_Conditional_Weigted_U_Cpp,
    method = "BFGS",
    pMat = pMat,
    pMatInv = pMatInv,
    prob_X_ast_YZ = probVec2,
    eduVec = dat$edu_star,
    ageVec = dat$age,
    yVec = dat$income2005,
    genderVec = dat$gender,
    rexpVec = rexpVec
  ) -> optimM2
  
  return(list(
    betaEff = optimEff$par,
    betaM1 = optimM1$par,
    betaM2 = optimM2$par
  ))
},
mc.cores = detectCores()) -> results_var

saveRDS(results_var, file = paste("dat_app/", "app_", p, ".RDS", sep = ""))
load("dat_app/", "app_", p, ".RDS", sep = "")

betaVal1
betaVal2
optimEff$par
optimM1$par
optimM2$par

sqrt(diag(vcov(oracle_glm$fit)))
sqrt(diag(vcov(naive_glm$fit)))
extract_info(results_var, "betaEff")
extract_info(results_var, "betaM1")
extract_info(results_var, "betaM2")
