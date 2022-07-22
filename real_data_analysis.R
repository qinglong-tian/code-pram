library(tidyverse)
library(readxl)
library(parallel)
source("application_functions.R")
Rcpp::sourceCpp("application_functions_fast.cpp")
#####################################
p <- 0.95
#####################################
## Level of Education: 1, 2, 3, 4, 5, 6, 7, 8
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
remove_bad_obs(dat, probVec1, pMatInv) -> dat_selected
dat1 <- dat_selected$dat1
index1 <- dat_selected$index

optim(
  par = betaVal1,
  fn = Compute_Conditional_U_Cpp,
  method = "BFGS",
  pMat = pMat,
  pMatInv = pMatInv,
  prob_X_ast_YZ = probVec1[index1],
  eduVec = dat1$edu_star,
  ageVec = dat1$age,
  yVec = dat1$income2005,
  genderVec = dat1$gender,
) -> optimM1

probVec2 <- Compute_X_ast_Given_YZ(dat, linkFunc = "logit")
remove_bad_obs(dat, probVec2, pMatInv) -> dat_selected
dat2 <- dat_selected$dat1
index2 <- dat_selected$index

optim(
  par = betaVal1,
  fn = Compute_Conditional_U_Cpp,
  method = "BFGS",
  pMat = pMat,
  pMatInv = pMatInv,
  prob_X_ast_YZ = probVec2[index2],
  eduVec = dat2$edu_star,
  ageVec = dat2$age,
  yVec = dat2$income2005,
  genderVec = dat2$gender,
) -> optimM2

oracle <- betaVal1
naive <- betaVal2
proposed <- optimEff$par
model1 <- optimM1$par
model2 <- optimM2$par

betaMat <- rbind(oracle, proposed, model1, model2, naive)
betaVec <- c(betaMat)

##
B1 <- 300
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
  remove_bad_obs(dat, probVec1, pMatInv) -> dat_selected
  dat1 <- dat_selected$dat1
  index1 <- dat_selected$index
  
  optim(
    par = betaVal1,
    fn = Compute_Conditional_Weigted_U_Cpp,
    method = "BFGS",
    pMat = pMat,
    pMatInv = pMatInv,
    prob_X_ast_YZ = probVec1[index1],
    eduVec = dat1$edu_star,
    ageVec = dat1$age,
    yVec = dat1$income2005,
    genderVec = dat1$gender,
    rexpVec = rexpVec
  ) -> optimM1
  
  probVec2 <- Compute_X_ast_Given_YZ(dat, linkFunc = "logit")
  remove_bad_obs(dat, probVec2, pMatInv) -> dat_selected
  dat2 <- dat_selected$dat1
  index2 <- dat_selected$index
  
  optim(
    par = betaVal1,
    fn = Compute_Conditional_Weigted_U_Cpp,
    method = "BFGS",
    pMat = pMat,
    pMatInv = pMatInv,
    prob_X_ast_YZ = probVec2[index2],
    eduVec = dat2$edu_star,
    ageVec = dat2$age,
    yVec = dat2$income2005,
    genderVec = dat2$gender,
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
results_var <-
  readRDS(paste("dat_app/", "app_", p, ".RDS", sep = ""))

oracle_sd <- sqrt(diag(vcov(oracle_glm$fit)))
naive_sd <- sqrt(diag(vcov(naive_glm$fit)))
proposed_sd <- extract_info(results_var, "betaEff")
model1_sd <- extract_info(results_var, "betaM1")
model2_sd <- extract_info(results_var, "betaM2")
sdMat <-
  rbind(oracle_sd, proposed_sd, model1_sd, model2_sd, naive_sd)
sdVec <- c(sdMat)

betaVec3 <- sprintf(betaVec, fmt = '%#.3f')
sdVec3 <- sprintf(sdVec, fmt = '%#.3f')
beta_sd <- paste(betaVec3, " (", sdVec3, ")", sep = "")
saveRDS(beta_sd, file = paste("p", p, ".RDS", sep = ""))

xtable::xtable(cbind(p0.75, p0.85, p0.95))
