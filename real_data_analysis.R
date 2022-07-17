library(tidyverse)
library(readxl)
source("application_functions.R")
#####################################
p <- 0.9
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
  ) -> dat

# Proposed Method

betaVal1 <- fit_oracle_estimator(dat)$coef
betaVal2 <- fit_naive_estimator(dat)$coef
pMat <- make_p_matrix(p, K)
pMatInv <- solve(pMat)

optim(
  par = betaVal1,
  fn = Compute_Efficient_Sum_App,
  method = "BFGS",
  pMatInv = pMatInv,
  yVec = dat$income2005,
  ageVec = dat$age,
  genderVec = dat$gender,
  eduVec = dat$edu_star,
  K = K
)

optim(
  par = betaVal1,
  fn = Compute_Conditional_U,
  method = "BFGS",
  dat = dat,
  pMat = pMat,
  pMatInv = pMatInv,
  K = K,
  linkFunc = "probit"
)
