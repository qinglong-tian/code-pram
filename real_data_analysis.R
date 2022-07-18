library(tidyverse)
library(readxl)
library(parallel)
source("application_functions.R")
#####################################
p <- 0.9
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

data_list_mc <- lapply(1:B, function(x)
{
  dat %>%
    mutate(edu_star = x_2_x_ast(edu_new_level, pMat)) %>%
    mutate(
      gender = factor(gender),
      edu_new_level = factor(edu_new_level),
      edu_star = factor(edu_star)
    ) -> dat1
  return(dat1)
})

# Proposed Method
betaVal1 <- fit_oracle_estimator(dat)$coef

pMat <- make_p_matrix(p, K)
pMatInv <- solve(pMat)

mclapply(data_list_mc, function(dat)
{
  betaVal2 <- fit_naive_estimator(dat)$coef
  
  yVec = dat$income2005
  ageVec = dat$age
  genderVec = dat$gender
  eduVec = dat$edu_star
  
  optim(
    par = betaVal1,
    fn = Compute_Efficient_Sum_App,
    pMatInv = pMat,
    yVec = yVec,
    ageVec = ageVec,
    genderVec = genderVec,
    eduVec = eduVec,
    K = 2
  )$par -> betaEff
  
  optim(
    par = betaVal1,
    fn = Compute_Conditional_U,
    dat = dat,
    pMat = pMat,
    pMatInv = pMatInv,
    K = 2,
    linkFunc = "logit"
  )$par -> betaM1
  
  return(list(
    betaEff = betaEff,
    betaM = betaM1,
    betaNaive = betaVal2
  ))
},
mc.cores = detectCores()) -> results

saveRDS(results, file = paste("application_results_p_", p, ".RDS", sep = ""))

