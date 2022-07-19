library(tidyverse)
library(readxl)
library(parallel)
source("application_functions.R")
Rcpp::sourceCpp("application_functions_fast.cpp")
#####################################
p <- .95
B <- 500
#####################################
## Level of Eduction: 1, 2, 3, 4, 5, 6, 7, 8
Sys.time() -> t0
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
  
  probVec <- Compute_X_ast_Given_YZ(dat, linkFunc = "probit")
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

Sys.time()-t0

#########
collect_results <- function(tmp, name)
{
  sapply(tmp, function(x)
  {
    x[[name]]
  }) -> outMat
  t(outMat)
}

compute_mse <- function(datMat, betaval)
{
  newDat <- datMat - matrix(betaval, nrow = nrow(datMat), ncol = ncol(datMat), byrow = T)
  newDat <- newDat*newDat
  colMeans(newDat)
}

compute_bias <- function(datMat, betaVal)
{
  colMeans(datMat)-betaVal
}
#########

p <- 0.95

tmp <- readRDS(paste("application_results_p_", p, ".RDS", sep = ""))
effMat <- collect_results(tmp, "betaEff")
effbias95 <- compute_bias(effMat, betaVal1)
effmse95 <- compute_mse(effMat, betaVal1)

mMat <- collect_results(tmp, "betaM")
mbias95 <- compute_bias(mMat, betaVal1)
mmse95 <- compute_mse(mMat, betaVal1)

nMat <- collect_results(tmp, "betaB")
nbias95 <- compute_bias(nMat, betaVal1)
nmse95 <- compute_mse(nMat, betaVal1)


p <- 0.85

tmp <- readRDS(paste("application_results_p_", p, ".RDS", sep = ""))
effMat <- collect_results(tmp, "betaEff")
effbias85 <- compute_bias(effMat, betaVal1)
effmse85 <- compute_mse(effMat, betaVal1)

mMat <- collect_results(tmp, "betaM")
mbias85 <- compute_bias(mMat, betaVal1)
mmse85 <- compute_mse(mMat, betaVal1)

nMat <- collect_results(tmp, "betaB")
nbias85 <- compute_bias(nMat, betaVal1)
nmse85 <- compute_mse(nMat, betaVal1)


p <- 0.75

tmp <- readRDS(paste("application_results_p_", p, ".RDS", sep = ""))
effMat <- collect_results(tmp, "betaEff")
effbias75 <- compute_bias(effMat, betaVal1)
effmse75 <- compute_mse(effMat, betaVal1)

mMat <- collect_results(tmp, "betaM")
mbias75 <- compute_bias(mMat, betaVal1)
mmse75 <- compute_mse(mMat, betaVal1)

nMat <- collect_results(tmp, "betaB")
nbias75 <- compute_bias(nMat, betaVal1)
nmse75 <- compute_mse(nMat, betaVal1)

##

effbias <- t(rbind(effbias75, effbias85, effbias95))
mbias <- t(rbind(mbias75, mbias85, mbias95))
nbias <- t(rbind(nbias75, nbias85, nbias95))

biasNew <- NULL
for(i in 1:4)
{
  biasNew <- rbind(biasNew, effbias[i,], mbias[i,], nbias[i,])
}
xtable(biasNew, digits = 3)

effmse <- t(rbind(effmse75, effmse85, effmse95))
mmse <- t(rbind(mmse75, mmse85, mmse95))
nmse <- t(rbind(nmse75, nmse85, nmse95))
mseNew <- NULL
for(i in 1:4)
{
  mseNew <- rbind(mseNew, effmse[i,], mmse[i,], nmse[i,])
}

latexOut <- NULL
for(i in 1:3)
{
  latexOut <- cbind(latexOut, biasNew[,i], mseNew[,i])
}
