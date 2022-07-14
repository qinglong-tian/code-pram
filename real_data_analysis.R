library(tidyverse)
library(fastDummies)
source("application_functions.R")
#####################################
p <- 0.9
#####################################
## Level of Eduction: 1, 2, 3, 4, 5, 6, 7, 8
read_excel("real_dat/dataforuse.xlsx") %>% mutate(
  income2005 = aux,
  edu_level = `level of education`-1
) %>% filter(
  edu_level != 98
) %>% select(
  age, edu_level, gender, income2005
) %>% na.omit %>% mutate(
  gender = factor(gender),
  edu_level = factor(edu_level)
) -> dat
K <- length(unique(dat$edu_level))
pMat <- make_p_matrix(p, K)
dat %>% mutate(edu_star = x_2_x_ast(edu_level, pMat)) -> dat
dat <- dummy_cols(dat, select_columns = c("edu_level", "edu_star"))

fit_oracle_estimator(dat)$coef
fit_naive_estimator(dat)$coef

# Proposed Method

