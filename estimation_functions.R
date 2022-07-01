# Fit X*~Y: for X is binary outcome
fit_x_star_y_logistic <- function(data)
{
  fit <- glm(X_star ~ Y, family = "binomial", data = data)
  
  return(fit)
}

# Solving P(X|Y): no adjustment
Compute_P_Inv <- function(p00, p11)
{
  P <- matrix(ncol = 2, nrow = 2)
  P[1, 1] <- p00
  P[1, 2] <- 1 - p11
  P[2, 1] <- 1 - p00
  P[2, 2] <- p11
  
  return(solve(P))
}

# Oracle

Compute_Oracle_Beta_Logistic <- function(data)
{
  fit <- glm(X_original ~ Y, family = binomial(link = "logit"), data = data)
  fit$coefficients
}

# Naive

Compute_Naive_Beta_Logistic <- function(data)
{
  fit <- fit_x_star_y_logistic(data)
  fit$coefficients
}
