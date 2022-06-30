#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix Compute_Q_CPP(double p00, double p11)
{
  double p10 = 1-p00, p01 = 1-p11;
  double det = p00*p11-p10*p01;
  
  NumericMatrix QMat(2,2);
  QMat(0,0) = p11/det;
  QMat(1,1) = p00/det;
  QMat(0,1) = -p10/det;
  QMat(1,0) = -p01/det;
  
  return QMat;
}

// [[Rcpp::export]]

NumericMatrix Compute_U_Logistic_CPP(NumericVector beta_hat, double yVal)
{
  NumericMatrix UMat(2,2);
  double expit = 1/(1+exp(-beta_hat(0)-beta_hat(1)*yVal));
  UMat(0,0) = 0-expit;
  UMat(0,1) = (0-expit)*yVal;
  
  UMat(1,0) = 1-expit;
  UMat(1,1) = (1-expit)*yVal;
  
  return UMat;
}

// [[Rcpp::export]]

NumericMatrix Compute_Eff_IF_Logistic_CPP(NumericVector beta_hat, NumericMatrix QMat, NumericVector yVec, NumericVector xVec)
{
  int num = xVec.length();
  double yVal;
  NumericMatrix UMat(2,2);
  NumericMatrix CMat(2,2);
  NumericMatrix IFMat(num, 2);
  for (int i=0; i<num; i++)
  {
    yVal = yVec(i);
    UMat = Compute_U_Logistic_CPP(beta_hat, yVal);
    IFMat(i,0) = QMat(xVec(i),0)*UMat(0,0)+QMat(xVec(i),1)*UMat(1,0);
    IFMat(i,1) = QMat(xVec(i),0)*UMat(0,1)+QMat(xVec(i),1)*UMat(1,1);
  }
  return IFMat;
}

// [[Rcpp::export]]

double Compute_Eff_IF_Logisitic_Sum_CPP(NumericVector beta_hat, NumericVector yVec, NumericVector xVec, NumericMatrix QMat)
{
  NumericMatrix IFMat = Compute_Eff_IF_Logistic_CPP(beta_hat, QMat, yVec, xVec);
  int num = IFMat.nrow();
  double sum1=0, sum2=0;

  for (int i=0; i<num; i++)
  {
    sum1 += IFMat(i,0);
    sum2 += IFMat(i,1);
  }
  sum1 /= num;
  sum2 /= num;
  double out = sum1*sum1+sum2*sum2;
  
  return out;
}

// [[Rcpp::export]]

double Compute_Eff_IF_Logisitic_Weighted_Sum_CPP(NumericVector beta_hat, NumericVector yVec, NumericVector xVec, NumericMatrix QMat, NumericVector rexpVec)
{
  NumericMatrix IFMat = Compute_Eff_IF_Logistic_CPP(beta_hat, QMat, yVec, xVec);
  int num = IFMat.nrow();
  double sum1=0, sum2=0;
  
  for (int i=0; i<num; i++)
  {
    sum1 += IFMat(i,0)*rexpVec(i);
    sum2 += IFMat(i,1)*rexpVec(i);
  }
  sum1 /= num;
  sum2 /= num;
  double out = sum1*sum1+sum2*sum2;
  
  return out;
}


