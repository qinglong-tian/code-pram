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

NumericMatrix Compute_U_Linear_CPP(NumericVector beta_hat, double yVal)
{
  NumericMatrix UMat(2,2);
  // For x=0
  UMat(0,0) = yVal-beta_hat(0);
  UMat(0,1) = 0;
  
  // For x=1
  UMat(1,0) = yVal-beta_hat(0)-beta_hat(1)*1;
  UMat(1,1) = yVal-beta_hat(0)-beta_hat(1)*1;
  
  return UMat;
}

// [[Rcpp::export]]

NumericMatrix Compute_Eff_IF_Logistic_CPP(NumericVector beta_hat, NumericMatrix QMat, NumericVector yVec, NumericVector xVec)
{
  int num = xVec.length();
  double yVal;
  NumericMatrix UMat(2,2);
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

NumericMatrix Compute_Eff_IF_Linear_CPP(NumericVector beta_hat, NumericMatrix QMat, NumericVector yVec, NumericVector xVec)
{
  int num = xVec.length();
  double yVal;
  NumericMatrix UMat(2,2);
  NumericMatrix IFMat(num, 2);
  for (int i=0; i<num; i++)
  {
    yVal = yVec(i);
    UMat = Compute_U_Linear_CPP(beta_hat, yVal);
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
double Compute_Eff_IF_Linear_Sum_CPP(NumericVector beta_hat, NumericVector yVec, NumericVector xVec, NumericMatrix QMat)
{
  NumericMatrix IFMat = Compute_Eff_IF_Linear_CPP(beta_hat, QMat, yVec, xVec);
  int num = IFMat.nrow();
  double sum1 = 0, sum2 = 0;
  
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

// [[Rcpp::export]]

NumericMatrix Compute_IF_Solving_CPP(NumericVector beta_hat, double p00, double p11, NumericMatrix pinv, NumericVector prob_x_star_1_y, NumericVector yVec, NumericVector xStar)
{
  int num = yVec.length();
  double pxs1y, pxs0y;
  double px0y, px1y;
  double px1yxs, px0yxs;
  NumericMatrix UMat(2,2);
  NumericMatrix IFMat(num,2);
  for(int i=0; i<num; i++)
  {
    pxs1y = prob_x_star_1_y(i);
    pxs0y = 1-pxs1y;
    px0y = pinv(0,0)*pxs0y+pinv(0,1)*pxs1y;
    px1y = pinv(1,0)*pxs0y+pinv(1,1)*pxs1y;
    if (xStar(i) == 1)
    {
      px1yxs = p11/pxs1y*px1y;
      px0yxs = (1-p00)/pxs1y*px0y;
    }
    else
    {
      px1yxs = (1-p11)/pxs0y*px1y;
      px0yxs = p00/pxs0y*px0y;
    }
    UMat = Compute_U_Logistic_CPP(beta_hat, yVec(i));
    IFMat(i,0) = UMat(0,0)*px0yxs+UMat(1,0)*px1yxs;
    IFMat(i,1) = UMat(0,1)*px0yxs+UMat(1,1)*px1yxs;
  }
  
  return IFMat;
}

// [[Rcpp::export]]

double Compute_IF_Solving_Sum_CPP(NumericVector beta_hat, double p00, double p11, NumericMatrix pinv, NumericVector prob_x_star_1_y, NumericVector yVec, NumericVector xStar)
{
  NumericMatrix nm = Compute_IF_Solving_CPP(beta_hat, p00, p11, pinv, prob_x_star_1_y, yVec, xStar);
  int num_of_row = nm.nrow();
  double tmp1 = 0, tmp2 = 0;
  for (int i=0; i< num_of_row; i++)
  {
    tmp1 += nm(i,0);
    tmp2 += nm(i,1);
  }
  tmp1 /= num_of_row;
  tmp2 /= num_of_row;
  return tmp1*tmp1+tmp2*tmp2;
}
