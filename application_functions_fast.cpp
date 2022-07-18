#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Compute_U_Mat_App_Cpp(NumericVector betaVal, double yVal, double age, double gender, int K)
{
  int i, j;
  double gender2 = gender-1;
  NumericMatrix eduMat(K, K-1);
  for(i=0; i<K-1; i++)
  {
    for(j=0; j<K-1; j++)
    {
      if (i == j)
      {
        eduMat(i+1,j) = 1;
      }
      else
      {
        eduMat(i+1,j) = 0;
      }
    }
  }
  
  NumericMatrix xMat(K, K+2);
  for(i=0; i<K; i++)
  {
    xMat(i,0) = 1;
    xMat(i,1) = age;
    xMat(i,K+1) = gender2;
  }
  for(i=0; i<K; i++)
  {
    for(j=0; j<K-1; j++)
    {
      xMat(i,j+2) = eduMat(i,j);
    }
  }
  NumericVector sumVec(K), resid(K);
  for(i=0; i<K; i++)
  {
    for(j=0; j<K+2; j++)
    {
      sumVec(i) += betaVal(j)*xMat(i,j);
    }
  }
  
  for(i=0; i<K; i++)
  {
    resid(i) = yVal-sumVec(i);
  }
  
  NumericMatrix uMat(K, K+2);
  for(i=0; i<K; i++)
  {
    for(j=0; j<K+2; j++)
    {
      uMat(i,j) = resid(i)*xMat(i,j);
    }
  }
  
  return uMat;
}

// [[Rcpp::export]]
NumericMatrix Compute_Efficient_Score_App_Cpp(NumericVector betaVal, NumericMatrix pMatInv, NumericVector yVec, NumericVector ageVec, NumericVector genderVec, NumericVector eduVec, int K)
{
  int i,j,k;
  int num_of_obs = yVec.length();
  NumericMatrix EffMat(num_of_obs, K+2), uMat(K, K+2);
  double age, gender, yVal, edu;
  for (i=0; i<num_of_obs; i++)
  {
    age = ageVec(i);
    gender = genderVec(i);
    yVal = yVec(i);
    edu = eduVec(i);
    
    uMat = Compute_U_Mat_App_Cpp(betaVal, yVal, age, gender, K);
    for (j=0; j<K+2; j++)
    {
      for (k=0; k<K; k++)
      {
        EffMat(i,j) += uMat(k,j)*pMatInv(k,edu-1);
      }
    }
  }
  
  return EffMat;
}

// [[Rcpp::export]]
double Compute_Efficient_Sum_App_Cpp(NumericVector betaVal, NumericMatrix pMatInv, NumericVector yVec, NumericVector ageVec, NumericVector genderVec, NumericVector eduVec, int K)
{
  NumericMatrix EffMat = Compute_Efficient_Score_App_Cpp(betaVal, pMatInv, yVec, ageVec, genderVec, eduVec, K);
  int num_of_rows = EffMat.nrow();
  int num_of_cols = EffMat.ncol();
  
  double _sum=0, tmp;
  for(int i=0; i<num_of_cols; i++)
  {
    tmp =0;
    for(int j=0; j<num_of_rows; j++)
    {
      tmp += EffMat(j,i);
    }
    tmp /= num_of_rows;
    _sum += tmp*tmp;
  }
  return _sum;
}

// [[Rcpp::export]]
NumericMatrix Compute_X_Given_YZ_Cpp(NumericMatrix pMatInv, NumericVector prob_X_ast_YZ)
{
  int num_of_obs = prob_X_ast_YZ.length();
  NumericMatrix outMat(num_of_obs,2);
  double tmp;
  
  for (int i=0; i<num_of_obs; i++)
  {
    tmp= prob_X_ast_YZ(i);
    outMat(i,0) = pMatInv(0,0)*(1-tmp)+pMatInv(0,1)*tmp;
    outMat(i,1) = pMatInv(1,0)*(1-tmp)+pMatInv(1,1)*tmp;
  }
  
  return outMat;
}

// [[Rcpp::export]]
NumericMatrix Compute_X_Given_X_star_YZ_Cpp(NumericVector prob_X_ast_YZ, NumericMatrix prob_x_yz, NumericMatrix pMat, NumericVector eduVec)
{
  int num_of_obs = prob_X_ast_YZ.length();
  NumericMatrix outMat(num_of_obs, 2);
  double eduStar, p_xast_y, p_xast_x, p_x_y, rawProb;
  for(int i=0; i<num_of_obs; i++)
  {
    eduStar = eduVec(i);
    if (eduStar == 1)
    {
      p_xast_y = 1-prob_X_ast_YZ(i);
    }
    else
    {
      p_xast_y = prob_X_ast_YZ(i);
    }
    
    p_xast_x = pMat(eduStar-1, 0);
    p_x_y = prob_x_yz(i, 0);
    rawProb = p_xast_x / p_xast_y *p_x_y;
    
    if (rawProb < 0)
    {
      rawProb = 0;
    }
    else if(rawProb > 1)
    {
      rawProb = 1;
    }
    
    outMat(i,0) = rawProb;
    outMat(i,1) = 1-rawProb;
  }
  
  return outMat;
}

// [[Rcpp::export]]
double Compute_Conditional_U_Cpp(NumericVector betaVal, NumericMatrix pMat,
                                 NumericMatrix pMatInv, NumericVector prob_X_ast_YZ,
                                 NumericVector eduVec, NumericVector ageVec,
                                 NumericVector yVec, NumericVector genderVec)
{
  int i,j;
  double _sum = 0, tmp;
  double age, gender, yVal;
  NumericMatrix prob_x_yz = Compute_X_Given_YZ_Cpp(pMatInv, prob_X_ast_YZ);
  NumericMatrix condProb = Compute_X_Given_X_star_YZ_Cpp(prob_X_ast_YZ, prob_x_yz, pMat, eduVec);
  NumericMatrix uMat(2,4);
  int num_of_obs = eduVec.length();
  NumericMatrix outU(num_of_obs, 4);
  
  for(i=0; i<num_of_obs; i++)
  {
    age = ageVec(i);
    gender = genderVec(i);
    yVal = yVec(i);
    uMat = Compute_U_Mat_App_Cpp(betaVal, yVal, age, gender, 2);
    
    for(j=0; j<4; j++)
    {
      outU(i,j) = condProb(i,0)*uMat(0,j)+condProb(i,1)*uMat(1,j);
    }
  }
  
  for(i=0; i<4; i++)
  {
    tmp = 0;
    
    for(j=0; j<num_of_obs; j++)
    {
      tmp += outU(j,i);
    }
    
    tmp /= num_of_obs;
    _sum += tmp*tmp;
  }
  
  return _sum;
}
