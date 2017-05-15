#include <vector>
#include <map>
#include <algorithm>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
arma::colvec PredictCpp(const arma::mat &x, const arma::mat &newx, const arma::colvec &y){

  int  n=x.n_rows, np = newx.n_rows, p=x.n_cols;
  double eps = 1e-12;
  double h=mean(mean(stddev(x,0,0)))/(pow(n,1/(p+4.0)));

  arma::colvec yp(np);
  for(int i=0;i<np;++i){
    arma::colvec dxi = sum(pow(x.each_row()-newx.row(i),2),1)/(2*h*h);
    arma::colvec k = exp(-dxi);
    if(sum(k)>eps){
      yp[i] = dot(k,y)/sum(k);
    }
    else{
      yp[i] = mean(y);
    }
  }

  return yp;
}
