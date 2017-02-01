
#include "MAVEfast.h"
#include "CVfast.h"
#include "CVfast_emxutil.h"
#include <vector>
#include <map>
#include <algorithm>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
List MAVEfastCpp(NumericVector x,NumericVector y,CharacterVector method) {

  int lenx=x.size(),leny=y.size();
  int nrow = leny, ncol = lenx/leny;
  //printf("%d %d\n",nrow,ncol);
  emxArray_real_T *emxArray_x,*emxArray_y,*emxArray_BB,*emxArray_ky;
  emxArray_char_T *emxArray_method;
  std::string std_method = Rcpp::as<std::string>(method);
  emxInit_real_T(&emxArray_x,2);
  emxInit_real_T(&emxArray_y,2);
  emxInit_real_T(&emxArray_BB,2);
  emxInit_real_T(&emxArray_ky,2);
  emxInit_char_T(&emxArray_method,2);
  emxArray_x->size[0]=nrow;
  emxArray_x->size[1]=ncol;
  emxArray_y->size[0]=nrow;
  emxArray_y->size[1]=1;
  emxArray_BB->size[0]=pow(ncol,3);
  emxArray_BB->size[1]=1;
  //emxArray_BB->size[2]=ncol;
  emxArray_method->size[0]=1;
  emxArray_method->size[1] = std_method.size();
  emxEnsureCapacity((emxArray__common *)emxArray_x, 0, (int)sizeof(double));
  emxEnsureCapacity((emxArray__common *)emxArray_y, 0, (int)sizeof(double));
  emxEnsureCapacity((emxArray__common *)emxArray_BB, 0, (int)sizeof(double));
  emxEnsureCapacity((emxArray__common *)emxArray_method, 0, (int)sizeof(char));

  for(int i=0;i<nrow;++i){
    for(int j=0;j<ncol;++j){
      emxArray_x->data[j*emxArray_x->size[0]+i] = x[j*emxArray_x->size[0]+i];
    }
  }
  //for(int i=0;i<5;++i) printf("%.3lf ",x[i]); printf("\n");
  for(int i=0;i<nrow;++i){
    emxArray_y->data[i] = y[i];
  }
  for(int i=0;i<std_method.size();++i){
    emxArray_method->data[i] = std_method[i];
  }
  MAVEfast(emxArray_x, emxArray_y, emxArray_method, emxArray_BB, emxArray_ky);
  NumericVector BB(pow(ncol,3));
  NumericVector ky(emxArray_ky->size[0]*emxArray_ky->size[1]);
  NumericVector nx(lenx);

  for(int i=0;i<lenx;++i){
    nx[i]=emxArray_x->data[i];
  }

  for(int i=0;i<pow(ncol,3);++i){
    BB[i]=emxArray_BB->data[i];
  }

  for(int i=0;i<emxArray_ky->size[0];++i){
    for(int j=0;j<emxArray_ky->size[1];++j){
      ky[j+i*emxArray_ky->size[1]] = emxArray_ky->data[j+i*emxArray_ky->size[1]];
    }
  }

  IntegerVector dimBB(3);
  for(int i=0;i<3;++i){
    dimBB[i] = ncol;
  }
  BB.attr("dim") = dimBB;

  IntegerVector dimky(2);
  dimky[0] = emxArray_ky->size[0];
  dimky[1] = emxArray_ky->size[1];
  ky.attr("dim") = dimky;

  IntegerVector dimnx(2);
  dimnx[0] = nrow;
  dimnx[1] = ncol;
  nx.attr("dim")=dimnx;
  Rcpp::List result = Rcpp::List::create(Rcpp::Named("BB")=BB,
                                         Rcpp::Named("ky")=ky,
                                         Rcpp::Named("x")=nx);
  result.attr("class")="mave";
  return result;

}
