#include <vector>
#include <map>
#include <algorithm>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

arma::mat tiedrank(arma::mat x){
  int nrow=x.n_rows,ncol=x.n_cols;
  arma::mat x_sorted=sort(x);
  arma::mat t_rank;t_rank.zeros(nrow,ncol);
  arma::mat rank;
  rank.zeros(x.n_rows,x.n_cols);
  for(int i=0;i<ncol;++i){

    map<double, double> val2t_rank;
    double flag=x_sorted(0,i);
    int ppos=0;
    for(int j=1;j<nrow;++j){
      if((abs(x_sorted(j,i)-flag)>1e-12)){
        for(int k=ppos;k<j;++k){
          rank(k,i)=(ppos+j-1.0)/2;
          val2t_rank.insert(make_pair(x_sorted(k,i),(ppos+j-1.0)/2));
        }
        flag=x_sorted(j,i);
        ppos=j;
      }
    }
    int j = x.n_rows;
    for(int k=ppos;k<j;++k){
      rank(k,i)=(ppos+j-1.0)/2;
      val2t_rank.insert(make_pair(x_sorted(k,i),(ppos+j-1.0)/2));
    }

    for(int j=0;j<nrow;++j){
      t_rank(j,i)=val2t_rank[x(j,i)];
    }
  }

  return t_rank;
}


bool cmp(pair<double,int>first,pair<double,int>second){
  if(first.first<second.first)
    return false;
  else
    return true;
}


arma::uvec max_num(arma::colvec x,int num) {
  int n=x.size();
  vector<pair<double,int> >win;
  for(int i=0;i<num;++i){
    win.push_back(make_pair(x[i],i));
  }

  make_heap(win.begin(),win.end(),cmp);

  for(int i=num;i<n;++i){

    if(x[i]>win[0].first){
      pop_heap(win.begin(),win.end(),cmp);
      win.pop_back();
      win.push_back(make_pair(x[i],i));
      push_heap(win.begin(),win.end(),cmp);
    }
  }
  arma::uvec maxIdx;
  maxIdx.resize(num);
  for(int i=0;i<num;++i){
    maxIdx[i]=win[i].second;
  }
  return maxIdx;
}


//[[Rcpp::export]]
arma::rowvec CVfast(arma::mat x,arma::mat ky,arma::colvec BB1D){
  int n = x.n_rows;
  int ny=ky.n_cols;
  int nD=x.n_cols;
  arma::rowvec cv;cv.zeros(nD);

  for(int p=0;p<nD;++p){

    arma::mat nBB;nBB.zeros(nD,p+1);
    for(int i=0;i<nBB.size();++i){
      nBB[i]=BB1D[pow(nD,2)*p+i];
    }
    arma::mat nx=x*nBB;
    nx=tiedrank(nx)/(n+1);
    double h=mean(mean(stddev(nx,0,0)))/(pow(n,1/(p+5.0)));
    double h2 = 2*h*h;

    for(int i=0;i<n;++i){
      arma::colvec dxi=sum(pow(nx-repmat(nx.row(i),n,1),2),1)/h2;
      arma::colvec k=(1-dxi)%(1>dxi);
      arma::uvec idx=max_num(k,2);
      arma::rowvec val;val.zeros(2);
      k.elem(idx)=val;
      k=k/(sum(k)+1e-6);
      arma::rowvec ye=k.t()*ky;
      cv[p]+=mean(abs(ye-ky.row(i)))/n;
    }
  }
  return cv;
}
