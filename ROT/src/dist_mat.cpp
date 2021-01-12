#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
//' @title distance between two samples(dimension=1)
//' @description Compute the distance matrix between samples of one dimention.
//' @param a sample a
//' @param b sample b
//' @param k the order of the norm (default k=1)
//' @return the kth power of the k-norm of the distance between a and b
//' @examples
//' distance=dists(1:100,1:100)
// [[Rcpp::export]]
arma::mat dists(arma::vec a,arma::vec b,int k=1,int dimention=1){
  int dim_a=a.size();
  int dim_b=b.size();
  mat dists(dim_a,dim_b);
  mat source(dim_a,dim_b,fill::zeros);
  mat target(dim_a,dim_b,fill::zeros);
  source.each_col()=a;
  target.each_row()=b.t();
  if(k==1){
    dists=abs(source-target);
  }
  if(k>1){
    dists=pow((source-target),k);
  }
  return dists;
}

//' @title distance between two samples(dimension>1)
//' @description Compute the distance matrix between two samples (dimensions>1)
//' @param a mat a(n,p) n:the number of samples,p:the dimensionnality
//' @param b mat b(n,p) n:the number of samples,p:the dimensionnality
//' @param k the order of the norm of the distance (default k=1)
//' @return the k-norm of the distance between a and b
//' @examples
//' a = cbind(1:10,1:10)
//' dist_mat(a,a)
// [[Rcpp::export]]
arma::mat dist_mat(arma::mat a,arma::mat b,int k=1,int square=0){
  int simple_a=a.n_rows;
  int simple_b=b.n_rows;
  int dimention=a.n_cols;
  mat dist(simple_a,simple_b,fill::zeros);
  mat tmp(simple_a,simple_b,fill::zeros);
  mat result(simple_a,simple_b,fill::zeros);
  int i=0;
  if(square==1) k=2;
  for(i=0;i<dimention;i++){
    vec v_a=a.col(i);
    vec v_b=b.col(i);
    mat distss(vec a,vec b,int k=1,int dimention=1);
    tmp=dists(v_a,v_b,k);
    dist=dist+tmp;
  }
  if(square==1) {
    result=dist;
    return result;}
  result=pow(dist,1/k);
  return result;
}

