#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

//'@title computer the plan of regularized unbalanced optimal transport
//'@description Solve the entropic regularization unbalanced optimal transport
//'@param a vec a is the the histogram of source distribution(dim_a=m)
//'@param b vec b is the histogram of target dirtribution (dim_b=n)
//'@param dab matrix dab is the cost matrix between a and b
//'@param reg entropic regularization parameter
//'@param reg_m marginal relaxation term > 0
//'@param numItemax the max iteration times
//'@param stopThr ending condition
//'@return a matrix that shape:m
//'@examples
//'a<-c(0.4,0.8)
//'b<-c(0.6,3)
//' m1<-c(0,1)
//' m2<-c(1,0)
//' m<-cbind(m1,m2)
//' sinkhorn_knopp_unbalanced(a,b,m,1,1)
// [[Rcpp::export]]
arma::mat sinkhorn_knopp_unbalanced(arma::vec a, arma::vec b, arma::mat dab, double reg, double reg_m, int numItermax=100, double stopThr=1e-6){
  int dim_a = dab.n_rows;
  int dim_b = dab.n_cols;
  arma::vec u(dim_a);
  u.fill(1.0/dim_a);
  arma::vec v(dim_b);
  v.fill(1.0/dim_b);
  
  arma::mat K(dim_a, dim_b);
  K = dab / (-reg);
  K = exp(K);

  double fi = reg_m / (reg_m + reg);
  double err = 1.0;
  
  arma::vec uprev;
  arma::vec vprev;
  arma::vec Kv;
  arma::vec Ktu;
  
  for(int i=0; i<numItermax; i++){
    uprev = u;
    vprev = v;
    Kv = K * v;
    u = arma::pow((a / Kv), fi);
    Ktu = K.t() * u;
    v = arma::pow((b / Ktu), fi);

    if(arma::any(Ktu == 0.0) or u.has_nan() or v.has_nan() or u.has_inf() or v.has_inf()){
      cout<<"Numerical errors at iteration"<<i<<endl;
      u = uprev;
      v = vprev;
      break;
    }
  }
  double m1 = (arma::max(arma::abs(u))>arma::max(arma::abs(uprev)))?arma::max(arma::abs(u)):arma::max(arma::abs(uprev));
  m1 = m1>1?m1:1;
  double m2 = (arma::max(arma::abs(v))>arma::max(arma::abs(vprev)))?arma::max(arma::abs(v)):arma::max(arma::abs(vprev));
  m2 = m2>1?m2:1;
  double err_u = arma::max(arma::abs(u - uprev)) / m1;
  double err_v = arma::max(arma::abs(v - vprev)) / m2;
  err = 0.5 * (err_u + err_v);
  mat U(dim_a, dim_b, fill::zeros);
  mat V(dim_a, dim_b, fill::zeros);
  U.each_col() = u;
  V.each_row() = v.t();
  mat res(dim_a, dim_b);
  res = U % K % V;
  return res;
}

//'@title sinkhorn unbalanced OT distance
//'@description Solve the entropic regularization unbalanced optimal transport problem and return the loss.
//'
//'@param a vector(dim_a). Unnormalized histogram of dimension dim_a.
//'@param b vector(dim_b). Unnormalized histogram of dimension dim_b.
//'@param dab matrix(dim_a, dim_b). loss matrix.
//'@param reg float. Entropy regularization term > 0.
//'@param reg_m float. Marginal relaxation term > 0.
//'
//'@return the OT distance between `a` and 'b`
//'
//'@examples
//'a<-c(0.4,0.8)
//'b<-c(0.6,3)
//'m1<-c(0,1)
//'m2<-c(1,0)
//'m<-cbind(m1,m2)
//'sinkhorn_knopp_unbalanced_was(a,b,m,1,1)
//'
//'@references
//'[1] Chizat, L., Peyre, G., Schmitzer, B., & Vialard, F. X. (2016). Scaling algorithms for unbalanced transport problems. arXiv preprint arXiv:1607.05816.
//'[2] Frogner C., Zhang C., Mobahi H., Araya-Polo M., Poggio T. Learning with a Wasserstein Loss,  Advances in Neural Information Processing Systems (NIPS) 2015.
// [[Rcpp::export]]
double sinkhorn_knopp_unbalanced_was(arma::vec a,arma::vec b,arma::mat dab,
                          double reg,double reg_m){
  double was;
  int dim_a = dab.n_rows;
  int dim_b = dab.n_cols;
  mat wasser(dim_a,dim_b);
  wasser=sinkhorn_knopp_unbalanced(a,b,dab,reg,reg_m)%dab;
  was=accu(wasser);
  return was;
}