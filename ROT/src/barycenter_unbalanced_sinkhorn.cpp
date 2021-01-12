#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

//'@title unbalanced barycenter
//'@description Compute the entropic unbalanced wasserstein barycenter of A.
//'
//'@param A matrix(dim, n). n training distributions a_i of dimension dim.
//'@param M matrix(dim, dim). ground metric matrix for OT.
//'@param weighs vector(n). Weight of each distribution (barycentric coodinates).
//'@param reg Entropy regularization term > 0.
//'@param reg_m Marginal relaxation term > 0.
//'@param numItermax int. Max number of iterations.
//'@param stopThr float. Stop threshol on error (> 0).
//'@param verbose bool. Print information along iterations.
//'
//'@return vector(dim). Unbalanced Wasserstein barycenter
//'
//'@examples
//'a=dnorm(1:100,20,10)
//'a=a/sum(a)
//'b=dnorm(1:100, mean=60, sd=2)
//'b=b/sum(b)*2
//'dab = matrix(c(a, b), 100, 2)
//'distance= dists(1:100,1:100)
//'distance = distance/max(distance)
//'res = barycenter_unbalanced_stabilized(dab, distance, reg = 1, reg_m = 1, weights = c(0.5, 0.5))
//'p = data.frame(ind = 1:100, a = a, b = b, c = res)
//'ggplot(data = p)+geom_line(aes(x = ind, y = a), color = "blue")+
//'geom_line(aes(x = ind, y = b), color = "red")+
//'geom_line(aes(x = ind, y = res), color = "green")
//'
//'@references
//'[1] Benamou, J. D., Carlier, G., Cuturi, M., Nenna, L., & Peyre, G.(2015). Iterative Bregman projections for regularized transportation problems. SIAM Journal on Scientific Computing, 37(2), A1111-A1138.
//'[2] Chizat, L., Peyre, G., Schmitzer, B., & Vialard, F. X. (2016). Scaling algorithms for unbalanced transport problems. arXiv preprin arXiv:1607.05816.
//'
// [[Rcpp::export]]
arma::mat barycenter_unbalanced_stabilized(arma::mat A, arma::mat M, arma::vec weights, double reg = 0.1, double reg_m = 1, int numItermax = 1000, double stopThr = 1e-8, bool verbose = false){
  
  int dim_a = A.n_rows;
  int dim_b = A.n_cols;
  arma::mat K(dim_a, dim_a);
  K = M / (-reg);
  K = arma::exp(K);
  double fi = reg_m / (reg_m + reg);
  double err = 1.0;
  
  arma::mat v(dim_a, dim_b, fill::ones);
  arma::mat u(dim_a, 1, fill::ones);
  arma::vec q(dim_a, fill::ones);
  arma::mat vprev(dim_a, dim_b, fill::ones);
  arma::mat uprev(dim_a, 1, fill::ones);
  arma::vec qprev(dim_a, fill::ones);
  arma::mat Q(dim_a, dim_b);
  arma::mat Kv(dim_a, dim_b);
  arma::mat Ktu(dim_a, 1);
  double m;
  for(int i=0; i<numItermax; i++){
    uprev = u;
    vprev = v;
    qprev = q;
    Kv = K * v;
    u = arma::pow((A / Kv), fi);
    Ktu = K.t() * u;
    q = (1 - fi) * Ktu * weights;
    q = arma::pow(q, 1 / (1 - fi));
    Q.each_col() = q;
    v = arma::pow((Q / Ktu), fi);
    if(Ktu.is_zero() or u.has_nan() or v.has_nan() or u.has_inf() or v.has_inf()){
      cout<<"Numerical errors at iteration"<<i<<endl;
      u = uprev;
      v = vprev;
      break;
    }
    m = (arma::max(arma::abs(q))>arma::max(arma::abs(qprev)))?arma::max(arma::abs(q)):arma::max(arma::abs(qprev));
    m = m>1?m:1;
    err = arma::max(arma::abs(q - qprev)) / m;
    if(err < stopThr && i > 10){
      break;
    }
    if(verbose){
      if(i % 10 == 0){
        cout<<"It.  Err    "<<endl;
      }
      cout<<i<<"   "<<err<<endl;
    }
  }
  
  return q;
}
