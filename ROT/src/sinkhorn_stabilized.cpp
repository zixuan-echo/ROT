#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;


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
//' sinkhorn_stabilized(a,b,m,1,1)
//[[Rcpp::export]]
arma::mat sinkhorn_stabilized(arma::vec a,arma::vec b,arma::mat dab,
    double reg,double reg_m,double tau=1e5,
    int numltermax=1000,double stopThr=1e-6){ //bool verbose,bool log
  int dim_a = dab.n_rows;
  int dim_b = dab.n_cols;
  //if(a.size()==0){
    
  //}
  //if(b.size()==0){
    
  //}
  vec u(dim_a);
  u.fill(1.0/dim_a);
  vec v(dim_b);
  v.fill(1.0/dim_b);
  mat K(dim_a, dim_b);
  K = dab / (-reg);
  K = exp(K); 
  double fi=reg_m/(reg_m+reg);
  double cpt=0;
  double error=1;
  vec alpha(dim_a);
  alpha.fill(0);
  vec beta(dim_b);
  beta.fill(0);
  vec uprev;
  vec vprev;
  vec Kv;
  vec Ktu;
  while(error>stopThr and cpt<numltermax){
    uprev=u;
    vprev=v;
    Kv = K * v;
    vec f_alpha = exp(- alpha / (reg + reg_m));
    vec f_beta = exp(- beta / (reg + reg_m));
    u = pow((a / (Kv + 1e-16)),fi) % f_alpha;
    Ktu = K.t()*u;
    v = pow((b / (Ktu + 1e-16)),fi)% f_beta;
    bool absorbing = 0;
    if(arma::any(u>tau) or arma::any(v>tau)){
    absorbing=1; 
    vec logu(dim_a);
    vec logv(dim_b);
    logu.fill(log(max(u)));
    logv.fill(log(max(v)));
    alpha=alpha+reg*logu;
    beta=beta+reg*logv;
    mat Alpha(dim_a, dim_b, fill::zeros);
    mat Beta(dim_a, dim_b, fill::zeros);
    Alpha.each_col()=alpha;
    Beta.each_row()= beta.t();
    K=exp((Alpha+Beta-dab)/reg);
    v=v.fill(1);
    }
    Kv=K*v;
    cpt=cpt+1;
      
      //if n_hists:
      //f_alpha = f_alpha[:, None]
      //f_beta = f_beta
    
  }
  vec log_u=alpha/reg+log(u);
  vec log_v=beta/reg+log(v);
  mat log_U(dim_a, dim_b, fill::zeros);
  mat log_V(dim_a, dim_b, fill::zeros);
  log_U.each_col() =log_u;
  log_V.each_row()= log_v.t();
  mat ot_mat(dim_a,dim_b);
  ot_mat=exp((log_U+log_V-dab)/reg);
  return ot_mat;
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
//'sinkhorn_stabilized_was(a,b,m,1,1)
//'
//'@references
//'[1] Chizat, L., Peyre, G., Schmitzer, B., & Vialard, F. X. (2016). Scaling algorithms for unbalanced transport problems. arXiv preprint arXiv:1607.05816.
//'[2] Frogner C., Zhang C., Mobahi H., Araya-Polo M., Poggio T. Learning with a Wasserstein Loss,  Advances in Neural Information Processing Systems (NIPS) 2015.
// [[Rcpp::export]]
double sinkhorn_stabilized_was(arma::vec a,arma::vec b,arma::mat dab,
                          double reg,double reg_m){
  double was;
  int dim_a = dab.n_rows;
  int dim_b = dab.n_cols;
  mat wasser(dim_a,dim_b);
  wasser=sinkhorn_stabilized(a,b,dab,reg,reg_m)%dab;
  was=accu(wasser);
  return was;
}