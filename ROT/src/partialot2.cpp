#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

//'@title computer the transport plan of partial optimal transport
//'@param a the the vector of source distribution
//'@param b the vector of target dirtribution 
//'@param mat the cost matrix
//'@param reg entropic regularization parameter
//'@param numItemax the max iteration times
//'@param stopThr ending condition
//'@examples
//'a<-c(0.4,0.8)
//'b<-c(0.6,3)
//' m1<-c(0,1)
//' m2<-c(1,0)
//' m<-cbind(m1,m2)
//' pat_uot(a,b,m,1)
// [[Rcpp::export]]
arma::mat pat_uot(arma::vec a,arma::vec b,arma::mat dab,double reg,
         int numItemax=1000,double stopThr=1e-20){
     int dim_a = dab.n_rows;
     int dim_b = dab.n_cols;
     vec sum_k(dim_a);
     vec sum_k1(dim_b);
     vec dx(dim_a);
     dx.fill(1.0);
     vec dy(dim_b);
     dy.fill(1.0);
    //double m=0.1;
    double m;
     if(accu(a)>accu(b))
       m=accu(b);
     else 
       m=accu(a);
     mat K(dim_a, dim_b);
     mat K1(dim_a,dim_b);
     mat K2(dim_a,dim_b);
     K = dab /(-reg);
     K = exp(K); 
     double c=accu(K);
     K=(m/c)*K;
     double cpt=0;
     double error=1;
     mat Kprep(dim_a,dim_b);
     
while(error>stopThr and cpt<numItemax){
     Kprep=K;
     vec sumk=(sum(K,1));
     sum_k=a/sumk;
      for(int i=0;i<dim_a;i++){
      if(sum_k[i]>dx[i])
      sum_k[i]=dx[i];
     }
      K1=diagmat(sum_k)*K;
      vec sumk1=(sum(K1,0)).t();
      sum_k1=b/sumk1;
      for(int j=0;j<dim_b;j++){
      if(sum_k1[j]>dy[j])
      sum_k1[j]=dy[j];
      }
      K2=diagmat(sum_k1)*K1;
      K=K2*(m/accu(K2));
      cpt++;
     }
     return Kprep;
   }


//'@title computer the wasserstein distance  of partial optimal transport
//'@param a the the vector of source distribution
//'@param b the vector of target dirtribution 
//'@param mat the cost matrix
//'@param reg entropic regularization parameter
//'@examples
//'a<-c(0.4,0.8)
//'b<-c(0.6,3)
//' m1<-c(0,1)
//' m2<-c(1,0)
//' m<-cbind(m1,m2)
//' pat_uot_was(a,b,m,1)
// [[Rcpp::export]]
double pat_uot_was(arma::vec a,arma::vec b,arma::mat dab,double reg){
   double was;
   int dim_a = dab.n_rows;
   int dim_b = dab.n_cols;
   mat wasser(dim_a,dim_b);
   wasser=pat_uot(a,b,dab,reg)%dab;
   was=accu(wasser);
   return was;
}
   