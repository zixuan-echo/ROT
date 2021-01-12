\name{ROT-package}
\alias{ROT-package}
\alias{ROT}
\docType{package}
\title{
\packageTitle{ROT}
}
\description{
\packageDescription{ROT}
}
\details{
This open source R package provides several solvers for  optimization problems related to unbalanced Optimal Transport.
}
\author{
\packageAuthor{ROT}

Maintainer: \packageMaintainer{ROT}
}
\references{
[1]Chizat, Lenaic, Peyre, Gabriel, Schmitzer, B. , & Vialard, Francois-Xavier. (2018). Scaling algorithms for unbalanced optimal transport problems. Mathematics of Computation, 87(314), 2563-2609.\cr
[2]Benamou, J. D. , Carlier, G. , Cuturi, M. , Nenna, L. , & Peyre, Gabriel. (2014). Iterative bregman projections for regularized transportation problems. Siam Journal on Scientific Computing, 37(2), A1111-A1138.\cr
[3]Peyre, Gabriel, Marco Cuturi, and Justin Solomon, Gromov-Wasserstein averaging of kernel and distance matrices. International Conference on Machine Learning (ICML).2016.\cr
[4]Chapel, L., Alaya, M., Gasso, G. (2020). Partial Optimal Transport with Applications on Positive-Unlabeled Learning. NeurIPS.
}

\keyword{ package }

\examples{
n=10
a=dnorm(1:100,20,10)
a=a/sum(a)
b=dnorm(1:100, mean=60, sd=2)
b=b/sum(b)*2
distance= dists(1:100,1:100)
distance= diatance/max(distance)
plan=sinkhorn_knopp_unbalanced(a,b,distance,1,1,1000,1e-6,TRUE,FALSE)
plot1D(a, b, plan, title = "UOT example")
}
