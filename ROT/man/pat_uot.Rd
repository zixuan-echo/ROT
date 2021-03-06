% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pat_uot}
\alias{pat_uot}
\title{computer the transport plan of partial optimal transport}
\usage{
pat_uot(a, b, dab, reg, numItemax = 1000L, stopThr = 1e-20)
}
\arguments{
\item{a}{the the vector of source distribution}

\item{b}{the vector of target dirtribution}

\item{reg}{entropic regularization parameter}

\item{numItemax}{the max iteration times}

\item{stopThr}{ending condition}

\item{mat}{the cost matrix}
}
\description{
computer the transport plan of partial optimal transport
}
\examples{
a<-c(0.4,0.8)
b<-c(0.6,3)
m1<-c(0,1)
m2<-c(1,0)
m<-cbind(m1,m2)
pat_uot(a,b,m,1)
}
