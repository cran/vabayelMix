\name{CostKL}
\alias{CostKL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal function for: Variational Bayesian Gaussian Mixture Model}
\description{
  Computes the value of the cost function. Used for monitoring convergence.
}
\usage{
CostKL(Ncat,data, m0, am0, aiv0, biv0, api0, m, am, aiv, biv, api, Catwm)
}
\arguments{
  \item{m0,am0,aiv0,biv0,api0}{Prior hyperparameters, see vabayelMix}.
  \item{m,am,aiv,biv}{Posterior parameters.}
  \item{lambda}{Categorical weight matrix, see References.}
  \item{s.lambda}{ Derived from lambda, see References.}
}

\value{ A list with the following components:
    \item{mean}{Means of gaussian posterior. Matrix of dimension Ncat x
      Ndim.}
    \item{ivarm}{Inverse variances of gaussian posterior. Matrix of
      dimension Ncat x Ndim.}
    \item{ivara,ivarb}{Parameters of gamma posterior. Matrices of
      dimension Ncat x Ndim.}
    \item{dapi}{Parameters of dirichlet posterior giving weights of
      components.}
  }
  
  \references{
  \item{1}{D.J.MacKay: Developments in probabilistic modelling with neural
    networks-ensemble learning. In Neural Networks: Artificial
    Intelligence and Industrial Applications. Proceedings of the 3rd
    Annual Symposium on Neural Networksm Nijmengen, Netherlands, Berlin
    Springer, 191-198 (1995).}
  \item{2}{J. W. Miskin : Ensemble Learning for Independent Component Analysis,
    PhD thesis University of Cambridge December 2000.}
  \item{3}{A. E. Teschendorff,...et al.: A variational bayesian mixture
    modelling framework for cluster analysis of gene expression
    data. Submitted to Bioinformatics.}
}

\author{Andrew Teschendorff\email{aet21@hutchison-mrc.cam.ac.uk}}


\keyword{internal}% at least one, from doc/KEYWORDS

