\name{MembProbFn2}
\alias{MembProbFn2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal function for: Variational Bayesian Gaussian Mixture Model}
\description{
  Computes cluster membership probabilities of of samples.
}
\usage{
MembProbFn2(data, NewVals, Nsamples)
}

\arguments{
  \item{data}{The data matrix}
  \item{NewVals}{Estimated parameter values}
  \item{Nsamples}{Number of samples}
}

\value{ A list with the following components:
    \item{wcl}{Integer vector of length Nsamples specifying cluster
      membership of sample using maximum probabilty criterion}
    \item{probs}{Matrix of dimension Ncat x Nsamples giving cluster
      membership probabilities of samples.}
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

