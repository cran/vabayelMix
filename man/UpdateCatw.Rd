\name{UpdateCatw}
\alias{UpdateCatw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal function for: Variational Bayesian Gaussian Mixture Model}
\description{
  Updates categorical weights
}
\usage{
UpdateCatw(Ncat, data, m, am, aiv, biv, api)
}
\arguments{
  \item{m,am,aiv,biv,api}{Posterior parameters.}
}

\value{ A list with the following components:
    \item{cwm}{Categorical weight matrix. See References}
    \item{scw}{Derive from above. See References.}
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

\keyword{internal}% __ONLY ONE__ keyword per line
