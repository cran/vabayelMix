\name{UseBasicPrior}
\alias{UseBasicPrior}

\title{Prior Function for Variational Gaussian Mixture Model}
\description{
  This function implements an uninformative prior distribution for the
  cluster centers and variances, but allows the user to define prior
  weights for the clusters.
}
\usage{
UseBasicPrior(data, weights.v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{ A matrix with columns representing variables
    and rows observations. Algorithm clusters observations.}
  \item{weights.v}{ A vector of relative prior weights for the clusters. }
}
\details{
 \code{weights.v} is a vector of length \code{Ncat}, the maximum number
 of clusters to look for. 
}
\value{A list with following components. The first four are matrices of dimension \code{Ncat} x \code{Ndim}, \code{dapi} is a vector of length \code{Ncat}.
  \item{mean}{the means of the cluster mean gaussian priors.}
  \item{varm}{the inverse variances for the cluster mean gaussian priors.}
  \item{ivara}{parameters for the gamma prior distribution of the inverse variances of the clusters. See references.}
  \item{ivarb}{parameters for the gamma prior distribution of the inverse variances of the clusters. See references.}
  \item{dapi}{weight vector specifying prior knowledge about the number of clusters.}
}

\references{
  \item{1}{D.J.MacKay: Developments in probabilistic modelling with neural
    networks-ensemble learning. In Neural Networks: Artificial
    Intelligence and Industrial Applications. Proceedings of the 3rd
    Annual Symposium on Neural Networksm Nijmengen, Netherlands, Berlin
    Springer, 191-198 (1995).}
  \item{2}{J.W.Miskin : Ensemble Learning for Independent Component Analysis,
    PhD thesis University of Cambridge December 2000.}
  \item{3}{A. E. Teschendorff,...et al.: A variational bayesian mixture
    modelling framework for cluster analysis of gene expression
    data. Submitted to Bioinformatics.}
}

\author{Andrew Teschendorff\email{aet21@hutchison-mrc.cam.ac.uk}}

\keyword{cluster}% at least one, from doc/KEYWORDS

