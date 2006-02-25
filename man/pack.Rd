\name{pack}
\alias{pack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Profile analysis using Clustering and Kurtosis}
\description{
 For a data matrix, selects features with kurtosis values in a specified
 range. Optionally, it further selects features based on whether their
 profiles are mixtures of gaussians or not.
}
\usage{
 pack(data.m,kurt.range=c(-2,0),cluster=T,method=c("bic","vb"));
}
\arguments{
  \item{data.m}{Data matrix with features along the rows.}
  \item{kurt.range}{Desired range of kurtosis values.}
  \item{cluster}{Logical, to indicate whether additional cluster
    learning step is desired.}
  \item{method}{Character specifying model selection to be used
    (bic=EM-algorithm + BIC, vb=variational Bayesian + evidence bound).}
}

\value{ A list with the following components:
    \item{out}{Matrix with rows labeling selected features and columns
      labeling kurtosis, cluster size and index position in data.m}
    \item{class}{A list with non-null elements giving the clustering
      classification of the selected features.}
  }
  
\author{Andrew Teschendorff \email{aet21@hutchison-mrc.cam.ac.uk}}


\keyword{cluster}% at least one, from doc/KEYWORDS

