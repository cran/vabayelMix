\name{vabayelMix}
\alias{vabayelMix}
\title{Variational Bayesian Gaussian Mixture Model}

\description{Learns a gaussian mixture model from data using an optimal
  separable approximation to the posterior density. The optimisation
  uses a variational procedure and implements an iterative ensemble
  learning algorithm. The algorithm gives a framework in which to infer
  the number of clusters in the data set. Prior information may be
  incorporated through specification of hyperparameters in a prior
  distribution. Current version implements a gaussian mixture model where
  the covariances matrices are diagonal.}

\usage{
vabayelMix(data, prior=NA, Ncat, nruns=100, npick=1, MaxIt=500,
conv.tol=0.001, nCVconv=10, verbatim=TRUE)
}

\arguments{
  \item{data}{A matrix of dimension Ns x Ndim containing the data to be clustered. Algorithm clusters rows of matrix and treats columns as dimensions.}
  \item{prior}{A list of various elements containing prior
    information as obtained for example by using \code{UseBasicPrior}. List elements are \code{prior$mean}, \code{prior$ivarm}, \code{prior$ivara},
    \code{prior$ivarb} and \code{prior$dapi}. The first four are matrices of dimension Ncat
  x Ndim, \code{prior$dapi} is a vector of length Ncat. \code{prior$mean} contains the
  means of the cluster mean gaussian priors. \code{prior$ivarm}
contains the inverse variances for the cluster mean gaussian priors. \code{prior$ivara} and
\code{prior$ivarb} contain the parameters for the gamma prior distribution of
the inverse variances of the clusters. \code{prior$dapi} is a weight vector
specifying prior knowledge about the number of clusters. If \code{prior}
is unspecified a complete uninformative prior is implemented that assumes
rows to be mean normalised to zero.}
  \item{Ncat}{The maximum number of clusters or categories to look for
    in the data set. Algorithm switches off clusters it doesn't
    need. See References.}
  \item{nruns}{Number of ensemble learning optimisation runs to be
    performed. Each optimisation run uses a different (random) starting
    point.}
  \item{npick}{The npick runs (out of nruns) that best optimise the cost
    function. See References.}
  \item{MaxIt}{Maximum number of iterations to be performed for a single
    optimisation run.}
  \item{conv.tol}{Threshold tolerance level for establishing
    convergence of iterations.}
  \item{nCVconv}{Number of consecutive iterations to consider in
    establishing convergence of the run at level \code{conv.tol}.}
  \item{verbatim}{Logical. If true prints out estimates and cost
    function value per iteration.}
}


 
\value{ A list with the following components:
  \item{estvals}{A list with components:
    \item{mean}{Means of gaussian posterior. Matrix of dimension Ncat x
      Ndim. A row containing all zeros means that component is absent.}
    \item{ivarm}{Inverse variances of gaussian posterior. Matrix of
      dimension Ncat x Ndim.}
    \item{ivara,ivarb}{Parameters of gamma posterior. Matrices of
      dimension Ncat x Ndim.}
    \item{dapi}{Parameters of dirichlet posterior giving weights of
      components. A value of 1 means that component is absent.}
  }
  \item{wcl}{A matrix of dimension npick x Ns. Each row gives cluster
    assignment of each row of data. Clusters are labeled by integers.
  }
  \item{probs}{A list of length npick, each list element is a matrix of
    dimension Ns x Ncat containing the probabilities of membership to
    clusters.}
  \item{costs}{A vector of length nruns specifying converged values of
    cost function.}
  \item{conv}{A binary vector of length nruns specifying if that run
    converged (0) or not (1).}
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

\examples{
NsTot <- 100; 
Nspg <- 50;
Ng <- 2;
deg.idx <- 1 ;
data <- matrix( nrow=NsTot, ncol=Ng);
for( s in 1:Nspg ){
  data[s,] <- rnorm(Ng,0,0.25);
}
for( s in (Nspg+1):NsTot){
  data[s,] <- rnorm(Ng,0,0.25);
  data[s,deg.idx] <- rnorm(1,2,0.25);
}
types.idx <- c(rep(1,50),rep(2,50));
useprior.l <- UseBasicPrior(data,rep(1,4));
vbmix <- vabayelMix(data, prior=NA, Ncat=4, nruns=10, npick=2,MaxIt=500, conv.tol=0.001, nCVconv=10);
# or could use
# vbmix <- vabayelMix(data, prior=useprior.l, Ncat=4, nruns=10, npick=2,MaxIt=500, conv.tol=0.001, nCVconv=10);
plot(1:NsTot,vbmix$wcl[1,],type="h",col=types.idx);
}


\keyword{cluster}

