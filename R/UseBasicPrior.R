"UseBasicPrior" <-
function(data,weights.v){

Ndim <- dim(data)[2];
Ncat <- length(weights.v);

# for means assume prior gaussian of means
m0 <- matrix(rep(0,Ncat*Ndim), nrow=Ncat, ncol=Ndim);
for( c in 1:Ncat){
 for( d in 1:Ndim){
  m0[c,d] <- runif(1,min=min(data[,d]),max=max(data[,d]));
 }
}
# and inv.variances 
am0 <- matrix(rep(0.001,Ncat*Ndim),nrow=Ncat, ncol=Ndim);

for( d in 1:Ndim){
 1/var(data[,d]) -> est.ivar;
 am0[,d] <- rep(est.ivar,Ncat);
}
 
# for inv. variances assume gamma distribution with parameters
aiv0 <- matrix(rep(0.001,Ncat*Ndim),nrow=Ncat, ncol=Ndim);
biv0 <- matrix(rep(0.001,Ncat*Ndim),nrow=Ncat, ncol=Ndim);

for( d in 1:Ndim){
 1/var(data[,d]) -> est.ivar;
 biv0[,d] <- est.ivar*aiv0[,d];
}

# parameter of dirichlet posterior
api0 <- weights.v ;
#################################################################
useprior.l <- list(mean=m0, ivarm=am0, ivara=aiv0, ivarb=biv0, dapi=api0);

return(useprior.l);

} # END OF FUNCTION

