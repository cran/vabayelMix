"MembProbFn2" <-
function(data,NewVals,Nsamples){
print("Entering MembProbFn2: Computing probabilities");

# constants
Ncat <- dim(NewVals$mean)[1];
Ndim <- dim(NewVals$mean)[2];
print(paste("Ncat=",Ncat, " Ndim=",Ndim," Nsamples=",Nsamples,sep=""));

# input
m   <- NewVals$mean;
am  <- NewVals$ivarm;
aiv <- NewVals$ivara;
biv <- NewVals$ivarb;
api <- NewVals$dapi;
# variables to be computed
pcat <- matrix(rep(0,times=Ncat*Nsamples),nrow=Ncat,ncol=Nsamples);
wCl <- vector(length=Nsamples);

# Weight of each cluster
prob.cat <- api/sum(api) ;

# Find which cluster each data point belongs to
pcat <- matrix(0,nrow=Ncat,ncol=Nsamples);
wCl <- vector(length=Nsamples);
for ( sn in 1:Nsamples){
  # p(c|d_s) = p(c) p(d_s|c);
  for ( c in 1:Ncat){

  pcat[c,sn] <- prob.cat[c]*exp(-0.5*sum( (biv[c,]/aiv[c,])*(data[sn,]-m[c,])*(data[sn,]-m[c,])))*sqrt(prod(biv[c,]/aiv[c,]));

  }

  spcat <- sum( pcat[,sn] );
  pcat[,sn] <- pcat[,sn]/spcat ;

  wCl[sn] <- which.max(pcat[,sn]);
}


 print("Finished membprobFn2");
 return(list( wcl=wCl, probs=pcat));


} ### END OF FUNCTION
