"vabayelMix" <-
function(data, prior=NA, Ncat, nruns=100, npick=1, MaxIt=500, conv.tol=0.001, nCVconv=10, verbatim=TRUE){

if( is.na(Ncat) ){
  stop("Must specify =maximum number of categories/clusters to look for.");
}

# Basic constants
Nsamples <- dim(data)[1];
Ndim <- dim(data)[2];

# Prior hyperparameters (labeled by 0) 

if( is.na(prior) ){ # then use broad priors
  print("Using broad hyperparameters in priors");
# for means assume prior gaussian of means
m0 <- matrix(rep(0,Ncat*Ndim), nrow=Ncat, ncol=Ndim);
# and inv.variances 
am0 <- matrix(rep(0.001,Ncat*Ndim),nrow=Ncat, ncol=Ndim);
# for inv. variances assume gamma distribution with parameters
aiv0 <- matrix(rep(0.001,Ncat*Ndim),nrow=Ncat, ncol=Ndim);
biv0 <- matrix(rep(0.001,Ncat*Ndim),nrow=Ncat, ncol=Ndim);
# parameter of dirichlet posterior
api0 <- rep(1,times=Ncat);
}
else { # have specified priors so use it instead
  
  if(length(prior) != 5 ){
   stop("prior must be a list of length 5");
  }
  if( (ncol(prior$mean)==Ndim) && (nrow(prior$mean)==Ncat) ){
    m0 <- prior$mean ;
  }
  else {
    stop("Wrong dimensions for the mean of prior mean");
  }
  if( (ncol(prior$ivarm)==Ndim) && (nrow(prior$ivarm)==Ncat) ){
    am0 <- prior$ivarm ;
  }
  else {
    stop("Wrong dimensions for inv variance of prior mean");
  }
  if( (ncol(prior$ivara)==Ndim) && (nrow(prior$ivara)==Ncat) ){
    aiv0 <- prior$ivara ;
  }
  else {
    stop("Wrong dimensions for parameter a in prior for inv. variance");
  }
  if( (ncol(prior$ivarb)==Ndim) && (nrow(prior$ivarb)==Ncat) ){
    biv0 <- prior$ivarb ;
  }
  else {
    stop("Wrong dimensions for parameter b in prior for inv. variance");
  }
  if( length(prior$dapi)==Ncat ){
    api0 <- prior$dapi ;
  }
  else {
    stop("Wrong dimensions for dirichlet prior ");
  }

}

  
##############################################################################

### VARIABLES NEEDED #################################################
wCl <- matrix(nrow=npick,ncol=Nsamples);
v.Costs <- vector() ; v.NC <- vector();
l.NewVals <- list(NULL);
#########################################################################

### START RUNS ##########################################################
for ( run in 1:nruns){

# ENSEMBLE INITIALISATION (only randomise means)#########################
m <- matrix(rep(0,Ncat*Ndim),nrow=Ncat,ncol=Ndim);
for ( i in 1:Ndim){
 m[,i] <- runif(Ncat,min=min(data[,i]),max=max(data[,i]));
}
am <- am0;
aiv <- aiv0;
biv <- biv0;
api <- api0;
############################################################################
############# STARTIN ITERATIONS ###########################################
print("About to enter iterations");
t <- 0; inc <-1 ;
Cost.old <- rep(10^6,times=nCVconv);
while( (t < MaxIt) && (inc==1)){
 print(c("Iteration",t)); 

 # Update Categorical weights for updating posterior parameters
  print("Calling UpdateCatw");
  Catw <- UpdateCatw(Ncat,data,m,am,aiv,biv,api);
 # Update posterior parameters 
  print("Calling UpdateMix");
  NewVals <- UpdateMix(Ncat,data,m0,am0,aiv0,biv0,api0,m,am,aiv,biv,Catw$cwm,Catw$scw);

  if( verbatim==TRUE){
  print("Means");
  print(NewVals$mean);
  print("Weights");
  print(NewVals$dapi);
  }
  m <- NewVals$mean;
  am <- NewVals$ivarm;
  aiv <- NewVals$ivara;
  biv <- NewVals$ivarb;
  api <-NewVals$dapi;

  Cost <- CostKL(Ncat,data,m0,am0,aiv0,biv0,api0,m,am,aiv,biv,api,Catw$cwm);
  Cost.new <- Cost$ckl ;

  if( verbatim==TRUE){
  print("Cost");
  print(Cost$ckl);
  }
  
  if ( mean(abs(Cost.new-Cost.old)) < conv.tol ){
    inc <- 0;
  }
  Cost.old <- sort( Cost.old );
  Cost.old[nCVconv] <- Cost.new ; 
  t <- t + 1;
} # iteration loop
######### END OF ITERATIONS ######################################
  v.Costs <- c(v.Costs, Cost$ckl);
  v.NC <- c(v.NC, inc);
  l.NewVals[[run]] <- NewVals;
} # matches runs loop  
######### END OF RUNS #############################################

# Finished runs, pick best runs ( minimise the Cost function ) and compute cluster membership probabilities and assign samples to clusters.
sv.Costs <- sort( v.Costs, decreasing=FALSE, index.return=TRUE);
EstVals <- list(NULL);
MembPr2 <- list(NULL);
for ( i in 1:npick){
  EstVals[[i]] <- l.NewVals[[sv.Costs$ix[i]]];
  # Compute membership probs. for these  
  MembPr2[[i]] <- MembProbFn2(data,EstVals[[i]],Nsamples);
  wCl[i,] <- MembPr2[[i]]$wcl ;
}


 return(list(estvals=EstVals,wcl=wCl,probs=MembPr2, costs=v.Costs[sv.Costs$ix], conv=v.NC[sv.Costs$ix]));


} ### END OF FUNCTION
