"UpdateMix" <-
function(Ncat,data,m0,am0,aiv0,biv0,api0,m,am,aiv,biv,lambda,s.lambda){

  # Ncat : Number of components to look for
  # Ndim : Number of dimensions of observation vector
  # m0   : initial hyperparameter for means, same dimension as m
  # m    : means matrix, of dimension Ncat x Ndim
  # am   : inverse variances of means, Ncat x Ndim. (am0)
  # aiv  : par for inverse variance posterior, Ncat x Ndim
  # biv  : par for inverse variance posterior, Ncat x Ndim
  Ndim <- dim(data)[2];
  #tmp var needed
  v <- rep(0,Ndim);
     
  # store old values
  m.o <- m ; am.o <- am ; biv.o <- biv ; aiv.o <- aiv ;

  for ( c in 1:Ncat){
   # Update of other posterior parameters
   # inverse variance of mean
   am[c,] <- am0[c,] + s.lambda[c]*biv.o[c,]/aiv.o[c,];
   # Update of mean
   for ( i in 1:Ndim){
     v[i] <- sum(lambda[c,]*data[,i]);
   }
   m[c,] <- ( m0[c,]+ v*biv.o[c,]/aiv.o[c,] )/am[c,] ;
     
   # par1 of inv. variance
   biv[c,] <- biv0[c,] + 0.5*s.lambda[c];


   for ( i in 1:Ndim){                                   # par2 of inv. variance
   aiv[c,i] <- aiv0[c,i] + 0.5*sum( lambda[c,]*( data[,i]*data[,i]-2*data[,i]*m[c,i]+ m[c,i]*m[c,i] + 1/am[c,i] ) ) ;

   }
  }

  # learn dirichlet posterior parameters
  api <- api0 + s.lambda ;

  return(list(mean = m, ivarm = am, ivarb = biv, ivara= aiv, dapi=api));

} # END of function

