"UpdateCatw" <-
function(Ncat,data,m,am,aiv,biv,api){
  Nsamples <- dim(data)[1];
  Ndim <- dim(data)[2];
  # Values needed
  s.api <- sum(api);
  t1 <- digamma(s.api);

  # temporary variables needed
  alpha <- matrix(nrow=Ncat, ncol=Nsamples);
  v <- rep(0,Ndim);
  lambda <- matrix(nrow=Ncat, ncol=Nsamples);
  s.lambda <- rep(0,Ncat);
  
  for ( c in 1:Ncat){
      
    value1 <- digamma(api[c]) - t1 + 0.5*( sum(-log(aiv[c,])+ digamma(biv[c,])) );
    vec1 <- m[c,]*m[c,] ;
    vec2 <- 1/am[c,];
    vec3 <- biv[c,]/aiv[c,] ;

   for ( n in 1:Nsamples){
#      print(c("Doing ",c,n));
      value2 <- sum( (data[n,]*data[n,]-2*data[n,]*m[c,]+ vec1 + vec2 )*vec3 );

      alpha[c,n] <- value1 - 0.5* value2;
   }
     
  }
     
  for ( c in 1:Ncat){
     s.lambda[c] <- 0;
   for( n in 1:Nsamples){

     lambda[c,n] <- 1/sum( exp( alpha[,n]-alpha[c,n] ) );
     s.lambda[c] <- s.lambda[c]+lambda[c,n];
   }
  }


   return(list(cwm=lambda, scw=s.lambda));


} # endoffunction

