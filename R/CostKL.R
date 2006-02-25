"CostKL" <-
function(Ncat,data,m0,am0,aiv0,biv0,api0,m,am,aiv,biv,api,Catwm){
 Nsamples <- dim(data)[1];
 Ndim <- dim(data)[2];

 sum1 <<- 0.5*sum( am0*(1/am + m*m) );
 sum2 <<- 0.5*sum( log(am) ) ;
 v <- as.vector( Catwm[Catwm > 0] );
 sum3 <<- sum( v*log(v) );
 sum4 <<- sum( lgamma(biv0)-lgamma(biv) + biv*log(aiv) - biv0*log(aiv0)) ;
 sum5 <<- sum( lgamma(api0)-lgamma(api)  );
 sum6 <<- lgamma(sum(api))-lgamma(sum(api0)) ;
 sum7 <<- -(0.5 + log(mean(am0)))*Ncat*Ndim ;
 sum8 <<- 0.5*Nsamples*log(2*pi);
 
 Ckl <- sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7 +sum8 ;

 return(list(ckl=Ckl));


}

