"unbiasedKurt" <-
function(v){
  z <- (v-mean(v,na.rm=T))/sqrt(var(v,na.rm=T));
  n <- length(which(is.na(v)==F));
  k <- n*(n+1)*sum(z^4)/((n-1)*(n-2)*(n-3)) - (3*(n-1)^2)/((n-2)*(n-3));
  return(k);
}

