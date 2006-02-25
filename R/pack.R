"pack" <-
function(data.m,kurt.range=c(-2,0),cluster=T,method=c("bic","vb")){

# compute kurtosis  
k.v <- vector();
for(r in 1:nrow(data.m)){
 k.v[r] <- unbiasedKurt(data.m[r,]);
}
feat.idx <- intersect(which(k.v >= kurt.range[1]),which(k.v <= kurt.range[2]));
print(paste("There are ",length(feat.idx)," genes in the specified kurtosis range",sep=""));

class.lv <- list();
sizeSubG.v <- vector();
if(cluster){
  # check for missing entries
  if( length(which(is.na(as.vector(data.m)))) > 0){
   print("Matrix has missing values. Please impute them and rerun.");
   break;
  }
  if(method=="bic"){
   library(mclust);
   sel.idx <- vector();
   i <- 1;
   for(r in feat.idx ){
    m <- matrix(data.m[r,],ncol=1,nrow=ncol(data.m));
    em.o <- Mclust(m,minG=1,maxG=2);
    ncl <- length(levels(factor(em.o$classification)));
    sizeCluster.v <- summary(factor(em.o$classification));
    if (ncl > 1){
     sel.idx <- c(sel.idx,r);
     sizeSubG.v[r] <- min(sizeCluster.v);
     class.lv[[r]] <- em.o$classification;
    }
    print(paste("Done ",i," out of ",length(feat.idx)," genes",sep=""));
    i <- i+1;
   }
  }
  else if (method=="vb"){
   sel.idx <- vector(); i <-1;
   for(r in feat.idx){
    m <- matrix(data.m[r,],ncol=1,nrow=ncol(data.m));
    prior <- UseBasicPrior(m,weights.v=c(1,1));
    vb.o <- vabayelMix(m,prior=prior,Ncat=2,nruns=5,npick=1,verbatim=F,MaxIt=100);
    ncl <- length(levels(factor(vb.o$wcl)));
    sizeCluster.v <- summary(factor(vb.o$wcl));
    if (ncl > 1){
     sel.idx <- c(sel.idx,r);
     sizeSubG.v[r] <- min(sizeCluster.v);
     class.lv[[r]] <- vb.o$wcl;
    }
    print(paste("Done ",i," out of ",length(feat.idx)," genes",sep=""));
    i <-i+1;
   }
 }
 out.m <- matrix(nrow=length(sel.idx),ncol=3);
 colnames(out.m) <- c("Kurtosis","Subgroup size","Index");
 k.s <- sort(k.v[sel.idx],decreasing=F,index.return=T);
 out.m[,1] <- k.s$x;
 out.m[,2] <- sizeSubG.v[sel.idx[k.s$ix]];
 out.m[,3] <- sel.idx[k.s$ix];
 rownames(out.m) <- rownames(data.m)[sel.idx[k.s$ix]];
}
else { # don't cluster
 sel.idx <- feat.idx;
 out.m <- matrix(nrow=length(sel.idx),ncol=2);
 colnames(out.m) <- c("Kurtosis","Index");
 k.s <- sort(k.v[sel.idx],decreasing=F,index.return=T);
 out.m[,1] <- k.s$x;
 out.m[,2] <- sel.idx[k.s$ix];
 rownames(out.m) <- rownames(data.m)[sel.idx[k.s$ix]];
  class.lv <- NULL;
}

 return(list(out=out.m,class=class.lv));
 
} # end of function

