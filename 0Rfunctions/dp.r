library(Rfast)
library(Rcpp)

haplotype.mosaic <- function( differences_from_aDNA, map_differences_from_aDNA, jumpcost=10, threshold=0.1 ) {

  chrs = sort(unique( map_differences_from_aDNA[,1]));
  optimal.path = list();
  for( chr in chrs ) {
    w = which( map_differences_from_aDNA[,1] == chr );
    bp = map_differences_from_aDNA[w,2];
    dist = differences_from_aDNA[,w];
    path = do.DP( dist, jumpcost=jumpcost, threshold=threshold );
    optimal.path[[chr]] = data.frame( chr=rep(chr,length(path)), bp=bp, haplotype=path)
  }
  return ( optimal.path )
}
      
    

do.DP <- function( dist, jumpcost, threshold ) {

  dist2 = matrix(0,nrow=nrow(dist)+1, ncol=ncol(dist));
  dist2[ 1:nrow(dist), ] = dist/2;
  dist2[nrow(dist2),] = threshold;
  dist2 = ifelse( is.na(dist2), threshold, dist2 )
  optimal.path = DP( dist2, jumpcost );
  return(optimal.path);
}
    
    

sourceCpp( code='
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]


NumericVector DP(NumericMatrix dist, double jumpcost) {

  int nstate = dist.nrow();
  int nsnp = dist.ncol();
  NumericMatrix score(nstate,nsnp);
  NumericMatrix path(nstate,nsnp);
  NumericVector optimalpath(nsnp);
  int i, j, k;

  // initialise path

  for(i=0;i<nstate;i++)  {
    for(j=0;j<nsnp;j++) {
      score(i,j) = 0;
    }
    score(i,0) = dist(i,0);
  }


  for(j=1;j<nsnp;j++) {
    for(i=0;i<nstate;i++) {
      double minpath = 1.0e10;
      for(k=0;k<nstate;k++) {
	double x = score(k,j-1) + dist(i,j) + (i!=k)*jumpcost;
	if ( x < minpath ) {
	  minpath = x;
	  score(i,j) = x;
	  path(i,j) = k;
	}
      }
    }
  }

  j=nsnp-1;
  double s = 1.0e10;
  k = -1;
  
  for(i=0;i<nstate;i++) {
    if ( score(i,j) < s ) {
      s = score(i,j);
      k = i;
    }
  }

  while( j >= 0 ) {
    optimalpath(j) = k+1;
    k = path(k,j);
    j--;
  }

  return(optimalpath);

}')
