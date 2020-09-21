### Pseudo-poisson informtion criterion (PPIC) for selecting the
### number of FPCs
PPIC<-function(K,f_locpoly,t,N,f_mu,G.eigen_v,xi,xgrid,delta){
  if(K==1){
    tmp     = f_mu + outer(G.eigen_v[,1],xi[1,])/delta # density functions
    # tmp2    = f_mu_d + outer(tmp2[,1],xi[1,])/delta # derivatives of density functions
  } else{
    tmp     = f_mu + {G.eigen_v[,1:K]/delta}%*%xi[1:K,] # density functions
    # tmp2    = f_mu_d + {tmp2[,1:K]/delta}%*%xi[1:K,] # derivatives of density functions
  }
  ## make adjustments to get valid density functions (nonnegative+integrate to 1)
  tmp     = apply(tmp,2,function(x){
    x[x<0] = 0     # non-negative
    x      = x/sum(x) # integrate to delta^2
    return(x)
  })
  tmp     = tmp/delta^2
  n       = length(N)
  tmp     = lapply(seq(1,n),function(i) approx(xgrid,tmp[,i],t[(sum(N[1:i])-N[i]+1):sum(N[1:i])])$y)
  f_locpoly = unlist(f_locpoly)
  if(sum(f_locpoly<=0)>0) f_locpoly[f_locpoly<1e-8] = 1e-8
  tmp       = unlist(tmp)
  tmp[tmp<1e-8] = 1e-8
  D_K       = 2*sum(f_locpoly*log(f_locpoly/tmp)-(f_locpoly-tmp))+2*K
  return(D_K)
}

### local linear regression
den_locpoly<-function(partition,t,N){
  n           = length(N)
  l           = length(partition)
  partition_x = {partition[-1]+partition[-l]}/2
  cumsumN     = c(0,cumsum(N))
  f_locpoly   = lapply(seq(1,n),function(i){
    a         = {cumsumN[i]+1}:cumsumN[i+1]
    ### the histogram density estimate for the ith subject
    f_H       = hist(t[a],breaks=partition,plot=FALSE)$counts/N[i]/diff(partition)
    ### leave one out cross validation, local linear regression
    GCV_loclinear = function(hh){
      D   = sapply(partition_x,function(x) partition_x-x)
      W   = dnorm(D/hh)
      S0  = apply(W,1,sum)
      S1  = apply(W*D,1,sum)
      S2  = apply(W*{D^2},1,sum)
      L   = W*{t(VTM(S2,l-1))-D*t(VTM(S1,l-1))}/t(VTM(S0*S2-S1^2,l-1))
      f_oneout = L%*%f_H
      return(sum({{f_H-f_oneout}/{1-mean(diag(L))}}^2))
    }
    ### find bandwidth that minimize generalized cross validation
    # here set lower limit to be the min(diff(partition_x)) to avoid error
    hh   = optimize(GCV_loclinear,lower=min(diff(partition_x)),upper=IQR(t)*2)$minimum
    ### use the obtained GCV bandwidth to get local linear regression estimates
    D   = t(sapply(t[a],function(x) partition_x-x))
    W   = dnorm(D/hh)
    S0  = apply(W,1,sum)
    S1  = apply(W*D,1,sum)
    S2  = apply(W*{D^2},1,sum)
    L   = W*{t(VTM(S2,l-1))-D*t(VTM(S1,l-1))}/t(VTM(S0*S2-S1^2,l-1))
    return(L%*%f_H)
  })
}

### local linear regression
den_locpoly2<-function(t,N){
  n           = length(N)
  nbreak      = numeric(n)
  cumsumN     = c(0,cumsum(N))
  f_locpoly   = lapply(seq(1,n),function(i){
    a         = {cumsumN[i]+1}:cumsumN[i+1]
    ### the histogram density estimate for the ith subject
    f_H       = hist(t[a],plot=FALSE)
    partition = f_H$breaks
    nbreak[i] = length(partition)
    f_H       = f_H$counts/N[i]/diff(partition)
    partition_x = {partition[-1]+partition[-length(partition)]}/2
    ### leave one out cross validation, local linear regression
    GCV_loclinear = function(hh){
      D   = sapply(partition_x,function(x) partition_x-x)
      W   = dnorm(D/hh)
      S0  = apply(W,1,sum)
      S1  = apply(W*D,1,sum)
      S2  = apply(W*{D^2},1,sum)
      L   = W*{t(VTM(S2,nbreak[i]-1))-D*t(VTM(S1,nbreak[i]-1))}/t(VTM(S0*S2-S1^2,nbreak[i]-1))
      f_oneout = L%*%f_H
      return(sum({{f_H-f_oneout}/{1-mean(diag(L))}}^2))
    }
    ### find bandwidth that minimize generalized cross validation
    # here set lower limit to be the min(diff(partition_x)) to avoid error
    hh   = optimize(GCV_loclinear,lower=min(diff(partition_x)),upper=IQR(t[a])*2)$minimum
    ### use the obtained GCV bandwidth to get local linear regression estimates
    D   = t(sapply(t[a],function(x) partition_x-x))
    W   = dnorm(D/hh)
    S0  = apply(W,1,sum)
    S1  = apply(W*D,1,sum)
    S2  = apply(W*{D^2},1,sum)
    L   = W*{t(VTM(S2,nbreak[i]-1))-D*t(VTM(S1,nbreak[i]-1))}/t(VTM(S0*S2-S1^2,nbreak[i]-1))
    return(L%*%f_H)
  })
}
