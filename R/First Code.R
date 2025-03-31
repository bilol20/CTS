######################################################
#####################  OUR TEST  #####################
######################################################
comp = function(x){
  ind = which(x[1]<x[-1])
  x1 = x[-1]
  y = rep(x[1], length.out = length(x)-1)
  y[ind] = x1[ind]
  return(y)
}

Stat = function(X,Y,Z){
  if(class(X)[1] == "numeric"){
    n = length(X)
  }else{
    n = nrow(X)
  }
  X1 = cbind(X,Z)
  Y1 = cbind(Y,Z)

  V = rbind(X1,Y1)
  m = 2*n
  D = as.matrix(dist(V))
  s = median(D)
  F = pexp


  T = sapply(1:n, function(i){
    M1 = cbind(D[1:m,i], D[1:m, 1:n])/s
    M2 = cbind(D[1:m,i], D[1:m, 1:n+n])/s
    M3 = cbind(D[1:m,i+n], D[1:m, 1:n+n])/s

    2*sum(F(apply(M2, 1, comp)))-sum(F(apply(M1, 1, comp)))-sum(F(apply(M3, 1, comp)))

  })
  return(sum(T)/(n^3))
}

Test.CTS = function(X,Y,Z, R){
  if(class(X)[1] == "numeric"){
    n = length(X)
  }else{
    n = nrow(X)
  }
  X1 = cbind(X,Z)
  Y1 = cbind(Y,Z)

  V = rbind(X1,Y1)
  m = 2*n
  D = as.matrix(dist(V))
  S = Stat_cpp(D)

  v = rbinom(R*n,1,0.5)
  S1 = lapply(1:R,function(l){
    index = which(v[1:n + l - 1]==1)
    pi = 1:(2*n)
    pi[c(index,n+index)] = pi[c(n+index,index)]
    D1 = D[pi,pi]
    return(Stat_cpp(D1))
  })

  pval = (1+sum(S1>S)+sum(S1==S))/(1+R)
  return(pval)
}
