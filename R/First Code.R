######################################################
#####################  OUR TEST  #####################
######################################################

# Resampling based test

# @param X, Y, Z sampled data
# @param R number of iterations of the resampling algorithm

Test.CTS = function(X, Y=NULL, Z=NULL, distance = FALSE, R = 200){
  if(distance  = TRUE){
    S = Stat_cpp(X)

    v = rbinom(R*n,1,0.5)
    S1 = lapply(1:R,function(l){
      index = which(v[1:n + l - 1]==1)
      pi = 1:(2*n)
      pi[c(index,n+index)] = pi[c(n+index,index)]
      D1 = X[pi,pi]
      return(Stat_cpp(D1))
    })

    pval = (1+sum(S1>S)+sum(S1==S))/(1+R)

  }else{
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
  }
  return(pval)
}

BD.ctest = function(X, Y = NULL,Z = NULL, distance = FALSE, R = 200){
  if(distance = TRUE){
    S = bd(X, distance = TRUE, size = c(n,n))

    v = rbinom(R*n,1,0.5)
    S1 = lapply(1:R,function(l){
      index = which(v[1:n + l - 1]==1)
      pi = 1:(2*n)
      pi[c(index,n+index)] = pi[c(n+index,index)]
      D1 = X[pi,pi]
      return(bd(D1, distance = TRUE, size = c(n,n)))
    })

    pval = (1+sum(S1>S)+sum(S1==S))/(1+R)
  }else{
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
    S = bd(D, distance = TRUE, size = c(n,n))

    v = rbinom(R*n,1,0.5)
    S1 = lapply(1:R,function(l){
      index = which(v[1:n + l - 1]==1)
      pi = 1:(2*n)
      pi[c(index,n+index)] = pi[c(n+index,index)]
      D1 = D[pi,pi]
      return(bd(D1, distance = TRUE, size = c(n,n)))
    })

    pval = (1+sum(S1>S)+sum(S1==S))/(1+R)
  }

  return(pval)

}
