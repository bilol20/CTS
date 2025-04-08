library(doParallel)
cl = makeCluster(detectCores())
registerDoParallel(cl)

#Example 1
n = 100
del = seq(0,0.5, length.out = 5)
dx = 10
dy = 10
dz = 10
var = diag(0.5,dx+dy+dz) + matrix(0.5,dx+dy+dz,dx+dy+dz)

dataset_ex1 = foreach(delta = del, .packages = c("mnormt"))%dopar%{
  mu = c(rep(0,dx), rep(delta,dy), rep(0,dz))
  U = rmnorm(n, mean = mu, varcov = var)

  X = U[,1:dx]
  Y = U[,dx+ 1:dy]
  Z = U[,dx+dy+1:dz]

  X1 = cbind(X,Z)
  Y1 = cbind(Y,Z)
  D = as.matrix(dist(rbind(X1,Y1)))


  return(list(dist_mat = D, X = X, Y = Y, Z = Z))
}


#Example 2
n = 100
del = seq(1,1.5, length.out = 5)
dx = 10
dy = 10
dz = 10

dataset_ex2 = foreach(delta = del, .packages = c("mnormt"))%dopar%{
  varz = diag(0.5,dz) + matrix(0.5,dz,dz)
  Z = rmnorm(n, varcov = varz)
  X = matrix(0,ncol = dx, nrow = n)
  Y = matrix(0,ncol = dy, nrow = n)
  for(i in 1:n){
    X[i,] = rmnorm(1, mean = Z[i,], varcov = varz)
    Y[i,] = rmnorm(1, mean = delta*Z[i,], varcov = delta*delta*varz)
  }

  X1 = cbind(X,Z)
  Y1 = cbind(Y,Z)
  D = as.matrix(dist(rbind(X1,Y1)))


  return(list(dist_mat = D, X = X, Y = Y, Z = Z))
}



#Example 3
n = 100
del = seq(1,10, length.out = 5)

dataset_ex3 = foreach(delta = del, .packages = c("mnormt", "transport"))%dopar%{
  Z = rnorm(n)
  X = matrix(0, ncol = 500, nrow = n)
  Y = matrix(0, ncol = 500, nrow = n)
  for (i in 1:n) {
    X[i,] = rnorm(500, Z[i], sd = 1)
    Y[i,] = rnorm(500, Z[i]+delta, sd = 1)
  }

  X1 = cbind(X,Z)
  Y1 = cbind(Y,Z)
  U = rbind(X1,Y1)

  k = 2*n
  D = matrix(0,k,k)
  for(i in 1:k){
    for(j in 1:k){
      D[i,j] = wasserstein1d(U[i,1:500], U[j,1:500], p = 1) + abs(U[i,501]-U[j,501])
    }
  }


  return(list(dist_mat = D))
}


#Example 4
n = 100
del = seq(1,2, length.out = 5)

dataset_ex4 = foreach(delta = del, .packages = c("mnormt", "transport"))%dopar%{
  Z = rnorm(n)
  X = matrix(0, ncol = 500, nrow = n)
  Y = matrix(0, ncol = 500, nrow = n)
  for (i in 1:n) {
    X[i,] = rnorm(500, Z[i], sd = 1)
    Y[i,] = rnorm(500, Z[i], sd = delta)
  }

  X1 = cbind(X,Z)
  Y1 = cbind(Y,Z)
  U = rbind(X1,Y1)

  k = 2*n
  D = matrix(0,k,k)
  for(i in 1:k){
    for(j in 1:k){
      D[i,j] = wasserstein1d(U[i,1:500], U[j,1:500], p = 1) + abs(U[i,501]-U[j,501])
    }
  }


  return(list(dist_mat = D))
}



#Example 5
n = 100
del = seq(0,1, length.out = 5)

dataset_ex5 = foreach(delta = del, .packages = c("mnormt", "igraph", "Matrix"))%dopar%{
  Z = rnorm(n)
  p1 = exp(Z)/(1+exp(Z))
  p2 = exp(Z-delta)/(1+exp(Z-delta))
  U = list()
  for(i in 1:(2*n)){
    if(i<=n){
      U[[i]] = list(G = erdos.renyi.game(n = 200, p = p1, directed = FALSE, loops = FALSE),
                    Z = Z[i])
    }else{
      U[[i]] = list(G = erdos.renyi.game(n = 200, p = p2, directed = FALSE, loops = FALSE),
                    Z = Z[i-n])
    }
  }

  k = 2*n
  D = matrix(0,k,k)
  for(i in 1:k){
    for(j in 1:k){
      L1 = as.matrix(laplacian_matrix(U[[i]]$G))
      L2 = as.matrix(laplacian_matrix(U[[j]]$G))
      D[i,j] = norm(L1 - L2, type = "F") + abs(U[[i]]$Z-U[[j]]$Z)
    }
  }


  return(list(dist_mat = D))
}


#Example 6
n = 100
del = seq(0,1, length.out = 5)

dataset_ex6 = foreach(delta = del, .packages = c("mnormt", "igraph", "Matrix"))%dopar%{
  Z = rnorm(n)
  p1 = pnorm(z)
  p2 = pnorm(z-delta)
  U = list()
  for(i in 1:(2*n)){
    if(i<=n){
      U[[i]] = list(G = erdos.renyi.game(n = 200, p = p1, directed = FALSE, loops = FALSE),
                    Z = Z[i])
    }else{
      U[[i]] = list(G = erdos.renyi.game(n = 200, p = p2, directed = FALSE, loops = FALSE),
                    Z = Z[i-n])
    }
  }

  k = 2*n
  D = matrix(0,k,k)
  for(i in 1:k){
    for(j in 1:k){
      L1 = as.matrix(laplacian_matrix(U[[i]]$G))
      L2 = as.matrix(laplacian_matrix(U[[j]]$G))
      D[i,j] = norm(L1 - L2, type = "F") + abs(U[[i]]$Z-U[[j]]$Z)
    }
  }


  return(list(dist_mat = D))
}



#Example 7
n = 100
del = seq(0,1, length.out = 5)

dataset_ex7 = foreach(delta = del, .packages = c("mnormt", "igraph", "Matrix"))%dopar%{
  Z = rnorm(n)
  U = list()
  for(i in 1:(2*n)){
    if(i<=n){
      U[[i]] = list(G = barabasi.game(n = 200, power = abs(Z[i]), directed = FALSE, zero.appeal = 0),
                    Z = Z[i])
    }else{
      U[[i]] = list(G = barabasi.game(n = 200, power = abs(Z[i-n]) + delta, directed = FALSE, zero.appeal = 0),
                    Z = Z[i-n])
    }
  }

  k = 2*n
  D = matrix(0,k,k)
  for(i in 1:k){
    for(j in 1:k){
      L1 = as.matrix(laplacian_matrix(U[[i]]$G))
      L2 = as.matrix(laplacian_matrix(U[[j]]$G))
      D[i,j] = norm(L1 - L2, type = "F") + abs(U[[i]]$Z-U[[j]]$Z)
    }
  }


  return(list(dist_mat = D))
}


stopCluster(cl)
