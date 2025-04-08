pval = function(X, Y, Z, R, alpha, t){
  p.CTS = CTS::Test.CTS(X,Y,Z, R)
  p.ECMMD = asy_test(list(X = X, Y = Y, Z = Z),
                              k = 10, kernel_choice = "Gaussian")$p_value
  p.bdc = CTS::BD.ctest(X,Y,Z, R)

  return(c(p.CTS, p.ECMMD, p.bdc))
}

library(doParallel)
cl = makeCluster(detectCores())
registerDoParallel(cl)



stopCluster(cl)
