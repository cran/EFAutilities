get.covR <-
function(GammaRho, p)
{
  c0 = seq(1,p^2,p+1)
  c1 = seq(1:p^2)[-c0]
  cov.list = cbind(expand.grid(1:ncol(GammaRho), 1:ncol(GammaRho)), c(GammaRho))
  cloc = diag(0, p, p); cloc[lower.tri(cloc)] = 1:ncol(GammaRho); cloc = cloc + t(cloc)
  ccloc = c(cloc)
  GammaR = matrix(NA, p^2, p^2)
  GammaR[c0,] = GammaR[,c0] = 0
  for (j in c1){
    for (k in c1){
      GammaR[j,k] = cov.list[which(cov.list[,1]==ccloc[j]&cov.list[,2]==ccloc[k]),3]
    }
  }
  return(GammaR)
}

