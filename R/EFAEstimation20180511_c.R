
fa.extract <- function(covmat,factors,extraction=NULL,start=NULL) {
  
  #......................................
  fnOLS <- function(Psi,S,m) {
    
    R.reduced = S - diag(Psi)
    
    E = eigen(R.reduced, symmetric = TRUE)
 
    

    if (m>1) {
      M.temp = diag(sqrt(ifelse(E$values[1:m]>=0,E$values[1:m],0)))
    } else {
      M.temp = matrix(sqrt(ifelse(E$values[1]>=0,E$values[1],0)),1,1)  
    }
    
    
    A = E$vectors[,1:m] %*% M.temp
    
    Residual = S - A %*% t(A) - diag(Psi)

#     Residual = Residual - diag(diag(Residual))
    
    sum(Residual ** 2) /2
    
  } # 
  #......................................
  fgOLS <- function(Psi,S,m) {
    
    R.reduced = S - diag(Psi)
    
    E = eigen(R.reduced, symmetric = TRUE)
    
    if (m>1) {
      M.temp = diag(sqrt(ifelse(E$values[1:m]>=0,E$values[1:m],0)))
    } else {
      M.temp = matrix(sqrt(ifelse(E$values[1]>=0,E$values[1],0)),1,1)  
    }
    
    A = E$vectors[,1:m] %*% M.temp
    
    Residual = S - A %*% t(A) - diag(Psi)
    
    - diag(Residual)
    
  } # fgOLS 
  
  #.........................................
  
  if (is.null(extraction)) extraction = 'ols'
  if (is.null(start)) start <- 1/diag(solve(covmat))
  
  m=factors
  p=ncol(covmat)
  
  if (extraction == 'ols') {
    tr = optim(start, fn=fnOLS, gr=fgOLS, method = c("L-BFGS-B"), lower = 0.0001, 
               upper = 1, S = covmat, m = factors)
    
    R.reduced = covmat - diag(tr$par)
    
    E = eigen(R.reduced, symmetric = TRUE)
    
    if (m>1) {
      M.temp = diag(sqrt(E$values[1:m]))
    } else {
      M.temp = matrix(sqrt(E$values[1]),1,1)  
    }
    
    A = E$vectors[,1:m] %*% M.temp
    
    Residual = covmat - A %*% t(A)
    Residual = Residual - diag(diag(Residual))
    
    Unrotated=A
    f = sum(Residual**2) / 2
    convergence = tr$convergence
    heywood = length(which(tr$par < 0.00011))
    
  } else if (extraction=='ml') {
    
    FA.ml = factanal(covmat = covmat, factors=factors, rotation="none",nstart=4)
    Unrotated = matrix(FA.ml$loadings[1:p,1:m],nrow=p,ncol=m)
    f = FA.ml$criteria[1]
    
    if (FA.ml$converged) {
      convergence = 0
    } else {
      convergence = 1
    }
    
    heywood = length(which(FA.ml$uniqueness < 0.00011))
    
  } # else if (extract=='ml')
  
 
  list(Unrotated=Unrotated, f = f, convergence = convergence, heywood = heywood)

  } # fa.extract
