## 2022-04-14, Thursday
## The task is to add the option of correlated residuals.


## 2018-05-21, Monday
## two updates: taking care of negative eigenvalues for the reduced correlation matrix
## and add multiple starting values for ML estimation


## 2016-06-02, Thursday, Guangjian Zhang
## the function fa.extract extracts factors from a correlation matrix
## using ols

### Guangjian Zhang, Thursday, 2016-08-11.
### The function fa.extract is modified to allow both ols and ml.



fa.extract <- function(covmat,factors,extraction=NULL,start=NULL,I.cr=NULL) {
  
## fa.extract contains seven internal functions, fnOLS, fgOLS, fnOLS.cr, fgOLS.cr, fnML.cr, fgML.cr, psi2A   
  

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
  
  #......................................
  # The beginning of fnML.cr
  fnML.cr <- function(psi.cr,S,m,I.cr) {
    
    p=ncol(S)
    i.cr = nrow(I.cr)  
    
    for (k in 1:i.cr) { # 2022-06-01, GZ, taking care of negative eigenvalues of Psi
      if ( (psi.cr[I.cr[k,1]] * psi.cr[I.cr[k,2]]) < psi.cr[p+k]^2) psi.cr[p+k] = sqrt(psi.cr[I.cr[k,1]] * psi.cr[I.cr[k,2]]*0.98)  * sign(psi.cr[p+k])    
    }
    
    
    Psi.cr = diag(psi.cr[1:p])
    
    for (i in 1:i.cr) {
      Psi.cr[I.cr[i,1],I.cr[i,2]] = psi.cr[p+i]
      Psi.cr[I.cr[i,2],I.cr[i,1]] = psi.cr[p+i]
    }      
    
    
    E = eigen(Psi.cr, symmetric = TRUE)
    
    Psi.half.inverse = E$vectors %*% diag(1/sqrt(E$values)) %*% t(E$vectors)    # values could be negative? GZ, 2022-04-12
    
    S.2 = Psi.half.inverse %*% S %*% Psi.half.inverse 
    
    E = eigen(S.2, symmetric = TRUE)
    
    return( sum(E$values[(m+1):p]) - sum(log(E$values)[(m+1):p]) - (p-m))
    
  } # fnML.cr
  
  # The end of fnML.cr
  #..............................................
  
  ###......................................
  ### The beginning of fgML.cr
  
  fgML.cr <- function(psi.cr,S,m,I.cr) {
    
    p=ncol(S)
    i.cr = nrow(I.cr)  
    
    for (k in 1:i.cr) { # 2022-06-01, GZ, taking care of negative eigenvalues of Psi
      if ( (psi.cr[I.cr[k,1]] * psi.cr[I.cr[k,2]]) < psi.cr[p+k]^2) psi.cr[p+k] = sqrt(psi.cr[I.cr[k,1]] * psi.cr[I.cr[k,2]]*0.98)  * sign(psi.cr[p+k])    
    }
    
    
    
    Psi.cr = diag(psi.cr[1:p])
    
    for (i in 1:i.cr) {
      Psi.cr[I.cr[i,1],I.cr[i,2]] = psi.cr[p+i]
      Psi.cr[I.cr[i,2],I.cr[i,1]] = psi.cr[p+i]
    }      
    
    
    E = eigen(Psi.cr, symmetric = TRUE)
    
    Psi.half.inverse = E$vectors %*% diag(1/sqrt(E$values)) %*% t(E$vectors)
    Psi.half = E$vectors %*% diag(sqrt(E$values)) %*% t(E$vectors)
    
    
    S.2 = Psi.half.inverse %*% S %*% Psi.half.inverse 
    
    E = eigen(S.2, symmetric = TRUE)
    
    
    if (m>1) {
      M.temp = diag(sqrt(ifelse(E$values[1:m]>=1,E$values[1:m]-1,0)))
    } else {
      M.temp = matrix(sqrt(ifelse(E$values[1]>=1,E$values[1]-1,0)),1,1)  
    }
    
    A = Psi.half %*% E$vectors[,1:m] %*% M.temp
    
    P.hat = A %*% t(A) + Psi.cr
    
    M.temp.2 = solve(P.hat) - solve(P.hat) %*% S %*% solve(P.hat)      
    
    result = psi.cr
    
    result[1:p] =  diag(M.temp.2)
    
    for (i in 1:i.cr) {
      result[p+i] = 2 * M.temp.2[I.cr[i,1],I.cr[i,2]]
    }      
    
    return(result)
    
  } # fgML.cr 
  
  ### The end of fgML.cr
  ###....................................................................
  
  
  #......................................
  fnOLS.cr <- function(psi.cr,S,m,I.cr) {
    
    p=ncol(S)
    i.cr = nrow(I.cr)  
    
    Psi.cr = diag(psi.cr[1:p])
    
    for (i in 1:i.cr) {
      Psi.cr[I.cr[i,1],I.cr[i,2]] = psi.cr[p+i]
      Psi.cr[I.cr[i,2],I.cr[i,1]] = psi.cr[p+i]
    }      
    
    R.reduced = S - Psi.cr
    
    E = eigen(R.reduced, symmetric = TRUE)
    
    if (m>1) {
      M.temp = diag(sqrt(ifelse(E$values[1:m]>=0,E$values[1:m],0)))
    } else {
      M.temp = matrix(sqrt(ifelse(E$values[1]>=0,E$values[1],0)),1,1)  
    }
    
    
    A = E$vectors[,1:m] %*% M.temp
    
    Residual = S - A %*% t(A) - Psi.cr
    
    #     Residual = Residual - diag(diag(Residual))
    
    return( sum(Residual ** 2) /2)
    
  } # fnOLS.cr
  
  #--------------------------------------------------------
  
  #......................................
  fgOLS.cr <- function(psi.cr,S,m,I.cr) {
    
    p=ncol(S)
    i.cr = nrow(I.cr)  
    
    Psi.cr = diag(psi.cr[1:p])
    
    for (i in 1:i.cr) {
      Psi.cr[I.cr[i,1],I.cr[i,2]] = psi.cr[p+i]
      Psi.cr[I.cr[i,2],I.cr[i,1]] = psi.cr[p+i]
    }      
    
    R.reduced = S - Psi.cr
    
    E = eigen(R.reduced, symmetric = TRUE)
    
    if (m>1) {
      M.temp = diag(sqrt(ifelse(E$values[1:m]>=0,E$values[1:m],0)))
    } else {
      M.temp = matrix(sqrt(ifelse(E$values[1]>=0,E$values[1],0)),1,1)  
    }
    
    A = E$vectors[,1:m] %*% M.temp
    
    Residual = S - A %*% t(A) - Psi.cr
    
    #   - diag(Residual) # Remove a bug on 2022-04-12, GZ
    
    result = psi.cr
    result[1:p] = - diag(Residual)
    
    for (i in 1:i.cr) {
      result[p+i] = -2 * Residual[I.cr[i,1],I.cr[i,2]]
    }      
    
    return(result)
    
  } # fgOLS.cr 
  
  
  ###------------------------------------------------------------
  
  ##----------------------------------------------------------------------------------
  ## The function psi2A compute the unrotated factor loading matrix and the model implied
  ## correlation matrix according to the EFA.cr model.
  
  psi2A <- function (S, m, I.cr=NULL, psi.cr=NULL,extraction='ml') {
    
    p=ncol(S)
    i.boundary.cr = NA
    
    
    if (is.null(I.cr)) 
    {i.cr = 0} else{
      i.cr = nrow(I.cr)
    }
    
    
    Psi.cr = diag(psi.cr[1:p])


if (i.cr>0) {            
  
  i.boundary.cr = rep(0, i.cr)
  
  for (k in 1:i.cr) { # 2022-06-01, GZ, taking care of negative eigenvalues of Psi
    if ( (psi.cr[I.cr[k,1]] * psi.cr[I.cr[k,2]]) < psi.cr[p+k]^2) { 
    psi.cr[p+k] = sqrt(psi.cr[I.cr[k,1]] * psi.cr[I.cr[k,2]]*0.98)  * sign(psi.cr[p+k])    
    i.boundary.cr[k] = 1 
      }
  }
  
      for (i in 1:i.cr) {
      Psi.cr[I.cr[i,1],I.cr[i,2]] = psi.cr[p+i]
      Psi.cr[I.cr[i,2],I.cr[i,1]] = psi.cr[p+i]
    }      
 
    
  } # if (i.cr>0)
  
  
       
    if (extraction=='ml') {
      
      E = eigen(Psi.cr, symmetric = TRUE)
      
      Psi.half.inverse = E$vectors %*% diag(1/sqrt(E$values)) %*% t(E$vectors)
      Psi.half = E$vectors %*% diag(sqrt(E$values)) %*% t(E$vectors)
      
      
      S.2 = Psi.half.inverse %*% S %*% Psi.half.inverse 
      
      E = eigen(S.2, symmetric = TRUE)
      
      
      if (m>1) {
        M.temp = diag(sqrt(ifelse(E$values[1:m]>=1,E$values[1:m]-1,0)))
      } else {
        M.temp = matrix(sqrt(ifelse(E$values[1]>=1,E$values[1]-1,0)),1,1)  
      }
      
      A = Psi.half %*% E$vectors[,1:m] %*% M.temp
      
      
    } else if (extraction=='ols') {
      
      R.reduced = S - Psi.cr
      
      E = eigen(R.reduced, symmetric = TRUE)
      
      if (m>1) {
        M.temp = diag(sqrt(ifelse(E$values[1:m]>=0,E$values[1:m],0)))
      } else {
        M.temp = matrix(sqrt(ifelse(E$values[1]>=0,E$values[1],0)),1,1)  
      }
      
      
      A = E$vectors[,1:m] %*% M.temp
      
      
    } else {
      stop ("wrong specification for the factor extraction method")
    }
    
    
    P.hat = A %*% t(A) + Psi.cr
    
    return(list(A=A, P.hat=P.hat,i.boundary.cr=i.boundary.cr,psi.cr=psi.cr))
    
  }
  # The end of the function psi2A
  #---------------------------------------------------------------------
  #---------------------------------------------------------------------
  
  # The start of the main part of fa.extract
  
  if (is.null(extraction)) extraction = 'ols'
  m=factors
  p=ncol(covmat)
  
  if (is.null(I.cr)) 
  {n.cr = 0} else{
    n.cr = nrow(I.cr)
  }
  
  
  if (is.null(start)) { 
    start <- rep(0,p + n.cr)
    start[1:p] <- 1/diag(solve(covmat))}
  
  
  if (n.cr == 0) { # regular EFA 
    
    if (extraction == 'ols')  {
      tr = optim(start, fn=fnOLS, gr=fgOLS, method = c("L-BFGS-B"), lower = 0.001, 
                 upper = 1, S = covmat, m = factors)
      
      
      A.lst = psi2A (S=covmat, m=factors, psi.cr = tr$par, extraction='ols')
      
    } # (extraction == 'ols')
    
    
    
    
    if (extraction=='ml') { # We use the R intrinsic function factanal for EFA with ML
      
      #    FA.ml = factanal(covmat = covmat, factors=factors, rotation="none") # 2018-05-11
      FA.ml = factanal(covmat = covmat, factors=factors, rotation="none",nstart=4,lower=0.001)
      Unrotated = matrix(FA.ml$loadings[1:p,1:m],nrow=p,ncol=m)
      f = FA.ml$criteria[1]
      
      if (FA.ml$converged) {
        convergence = 0
      } else {
        convergence = 1
      }
      
      heywood = length(which(FA.ml$uniqueness < 0.00101))
      
    } # ((extraction=='ml')
    
    
    
  } else if (n.cr > 0) { # EFA with correlated residuals
    
    lower.v = rep(0.0001, p + n.cr)
    upper.v = rep(1, p + n.cr)
    lower.v[(p+1):(p+n.cr) ] = -0.99
    upper.v[(p+1):(p+n.cr) ] = 0.99
    
    tr = optim(start, fn=fnOLS.cr, gr=fgOLS.cr, method = c("L-BFGS-B"), lower = lower.v, 
               upper = upper.v, S = covmat, m = factors, I.cr = I.cr)
    
    
    if (extraction=='ml') {
      start = tr$par
      tr = optim(start, fn=fnML.cr, gr=fgML.cr, method = c("L-BFGS-B"), lower = lower.v, 
                 upper = upper.v, S = covmat, m = factors, I.cr = I.cr)
      
    } # (extraction=='ml')
    
    A.lst = psi2A (S=covmat, m=factors, I.cr = I.cr, psi.cr = tr$par, extraction = extraction)
    
    
  }
  
  
  

  # 2022-04-14, 1:14pm
  # 
    
  if (! ( (extraction=='ml') & (n.cr == 0) )) { # We do not need to compute unrotated factor loadings
                                                # and other summary statistics for ML with EFA
 
  if (is.na(A.lst$i.boundary.cr[1])) {     

        Unrotated = A.lst$A
        f = tr$value
        convergence = tr$convergence
        heywood = length(which(tr$par[1:p] < 0.00101))
        psi.cr = tr$par

      } else {

        Unrotated = A.lst$A
        convergence = tr$convergence
        heywood = length(which(tr$par[1:p] < 0.00101))
        psi.cr = A.lst$psi.cr
        
        if (extraction=='ols') f = fnOLS.cr(psi.cr=psi.cr,S=covmat,m=factors,I.cr=I.cr)
        if (extraction=='ml') f = fnML.cr(psi.cr=psi.cr,S=covmat,m=factors,I.cr=I.cr)
            
  }
    
    }
  
  
  
  
  
  compsi = matrix(0,p,4)
  compsi[1:p,1] = eigen(covmat)$values
  compsi[1:p,2] = 1 - start[1:p] 
  compsi[1:p,3] = diag(Unrotated %*% t(Unrotated))
  compsi[1:p,4] = 1 - compsi[1:p,3] 
  
  colnames(compsi) = c('Eval','SMC','Comm','UniV')
  names(f) = 'Discrepancy'
  
  if (is.null(I.cr)) { 
  result = list(Unrotated=Unrotated, f = f, convergence = convergence, heywood = heywood, compsi=compsi,psi.cr=compsi[1:p,4]) # Add psi.cr to be compatible with bootjack
  } else {
    result = list(Unrotated=Unrotated, f = f, convergence = convergence, heywood = heywood, compsi=compsi,psi.cr=psi.cr, i.boudary.cr=A.lst$i.boundary.cr)
  }
  
  
    
  return(result)

  } # fa.extract
