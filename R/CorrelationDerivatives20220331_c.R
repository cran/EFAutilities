###  2016-06-02, Thursday, Guangjian Zhang
### Two external functions: EliU, D.g.2.r
### The function D.g.2.r was modified to allow correlated residuals, 2022-03-31, Guangjian Zhang. 


EliU <- function(MP, eta = 1) {

### Coded by Guangjian Zhang, 2016-05-20

## Browne, M. W. & Shapiro, A. (1986). The asymptotic Covariance matrix of 
## sample correlation coefficients under general conditions. Linear Algebra
## and its applications, 82, 169-176.
## Equations (4.1) and (4.3)


## EliU -> The asymptotic covariance matrix of sample correlations if manifest variables are
## of an elliptical distribution.

# It does not require any external functions.

p = dim(MP)[1]

Ms= matrix(0,p*p,p*p)

for (j in 1:p) {
 for (i in 1:p)  {
 
   if (j==i) {
   ii = (i-1)*p + i
   Ms[ii,ii] = 1
   } else
   {
   ij = (j-1)*p + i
   ji = (i-1)*p + j
   Ms[ij,ij] = 0.5
   Ms[ij,ji] = 0.5
   }

  } # i
} # j


Kd = matrix(0,p*p,p)
for (i in 1:p) {
 ii = (i-1) * p + i
 Kd[ii,i] = 1 
}


A = Ms %*% (MP %x% diag(p)) %*% Kd

Gamma = 2 * Ms %*% (MP %x% MP)

if (eta != 1 ) {
MP.v = array(MP)
Gamma = eta * Gamma + (eta -1) * outer(MP.v, MP.v)
}

B = Gamma %*% Kd
G = t(Kd) %*% Gamma %*% Kd

Cov.r = Gamma - A %*% t(B) - B %*% t(A) + A %*% G %*% t(A)

} # EliU


#################################################################
 
D.g.2.r <- function(Lambda, Phi, extraction=NULL,I.cr=NULL, psi.cr=NULL) {

   # D.g.2.r computes the partial derivatives of the gradient function WRT to correlations 
   # Lambda <- (p,m)
   # Phi <- (m,m), a symmetric matrix
   # I.cr <- (n.cr,2), a two-dimensional matrix of the locations of correlated residuals, 2022-03-30, GZ
   # psi.cr <- (p + n.cr), a vector of unique variances and correlated residuals   
   # Result (D.Gradient.2.R.vector) <- (p^2, q), a matrix of derivatives of the gradient WRT to variable correlations
   # or the partial derivatives of a discrepancy function WRT parameters and variable correlations
   

#-----------------------------------------------------------------------------------

   # DifS2LPhiPsi computes the partial derivatives of the manifest variable
   # covariance matrix WRT factor loadings, factor correlations, and unique variances.
   # Lambda <- (p,m)
   # Phi <- (m,m), a symmetric matrix
   # I.cr <- (i.cr,2), a two-dimensional matrix of the locations of correlated residuals, 2022-03-30, GZ
   # Result <- (p,p,(p*m + m*(m-1)/2 + p))
   
   DifS2LPhiPsi <- function(Lambda, Phi, I.cr=NULL){
      
      # The function invokes no other external functions.
      
      p = dim(Lambda)[1]
      m = dim(Lambda)[2]
      
      if (is.null(I.cr)) 
      {n.cr = 0} else{
         n.cr = nrow(I.cr)
      }
      
      
      
      Result = array( rep(0, p*p*( p*m + m*(m-1)/2 + p + n.cr)), dim=c(p,p,( p*m + m*(m-1)/2 + p + n.cr) ))
      
      ## Factor loadings
      
      ij = 0
      LPhi = Lambda %*% Phi
      
      for (j in 1:m) {
         for (i in 1:p) {
            ij = ij + 1
            Result[i,1:p,ij] = LPhi[1:p,j]
            Result[1:p,i,ij] = LPhi[1:p,j] 
         } # i
      } # j
      
      ## Factor correlations
      
      ij = p*m
      
      if (m>1) {
         for (j in 2:m) {
            for (i in 1:(j-1)) {
               ij = ij + 1    
               Temp = Lambda[1:p,i] %*% t(Lambda[1:p,j]) 
               Result[1:p,1:p,ij] = Temp + t(Temp)
            } # i
         } # j
      }
      
      ## unique variances
      ij = p*m + m*(m-1)/2
      
      for (i in 1:p) {
         ij = ij + 1
         Result[i,i,ij] = 1
      }
      
      
      ## Correlated residuals
      
      if (!(is.null(I.cr))) {
         
         for (i in 1:n.cr) {
            ij = ij + 1
            Result[I.cr[i,1],I.cr[i,2],ij] = 1
            Result[I.cr[i,2],I.cr[i,1],ij] = 1
         }
         
         
      }
      
      
      ## Output
      
      return(Result)
      
   } # DifS2LPhiPsi
   
   
#-------------------------------------------------------------------------------------


if (is.null(extraction)) extraction='ml'
p = dim(Lambda)[1]
m = dim(Lambda)[2]

if (is.null(I.cr)) 
{n.cr = 0} else{
   n.cr = nrow(I.cr)
}



Nq = p * m + m * (m - 1) / 2 + p + n.cr


PM = Lambda %*% Phi %*% t(Lambda)
for (i in 1:p) {
PM[i,i] = 1
}


if (n.cr>0) {
   for (i in 1:n.cr) {
      PM[I.cr[i,1],I.cr[i,2]] = PM[I.cr[i,1],I.cr[i,2]] + psi.cr[p+i] # correct a bug, 2022-04-05, GZ
      PM[I.cr[i,2],I.cr[i,1]] = PM[I.cr[i,2],I.cr[i,1]] + psi.cr[p+i]
   }
}


PInverse = solve(PM)

Gradient = DifS2LPhiPsi (Lambda, Phi, I.cr)

D.Gradient.2.R.vector = matrix(0,p*p,Nq)

if (extraction=='ml') {

for (i in 1:Nq) {
D.Gradient.2.R.vector[, i] = array((-PInverse %*% Gradient[1:p,1:p,i] %*% PInverse), dim=c(p*p,1)) 
}

} else if (extraction=='ols') {

for (i in 1:Nq) {
D.Gradient.2.R.vector[, i] = - 2 * array( Gradient[1:p,1:p,i], dim=c(p*p,1)) 
} 

} else {
  stop ("wrong specification for the factor extraction method")
}

return(D.Gradient.2.R.vector)

} # D.g.2.r

###############################################################################################