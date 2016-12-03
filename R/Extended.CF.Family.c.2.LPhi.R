Extended.CF.Family.c.2.LPhi <-
function(Lambda,Phi,rotation=NULL,normalize=FALSE) {

## ### Internal functions: The.M.Matrix, Derivative.Constraints.Loadings, Derivative.Constraints.Phi

#--------------------------------------------------------------------------------------

### The function The.M.Matrix compute the M matrix described
### in Jennrich's 1973 Psychometrika paper on oblique rotation

 The.M.Matrix <- function (k1,k2,k3,k4,Lambda) {

# The function invokes no other external functions.

 
 Lambda2 = Lambda * Lambda
 
 sum.total = sum(Lambda2)
 sum.row = apply(Lambda2,1,sum) # p by 1 vector
 sum.column = apply(Lambda2,2, sum) # m by 1 vector
 
 p = dim(Lambda)[1]
 m = dim(Lambda)[2]
 
 Result = Lambda # Give the dimension of Result
 
 for (j in 1:m) {
   for (i in 1:p) {
   Result[i,j] = k1 * sum.total + k2 * sum.row[i] + k3 * sum.column[j] + k4 * Lambda2[i,j]
 } # i
 } # j
 
 Result
 
 } # The.M.Matrix <- function


#---------------------------------------------------------------------------------------

### The function Derivative.Constraints.Loadings implements Equation 53 
### of Jennrich, 1973, Psychometrika, 38, 593-604
### Guangjian Zhang, 2009-06-17, Wednesday

Derivative.Constraints.Loadings <- function (k1,k2,k3,k4,Lambda,Phi, M.Matrix ) {

# The function invokes no other external functions.

#library(MASS)

 p = dim(Lambda)[1]
 m = dim(Lambda)[2]


Phi.inverse = solve(Phi)

Term1 = (Lambda * M.Matrix) %*% Phi.inverse
Term3 = t(Lambda) %*% Lambda %*% Phi.inverse
Term4 = Lambda %*% Phi.inverse
Term5 = t(Lambda) %*% Lambda

Result.Loadings <- array (rep(0,m*m*p*m), dim=c(m,m,p,m))

for (v in 1:m) {
  for (u in 1:m) {
     for (r in 1:m) {
         for (i in 1:p) {
                   if (v != u) { 
                   delta.ur=0
                   if (u==r) delta.ur = 1
                   Result.Loadings[u,v,i,r] = delta.ur * Term1[i,v] +
                                           M.Matrix[i,r] * Lambda[i,u] * Phi.inverse[r,v] +
                                           2 * k1 * Lambda[i,r] * Term3[u,v] +
                                           2 * k2 * Lambda[i,r] * Lambda[i,u] * Term4[i,v] +
                                           2 * k3 * Lambda[i,r] * Term5[u,r] * Phi.inverse[r,v] +
                                           2 * k4 * Lambda[i,r] * Lambda[i,r] * Lambda[i,u] * Phi.inverse[r,v]  
                               } # if (v != u) 
                        } # (for i in 1:p)
                    } # (r in 1:m)
                 } # (u in 1:m)
               } # (v in 1:m)

Result.Loadings

} # Derivative.Constraints.Loadings

#---------------------------------------------------------------------------

### The function Derivative.Constraints.Phi implements Equation 54 
### of Jennrich, 1973, Psychometrika, 38, 593-604
### Guangjian Zhang, 2009-06-17, Wednesday

Derivative.Constraints.Phi <- function (k1,k2,k3,k4,Lambda,Phi, M.Matrix ) {

# The function invokes no other external functions.

#library(MASS)

 p = dim(Lambda)[1]
 m = dim(Lambda)[2]


Phi.inverse = solve(Phi)

Constraint = t(Lambda) %*% (Lambda * M.Matrix) %*% Phi.inverse

Result.Phi <- array (rep(0,m*m*m*m), dim=c(m,m,m,m))

for (v in 1:m) {
  for (u in 1:m) {
     for (y in 1:m) {
         for (x in 1:m) {
                   if ( v != u && x != y )  { 
                               delta.ux =0
                               delta.uy =0
                               if (u==x) delta.ux=1
                               if (u==y) delta.uy=1
                               Result.Phi[u,v,x,y] = - (delta.ux * Phi.inverse[y,v] + delta.uy * Phi.inverse[x,v]) * Constraint[u,u]                                
                               } # if ( v != u && x != y )
                        } # (for x in 1:m)
                    } # (y in 1:m)
                 } # (u in 1:m)
               } # (v in 1:m)

Result.Phi

} # Derivative.Constraints.Phi

#--------------------------------------------------------------------------------------------


p = dim(Lambda)[1]
m = dim(Lambda)[2]
Nc = m * (m-1) # the number of constraints
Nq = p * m + m * (m - 1) / 2 + p  # the number of parameters


if (normalize) {
Lambda0 = Lambda
Lambdao = tcrossprod (Lambda, Phi)
   h = rowSums(Lambda * Lambdao)

# h = rowSums(Lambda0 **2)
h.half.inverse =  1 / sqrt(h)
h.onehalf.inverse = (h.half.inverse)**3
 for (k in 1:p) {
  Lambda[k,] = Lambda[k,] * h.half.inverse[k]
  } # (k in 1:p)

} # (normalize)


if (is.null(rotation)) rotation='CF-varimax'

if (rotation=='CF-varimax') {

k1=0
k2=  (p-1)/p
k3=  1/p
k4= - 1


} else if (rotation=='CF-quartimax') {

k1=0
k2=  1
k3=  0
k4= - 1

} else {
  stop ("wrong specification for the factor rotation criterion")
}


M.Temp = The.M.Matrix (k1,k2,k3,k4,Lambda)

d.Con.Loading = Derivative.Constraints.Loadings (k1,k2,k3,k4, Lambda, Phi, M.Temp )
d.Con.Phi = Derivative.Constraints.Phi (k1,k2,k3,k4,Lambda,Phi, M.Temp )


d.Con.Parameters = array ( rep(0, Nc * Nq), dim = c (Nc, Nq) )

ICon = 0
for (j in 1:m) {
  for (i in 1:m)  {

     if (i != j) { # We need to worry about the off-diagonal elements of the Constraint matrix
       
      ICon = ICon + 1
      
      # Factor loadings
      tempc2L = d.Con.Loading[i,j,1:p,1:m]

if (normalize) {
      
      temp = rowSums(tempc2L * Lambda0)     
      for (k in 1:p) {
      tempc2L[k,] = tempc2L[k,] * h.half.inverse[k] - temp[k] * h.onehalf.inverse[k] * Lambdao[k,]
      } # k

} # (normalize)


      for (l in 1:m) {
         d.Con.Parameters[ICon,((l-1)*p+1):(l*p)] = tempc2L[1:p,l]            
      } # (l in 1:m)
 
      # Factor correlations
      Itemp = p * m      

      tempc2Phi = d.Con.Phi[i,j,1:m,1:m]

 if (normalize) {
 tempc2L = d.Con.Loading[i,j,1:p,1:m]
 M1 = tempc2L * Lambda0 
 for (ii in 1:p) {
  M1[ii,] = M1[ii,] * h.onehalf.inverse[ii]
  }


 } # (normalize)

      for (l in 2:m) {
         for (k in 1:(l-1)) {
           Itemp = Itemp + 1

if (normalize) {
  Mt = matrix(0,p,m)
  for (ii in 1:p) {  
  Mt[ii,] = M1[ii,] * Lambda0[ii,l] * Lambda[ii,k]
  }
  tempc2Phi[k,l] = - sum(Mt) + tempc2Phi[k,l]  
}  ####  normalize 


           d.Con.Parameters[ICon,Itemp] =  tempc2Phi[k,l] 

        } # (k in 1:l)
      } # (l in 2:m)

   } # (i != j)

  } # (i in 1:m)
} # (j in 1:m)



d.Con.Parameters

} # Extended.CF.Family.c.2.LPhi

