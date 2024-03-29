#### 2016-06-24, Friday, Guangjian Zhang. Non-normal distributions
#### 2016-06-02, Thursday, Guangjian Zhang


oblq.se.augmt <- function(Lambda, Phi, Rsample, N, extraction, rotation, normalize=FALSE, modelerror, geomin.delta=NULL, MTarget=NULL, MWeight=NULL,
                          PhiWeight = NULL, PhiTarget = NULL, u.r = NULL, acm.type, wxt2=1e0, I.cr=NULL, psi.cr=NULL) {

# It invokes external functions (EliU,D.g.2.r, EFA.Hessian, Extended.CF.Family.c.2.LPhi,Derivative.Constraints.Numerical) 



# We assume that the manifest variables are normally distributed for the time being.

if (is.null(rotation)) stop ("No rotaton criterion is specified for numberical approximation")
if ((rotation=='geomin') & (is.null(geomin.delta))) geomin.delta = 0.01
if ((rotation=='target') & ((is.null(MWeight)) | (is.null(MTarget)))) stop ("MWeight or MTarget is not specified for target rotation")
if ((rotation=='xtarget') & ((is.null(MWeight)) | (is.null(MTarget)) | (is.null(PhiWeight)) | (is.null(PhiTarget)) )) stop ("MWeight or MTarget is not specified for target rotation") # 2016-06-03, GZ
if (is.null(modelerror)) modelerror='YES' 

p = dim(Lambda)[1]
m = dim(Lambda)[2]

if (is.null(I.cr)) 
{n.cr = 0} else{
  n.cr = nrow(I.cr)
}


Nc = m * (m-1) # the number of constraints

Nq = p * m + m * (m - 1) / 2 + p + n.cr  # the number of parameters


PM = Lambda %*% Phi %*% t(Lambda) 
PM = PM - diag(diag(PM)) + diag(p)

if (n.cr>0) {
  for (i in 1:n.cr) {
    PM[I.cr[i,1],I.cr[i,2]] = PM[I.cr[i,1],I.cr[i,2]] + psi.cr[p+i]
    PM[I.cr[i,2],I.cr[i,1]] = PM[I.cr[i,2],I.cr[i,1]] + psi.cr[p+i]
  }
}


##

if (modelerror== 'NO') Rsample = PM

if (is.null(u.r)) u.r = EliU(Rsample)

# if (is.null(u.r)) {

# if (modelerror== 'YES') {
#  u.r = EliU(Rsample) 
#} else if (modelerror == 'NO') {
#  u.r = EliU(PM)
#} else {
# stop("Model Error option is inappropriately specified.")
#}
# } ## (u.r == NULL)

##

if (extraction == 'ml') {

dg2r = D.g.2.r (Lambda, Phi, extraction='ml',I.cr, psi.cr)
Hessian = EFA.Hessian (Lambda, Phi, Rsample, extraction='ml',I.cr, psi.cr) ####


} else if (extraction == 'ols') {

dg2r = D.g.2.r (Lambda, Phi, extraction='ols',I.cr, psi.cr)
Hessian = EFA.Hessian (Lambda, Phi, Rsample, extraction='ols',I.cr, psi.cr) #####


} else {

stop ('Factor extraction method is incorrectly specified.')
}


## start of 2020-05-13, GZ

if (acm.type==1) {
  Y.Hat = u.r
} else {
  
  
  Y.Hat = matrix(0, (p*(p-1)/2), (p*(p-1)/2))
  u.r.col = matrix(0, (p*p), (p*(p-1)/2))
  
  ij = 0
  for (j in 2:p) {
    for (i in 1:(j-1)) {
      ij = ij + 1
      u.r.col[,ij] = u.r[,((j-1)*p + i)]
    }
  }
  
  
  ij = 0
  for (j in 2:p) {
    for (i in 1:(j-1)) {
      ij = ij + 1
      Y.Hat[ij,] = u.r.col[((j-1)*p + i) , ]
    }
  }
  
}

dg2r.upper = matrix(0, (p*(p-1)/2), Nq)


ij = 0
for (j in 2:p) {
  for (i in 1:(j-1)) {
    ij = ij + 1
    dg2r.upper[ij,] = dg2r[((j-1)*p + i) , ]
  }
}

Ham = (t(dg2r.upper) %*% Y.Hat %*% dg2r.upper) * 4 # Multiplying by 4 to accommodate the partial derivatives

## end of 2020-05-13, GZ

# Ham = t(dg2r) %*% u.r %*% dg2r

###

if (rotation=='CF-varimax') {

Olq.Con.Parameters = Extended.CF.Family.c.2.LPhi (Lambda, Phi, 'CF-varimax',normalize)

} else if (rotation=='CF-quartimax') {

Olq.Con.Parameters = Extended.CF.Family.c.2.LPhi (Lambda, Phi, 'CF-quartimax',normalize)

} else if (rotation=='CF-facparsim') {
  
  Olq.Con.Parameters = Extended.CF.Family.c.2.LPhi (Lambda, Phi, 'CF-facparsim',normalize)
  
} else if (rotation=='CF-equamax') {
  
  Olq.Con.Parameters = Extended.CF.Family.c.2.LPhi (Lambda, Phi, 'CF-equamax',normalize)
  
} else if (rotation=='CF-parsimax') {
  
  Olq.Con.Parameters = Extended.CF.Family.c.2.LPhi (Lambda, Phi, 'CF-parsimax',normalize)
  
} else if (rotation=='geomin') {

if (is.null(geomin.delta)) geomin.delta = 0.01

Olq.Con.Parameters = Derivative.Constraints.Numerical(Lambda,Phi,'geomin',normalize,geomin.delta)

} else if (rotation=='target') {

Olq.Con.Parameters = Derivative.Constraints.Numerical(Lambda,Phi,'target',normalize,MWeight=MWeight, MTarget=MTarget)

} else if (rotation=='xtarget') {

Olq.Con.Parameters = Derivative.Constraints.Numerical(Lambda,Phi,'xtarget',normalize,MWeight=MWeight, MTarget=MTarget,PhiWeight = PhiWeight, PhiTarget = PhiTarget, wxt2 = wxt2) # 2017-11-28, GZ!


} else {

  stop ("wrong specification for the factor rotation criterion")
}

### 

nlphi = p * m + m *(m-1)/2 # the number of factor loadings and factor correlations
Temp.Bigger = matrix(0,(Nq+Nc),(Nq+Nc))
Temp.Bigger[1:Nq,1:Nq] = Hessian
Temp.Bigger[1:nlphi,(Nq+1):(Nq+Nc)] = t(Olq.Con.Parameters[1:Nc, 1:nlphi])
Temp.Bigger[(Nq+1):(Nq+Nc),1:nlphi] = Olq.Con.Parameters[1:Nc, 1:nlphi]


Temp.Bigger.inverse = solve(Temp.Bigger)

Sandwich = Temp.Bigger.inverse[1:Nq,1:Nq] %*% Ham %*% Temp.Bigger.inverse[1:Nq,1:Nq]

###

SE = sqrt(diag(Sandwich)/(N-1))

Lambda.se <- array(SE[1: (p*m) ],dim=c(p,m))

Phi.se <- matrix(0,m,m)
ij=0
for (j in 2:m) {
 for (i in 1:(j-1)) {
    ij = ij + 1
    Phi.se[i,j] = SE[p*m+ij]
  }
}

Phi.se = Phi.se + t(Phi.se)

Psi.se <- SE[(p*m+m*(m-1)/2+1):Nq]

##
## list(Lambda.se = Lambda.se, Phi.se = Phi.se, Psi.se = Psi.se, Temp.Bigger = Temp.Bigger, Hessian = Hessian)

list(Lambda.se = Lambda.se, Phi.se = Phi.se, Psi.se = Psi.se)

} # oblq.se.augmt

