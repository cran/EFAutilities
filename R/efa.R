#'@importFrom mvtnorm pmvnorm

#'@importFrom mvtnorm dmvnorm

#'@importFrom GPArotation cfQ

#'@importFrom GPArotation geominQ

#'@importFrom GPArotation pstQ

#'@importFrom GPArotation cfT

#'@importFrom GPArotation geominT

#'@importFrom GPArotation pstT

#'@importFrom stats acf

#'@importFrom stats cor

#'@importFrom stats dnorm

#'@importFrom stats factanal

#'@importFrom stats optim

#'@importFrom stats pnorm

#'@importFrom stats qnorm

#'@importFrom stats quantile

#'@importFrom stats sd

#'@importFrom stats var

#'@importFrom utils combn

#'@importFrom utils head

#'@export efa


efa <- function(x=NULL, factors=NULL, covmat=NULL, n.obs=NULL, dist='normal', fm='ols', rtype='oblique', rotation='CF-varimax', normalize=FALSE, geomin.delta=NULL, MTarget=NULL, MWeight=NULL,
                          PhiWeight = NULL, PhiTarget = NULL, useorder=FALSE, se='information', Ib=2000, mnames=NULL, fnames=NULL, merror='YES') {

## Internal functions: Make.Rot.Args

# library(GPArotation)
# library(polycor)
# source('E:/CurrentSimulation/RCPhi/EFAEstimation20160811.R')
# source('E:/CurrentSimulation/RCPhi/oblqSE20160622.R')
# source('E:/CurrentSimulation/RCPhi/orthSE20160624.R')
# source('E:/CurrentSimulation/RCPhi/AsyCovCorr.R')
# source('E:/CurrentSimulation/RCPhi/TSCovCorr.R')
# source('E:/CurrentSimulation/RCPhi/Polychoric/functions20160810.R')
# source('E:/CurrentSimulation/RCPhi/BootJack20160816.R')


##--------------------------------------------------------------------------------------------------------
Make.Rot.Args <- function(rtype,rotation,normalize,p,m, geomin.delta, MTarget, MWeight,
                          PhiWeight, PhiTarget,transformation=NULL) {


if ((rotation=='geomin') & (is.null(geomin.delta))) geomin.delta = 0.01
if ((rotation=='target') & ((is.null(MWeight)) | (is.null(MTarget)))) stop ("MWeight or MTarget is not specified for target rotation")
if ((rotation=='xtarget') & ((is.null(MWeight)) | (is.null(MTarget)) | (is.null(PhiWeight)) | (is.null(PhiTarget)) )) stop ("MWeight or MTarget is not specified for xtarget rotation") # 2016-06-03, GZ


if (rtype=='oblique') {

  if (rotation=='CF-varimax') {

  fnames = 'cfQ'
  Rot.Args <- list(Tmat=diag(m),kappa = 1/p, normalize=normalize, eps=1e-6, maxit=1000)   

  } else if (rotation=='CF-quartimax') {

  fnames = 'cfQ'
  Rot.Args <- list(Tmat=diag(m),kappa = 0, normalize=normalize, eps=1e-6, maxit=1000)   

  } else if (rotation=='geomin') {


  fnames = 'geominQ'
  Rot.Args <- list(Tmat=diag(m),delta = geomin.delta, normalize=normalize, eps=1e-6, maxit=1000)   


  } else if (rotation=='target') {

  fnames = 'pstQ'
  Rot.Args <- list(Tmat=diag(m), W = MWeight, Target=MTarget, normalize=normalize, eps=1e-6, maxit=1000)   


  } else if (rotation=='xtarget') {

  
    fnames = 'xpstQ'
    Rot.Args <- list(Tmat=transformation, normalize=FALSE, eps=1e-6, maxit=1000,
                     method="pst",methodArgs = list(W = MWeight, Target = MTarget),PhiWeight = PhiWeight, PhiTarget = PhiTarget)   
    

  } else {
  stop (paste(rotation, ' has not been implemented yet.'))
  }



} else {   ### orthogonal rotations


  if (rotation=='CF-varimax') {

  fnames = 'cfT'
  Rot.Args <- list(Tmat=diag(m),kappa = 1/p, normalize=normalize, eps=1e-6, maxit=1000)   

  } else if (rotation=='CF-quartimax') {

  fnames = 'cfT'
  Rot.Args <- list(Tmat=diag(m),kappa = 0, normalize=normalize, eps=1e-6, maxit=1000)   

  } else if (rotation=='geomin') {


  fnames = 'geominT'
  Rot.Args <- list(Tmat=diag(m),delta = geomin.delta, normalize=normalize, eps=1e-6, maxit=1000)   


  } else if (rotation=='target') {

  fnames = 'pstQT'
  Rot.Args <- list(Tmat=diag(m), W = MWeight, Target=MTarget, normalize=normalize, eps=1e-6, maxit=1000)   

  } else {
  stop (paste(rotation, ' has not been implemented yet.'))
  }


} # End of the orthogonal rotations


list(fnames = fnames, Rot.Args = Rot.Args)

} # Make.Rot.Args

#-------------------------------------------------------------------------------------------------------

### 2016-08-12, Friday!
### The task is to plan the function Make.se.Args
### Would it be a better idea to compute SE in this function?
### How do you handle the two rotation types and different methods of compute SEs?
### Should I add a parametric bootstrap?
### How do you handle different distributions? normal, continuous, ts, and ordinal?


### How about treat se methods as the outmost cycle, then rtype, and then Dist?

# Make.se.Args <- function (x, factors, R0, n, dist, fm, rtype, rotation, normalize, geomin.delta, MTarget, MWeight,
#                           PhiWeight, PhiTarget, se, confid) {

#  if ( ( ! ( dist=='normal')) & (se =='information') ) {
#  se == 'sandwich'
#  message('The fisher information SE estimates are only for normal data; Sandwich SE are used for non-normal data.')
#  } 

#  f.se.name = 'oblq.se.augmt'
#  if (rtype=='orthogonal') f.se.name = 'orth.se.augmt'


#  if (se=='information') {
#  se.Args = list()
#  } else if (se=='sandwich') {

#  } 


#  } # Make.se.Args


#-------------------------------------------------------------------------------------------------------

### 2016-08-12, Friday!
### How do you handle the two rotation types and different methods of compute SEs?
### Should I add a parametric bootstrap?
### How do you handle different distributions? normal, continuous, ts, and ordinal?


### How about treat se methods as the outmost cycle, then rtype, and then Dist?

Compute.se <- function (x,R0, n, rotated, phi, dist, fm, rtype, rotation, normalize, geomin.delta, MTarget, MWeight,
                          PhiWeight, PhiTarget, se, confid,Ib, FE.Arg, Rot.Controls) {


if (se=='information') {

  # information deals with only normal variables and correctly specified models

  if (rtype=='oblique') {

  analytic.se = oblq.se.augmt(Lambda = rotated, Phi = phi, Rsample=R0, N=n, extraction=fm, 
                              normalize=normalize, rotation=rotation, modelerror='NO', geomin.delta = geomin.delta,
                              MTarget=MTarget, MWeight=MWeight, PhiWeight=PhiWeight, PhiTarget=PhiTarget) 
  } else { # orthogonal

  analytic.se = orth.se.augmt(Lambda = rotated, Rsample=R0, N=n, extraction=fm, 
                             normalize=normalize, rotation=rotation, modelerror='NO', geomin.delta = geomin.delta,
                             MTarget=MTarget, MWeight=MWeight) 
  } # orthogonal



} else if (se == 'sandwich') {

  # se == Sandwich deals with non-normal distributions 

  if (dist=='continuous') {
    u.r = AsyCovCorr(x)$asc

  } else if (dist=='ts') {
    u.r = TSCovCorr(x)$asc

  } else if (dist=='ordinal') {
    u.r = get.RGamma(x, gamma=TRUE)$GammaR
  }


  if (rtype=='oblique') {
 
  analytic.se = oblq.se.augmt(Lambda = rotated, Phi = phi, Rsample=R0, N=n, extraction=fm, 
                              normalize=normalize, rotation=rotation, modelerror=merror, geomin.delta = geomin.delta,
                              MTarget=MTarget, MWeight=MWeight, PhiWeight=PhiWeight, PhiTarget=PhiTarget,u.r=u.r) 
  } else { # orthogonal

  analytic.se = orth.se.augmt(Lambda = rotated, Rsample=R0, N=n, extraction=fm, 
                             normalize=normalize, rotation=rotation, modelerror=merror, geomin.delta = geomin.delta,
                             MTarget=MTarget, MWeight=MWeight, u.r = u.r) 
  } # orthogonal


} else if (se == 'jackknife') {

Jack.Arg <- list(bj='jackknife',Ib=n, rtype=rtype,dist=dist, Level.Confid=confid,FE.Arg=FE.Arg, fnames = Rot.Controls$fnames, Rotation.Arg=Rot.Controls$Rot.Args) # Correct a bug, 2016-08-26, GZ
Jack = BootJack(x,rotated,Jack.Arg) # Remove a bug 2016-08-27, GZ

} else if (se =='bootstrap') {

Boot.Arg <- list(bj='bootstrap',Ib=Ib, rtype=rtype,dist=dist, Level.Confid=confid,FE.Arg=FE.Arg, fnames = Rot.Controls$fnames, Rotation.Arg=Rot.Controls$Rot.Args)
Boot = BootJack(x,rotated,Boot.Arg) # Remove a bug 2016-08-27, GZ

} # bootstrap


# output

# analytic.se

} # Compute.se


#-------------------------------------------------------------------------------------------------------

# check input arguments
# external functions: get.RGamma, fa.extract
confid = 0.95
if ( (is.null(x)) & (is.null(covmat)) ) stop ("Neither raw data nor the correlation matrix is provided!")
if (! (dist == 'normal') & (is.null(x)) ) stop ("Raw data are required for non-normal distributions!")
if ( (is.null(x)) & (is.null(n.obs)) ) stop ("The sample size is not provided for the correlation matrix!")
###------------------------------------------------------
if (!(is.null(x))) {

p = ncol (x)
n = nrow (x)

if (n <p) stop ("The sample size is less than the number of manifest variables!")

if (dist=='ordinal') {
polychor= get.RGamma(x)
R0 = polychor$R
} else {
R0 = cor(x)
}

if (!(is.null(covmat))) {

if( min(abs(R0 - covmat)) > 0.0001) message ('covmat is different from the one computed from the raw data! The one computed from raw data is used for EFA!') 

 } # if (!(is.null(covmat)))

} else { ## (!(is.null(x)))

  p = ncol(covmat)
  n = n.obs

  mvariance <- diag(covmat)
  
if (all(mvariance==1)) {
 R0 = covmat
} else {
 message('covmat is not a correlation matrix; EFA is conducted with the corresponding correlation matrix.')
 R0 = covmat
 msd = sqrt(mvariance)
 for (i in 1:p) {
  R0[i,1:p] = R0[i,1:p] / msd[i]
  R0[1:p,i] = R0[1:p,i] / msd[i] 
  }
} # if the input matrix is a covariance matrix 
  
    
} ## is.null(x)

###### -------------------------------------------------

max.factors = floor(((2*p + 1) - sqrt(8*p+1))/2)

ev = eigen(R0)$values

  if (is.null(factors)) {factors = length(which(ev > 1))
} else {
  if (factors > max.factors) stop (paste(factors, "factors is too many for",p,"variables."))
}

## manifest variable names and factor names

if (is.null(mnames)) {
mnames = rep(" ", p)  
  for (i in 1:p) {
    mnames[i] = paste("MV",i,sep="")
  }
} else{
  if (length(mnames) != p) stop("The number of MV names is different from the number of MVs!")
  }


if (is.null(fnames)) {
  fnames = rep(" ", factors)  
  for (i in 1:factors) {
    fnames[i] = paste("F",i,sep="")
  }
} else{
  if (length(fnames) != factors) stop("The number of factor names is different from the number of factors!")
}

#### 

# extract m factors
A.lst = fa.extract(R0,factors, extraction = fm)

FE.Arg <- list(factors=factors, extraction = fm)

######


# factor rotation
# I need an envelope function to handle multiple rotation methods.

if (factors > 1) {

transformation = NULL
if (rotation=='xtarget') {
  if(rtype=='orthogonal') { 
    rtype = 'oblique'
    message('xtarget requires oblique rotation.')
    }
rotation = 'target'  
Rot.Controls <- Make.Rot.Args(rtype,rotation,normalize,p,factors,geomin.delta,MTarget, MWeight,PhiWeight, PhiTarget)
Lambda.lst = do.call (Rot.Controls$fnames, append(list(A.lst$Unrotated),Rot.Controls$Rot.Args))
transformation = (t(A.lst$Unrotated) %*% A.lst$Unrotated) %*% solve(t(Lambda.lst$loadings) %*% A.lst$Unrotated)

rotation = 'xtarget'

} # rotation = 'target'

Rot.Controls <- Make.Rot.Args(rtype,rotation,normalize,p,factors,geomin.delta,MTarget, MWeight,PhiWeight, PhiTarget,transformation)
Lambda.lst = do.call (Rot.Controls$fnames, append(list(A.lst$Unrotated),Rot.Controls$Rot.Args))



if ((rotation == 'target') | (rotation =='xtarget')) useorder = FALSE

if (useorder) {

if ( ( is.null(MWeight) ) | ( is.null(MTarget) ) ) stop ("MWeight or MTarget is not specified when ordering is request.")

M.in.temp = matrix(0,(p+factors),factors)
M.in.temp [1:p,1:factors] = Lambda.lst$loadings

if (rtype=='oblique') {
M.in.temp [(p+1):(p+factors),1:factors] = Lambda.lst$Phi
} else{
M.in.temp [(p+1):(p+factors),1:factors] = diag(factors)
}

M.out.temp = Align.Matrix (MTarget, M.in.temp, MWeight)

Lambda.lst$loadings = M.out.temp[1:p,1:factors]
if (rtype=='oblique') Lambda.lst$Phi = M.out.temp[(p+1):(p+factors),1:factors]

} # (useorder)

} # if (factors > 1)


########


# standard errors
# I need an envelope function to handle multiple standard error procedures.

## modify se if necessary
if ( ( ! ( dist=='normal')) & (se =='information') ) {
  se = 'sandwich'
  message('The fisher information SE estimates are only for normal data; Sandwich SE are used for non-normal data.')
} 


if (  ( dist=='normal') & (se =='sandwich') ) {
  se = 'information'
  message('The fisher information SE estimates are computed for normal data.')
} 

if (  ( (dist=='ordinal') | (dist=='ts')) & (se =='jackknife') ) {
  se = 'bootstrap'
  message('The jackknife SE estimates for ordinal data and time series data have not been developped yet; bootstrap SE estimates are computed instead.')
} 
###

Phi = diag(factors)


if (factors > 1) {
  if (rtype == 'oblique') Phi = Lambda.lst$Phi
  
  SE = Compute.se (x, R0, n, rotated=Lambda.lst$loadings, phi=Phi, dist, fm, rtype, rotation, normalize, geomin.delta, MTarget, MWeight,
                          PhiWeight, PhiTarget, se, confid, Ib, FE.Arg, Rot.Controls)
} else {
SE = Compute.se (x, R0, n, rotated=A.lst$Unrotated, phi=Phi, dist, fm, rtype='orthogonal', rotation='unrotated', normalize, geomin.delta, MTarget, MWeight,
                          PhiWeight, PhiTarget, se, confid, Ib, FE.Arg, Rot.Controls)
}

## Outputs
details = list(manifest=p,factors=factors, n.obs=n, dist=dist, fm=fm, rtype=rtype, rotation=rotation, normalize=normalize, 
           geomin.delta=geomin.delta, MTarget=MTarget, MWeight=MWeight,
                          PhiWeight = PhiWeight, PhiTarget = PhiTarget, se=se, confid=confid, Ib=Ib)

### 

unrotated = A.lst$Unrotated
fdiscrepancy = A.lst$f
convergence = A.lst$convergence
heywood = A.lst$heywood

if (factors > 1) {
rotated = Lambda.lst$loadings
} else {
rotated = A.lst$Unrotated} 

Phi = Phi

if ((se=='bootstrap') | (se=='jackknife') ) {
  rotatedse = SE$SE.LPhi[1:p,1:factors]
  Phise = SE$SE.LPhi[(p+1):(p+factors),1:factors]
}


if ((se=='information') | (se=='sandwich') ) {
  rotatedse = SE$Lambda.se[1:p,1:factors]
  Phise = matrix(0, factors, factors)
  if (rtype=='oblique') Phise =  SE$Phi.se[1:factors,1:factors]
}

if (factors==1) {
  rotatedse = matrix(rotatedse,nrow=p,ncol=factors)
  Phise = matrix(0, factors, factors)
}

dimnames(unrotated) = list(mnames,fnames)
dimnames(rotated) = list(mnames,fnames)
dimnames(rotatedse) = list(mnames,fnames)

dimnames(Phi) = list(fnames,fnames)
dimnames(Phise) = list(fnames,fnames)

efaout = list(details = details, unrotated=unrotated,fdiscrepancy=fdiscrepancy,convergence=convergence,heywood=heywood,
              rotated=rotated,Phi=Phi, rotatedse=rotatedse, Phise= Phise)

class(efaout) <- "efa" 

return(efaout)

# confidence intervals


# output

# unroated
# f
# rotated
# phi
# rotated.se
# phi.se
# rotated.confid
# phi.confid

} # efa

print.efa <- function(x, ...) {

  cat("\nAnalysis Details: \n")
  cat("   Number of Manifest Variable:  ", x$details$manifest, "\n")
  cat("   Number of Factors:  ",x$details$factors,"\n")
  cat("   Sample Size:  ",x$details$n.obs,"\n")
  cat("   Manifest Variable Distribution:  ",x$details$dist,"\n")  
  cat("   Estimation Method:  ",x$details$fm,"\n")
  cat("   Rotation Type:  ",x$details$rtype,"\n")  
  cat("   Rotation Criterion:  ",x$details$rotation,"\n") 
  cat("   Rotation Standardization:  ",x$details$normalize,"\n")   
  cat("   Standard Error:  ",x$details$se,"\n")  

  cat("\nUnrotated Factor Loadings: \n")
  print(round(x$unrotated,3))
  
  cat("\nRotated Factor Loadings: \n")
  print(round(x$rotated,3))
  
  cat("\nFactor Correlations: \n")
  print(round(x$Phi,3))
  
  cat("\nSE for rotated Factor Loadings: \n")
  print(round(x$rotatedse,3))
  
  cat("\nSE for Factor Correlations: \n")
  print(round(x$Phise,3))
  
} # print.efa 


summary.efa <- function(object, ...){

res <- object

class(res) <- "summary.efa"

return(res)

} # summary.efa


print.summary.efa <- function(x, ...) {

  cat("\nAnalysis Details: \n")
  cat("   Number of Manifest Variable:  ", x$details$manifest, "\n")
  cat("   Number of Factors:  ",x$details$factors,"\n")
  cat("   Sample Size:  ",x$details$n.obs,"\n")
  cat("   Manifest Variable Distribution:  ",x$details$dist,"\n")  
  cat("   Estimation Method:  ",x$details$fm,"\n")
  cat("   Rotation Type:  ",x$details$rtype,"\n")  
  cat("   Rotation Criterion:  ",x$details$rotation,"\n") 
  cat("   Rotation Standardization:  ",x$details$normalize,"\n")   
  cat("   Standard Error:  ",x$details$se,"\n")  

  cat("\nUnrotated Factor Loadings: \n")
  print(round(x$unrotated,3))
  
  cat("\nRotated Factor Loadings: \n")
  print(round(x$rotated,3))
  
  cat("\nFactor Correlations: \n")
  print(round(x$Phi,3))
  
  cat("\nSE for rotated Factor Loadings: \n")
  print(round(x$rotatedse,3))
  
  cat("\nSE for Factor Correlations: \n")
  print(round(x$Phise,3))

} # print.summary.efa

UseMethod('summary',efa)