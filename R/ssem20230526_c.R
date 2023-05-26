ssem <- function(x=NULL, factors=NULL, exfactors=1, covmat=NULL, acm=NULL, n.obs=NULL, dist='normal', fm='ml', mtest = TRUE, rotation='semtarget', normalize=FALSE, maxit=1000, geomin.delta=NULL, MTarget=NULL, MWeight=NULL,
                          BGWeight = NULL, BGTarget = NULL, PhiWeight = NULL, PhiTarget = NULL, useorder=TRUE, se='sandwich', LConfid=c(0.95,0.90), CItype='pse', Ib=2000, mnames=NULL, fnames=NULL, merror='YES', wxt2 = 1e0) {

## Internal functions: Make.Rot.Args
## Extenral functions:


##--------------------------------------------------------------------------------------------------------
Make.Rot.Args <- function(rtype,rotation,normalize,p,m, geomin.delta, MTarget, MWeight, BGTarget, BGWeight,
                          PhiTarget,PhiWeight,wxt2=1e0,transformation=NULL) {



if (rtype=='oblique') {

  if (rotation=='CF-varimax') {

  fnames = 'cfQ'
  Rot.Args <- list(Tmat=diag(m),kappa = 1/p, normalize=normalize, eps=1e-6, maxit=maxit)

  } else if (rotation=='CF-quartimax') {

  fnames = 'cfQ'
  Rot.Args <- list(Tmat=diag(m),kappa = 0, normalize=normalize, eps=1e-6, maxit=maxit)

  } else if (rotation=='CF-facparsim') {

    fnames = 'cfQ'
    Rot.Args <- list(Tmat=diag(m),kappa = 1, normalize=normalize, eps=1e-6, maxit=maxit)


  } else if (rotation=='CF-equamax') {

    fnames = 'cfQ'
    Rot.Args <- list(Tmat=diag(m),kappa = m/(2*p), normalize=normalize, eps=1e-6, maxit=maxit)


  } else if (rotation=='CF-parsimax') {

    fnames = 'cfQ'
    Rot.Args <- list(Tmat=diag(m),kappa = (m-1)/(p+m-2), normalize=normalize, eps=1e-6, maxit=maxit)


  } else if (rotation=='geomin') {


  fnames = 'geominQ'
  Rot.Args <- list(Tmat=diag(m),delta = geomin.delta, normalize=normalize, eps=1e-6, maxit=maxit)


  } else if (rotation=='target') {

  fnames = 'pstQ'
  Rot.Args <- list(Tmat=diag(m), W = MWeight, Target=MTarget, normalize=normalize, eps=1e-6, maxit=maxit)


  } else if (rotation=='semtarget') {


    fnames = 'ESEMpstQ'
    Rot.Args <- list(Tmat=transformation, normalize=normalize, eps=1e-6, maxit=maxit,
                     method="pst",methodArgs = list(W = MWeight, Target = MTarget),BGTarget=BGTarget, BGWeight = BGWeight, PhiTarget = PhiTarget, PhiWeight = PhiWeight,  wxt2 = wxt2)   ## 2018-08-14, GZ ####



  } else {
  stop (paste(rotation, ' has not been implemented yet.'))
  }



} else {   ### orthogonal rotations


  stop ('Only oblique rotation is allowed in sSEM.')


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
                          BGTarget, BGWeight, PhiTarget, PhiWeight, se, confid,Ib, FE.Arg, Rot.Controls,acm.type,wxt2=1e0) {


if (se=='information') {

  # information deals with only normal variables and correctly specified models

  if (rtype=='oblique') {

  analytic.se = ssem.se.augmt(Lambda = rotated, Phi = phi, Rsample=R0, N=n, extraction=fm,
                              normalize=normalize, rotation=rotation, modelerror='NO', geomin.delta = geomin.delta,
                              MTarget=MTarget, MWeight=MWeight, BGTarget=BGTarget, BGWeight=BGWeight, PhiTarget=PhiTarget, PhiWeight=PhiWeight, acm.type=acm.type, wxt2=wxt2) # 2018-08-14, GZ
  } #



} else if (se == 'sandwich') {

  # se == Sandwich deals with non-normal distributions

  if (!(is.null(acm))) u.r = acm  # 2020-05-12, GZ



  if (dist=='continuous') {
    if (is.null(acm)) u.r = AsyCovCorr(x)$asc # 2020-05-12, GZ

  } else if (dist=='ts') {
    if (is.null(acm)) u.r = TSCovCorr(x)$asc # 2020-05-12, GZ

  } else if (dist=='ordinal') {
#    if (is.null(acm)) u.r = get.RGamma(x, gamma=TRUE)$GammaR # 2020-05-12, GZ

  } else if (dist=='normal') {

    if (is.null(acm)) u.r = EliU(R0) # 2020-05-12, GZ
  } else {
    stop (paste(dist, ' is not recognized as a data type. Four types of data are allowed: continuous, ts, ordinal, and normal.'))
  }


  if (rtype=='oblique') {

  analytic.se = ssem.se.augmt(Lambda = rotated, Phi = phi, Rsample=R0, N=n, extraction=fm,
                              normalize=normalize, rotation=rotation, modelerror=merror, geomin.delta = geomin.delta,
                              MTarget=MTarget, MWeight=MWeight, BGTarget=BGTarget, BGWeight=BGWeight, PhiWeight=PhiWeight, PhiTarget=PhiTarget,u.r=u.r, acm.type=acm.type, wxt2=wxt2)
  }


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

rtype='oblique' # only oblique rotation allowed in SSEM.

if (factors < 2) {
  stop ('SSEM requires at least two factors.')

}


confid = LConfid[1]

if ( (is.null(x)) & (is.null(covmat)) ) stop ("Neither raw data nor the correlation matrix is provided!")
if ( ! (is.null(acm)) & (is.null(covmat)) ) stop ("specifying the acm requires the correlation matrix.")
if (! (dist == 'normal') & (is.null(x)) ) stop ("Raw data are required for non-normal distributions!")
if ( (is.null(x)) & (is.null(n.obs)) ) stop ("The sample size is not provided for the correlation matrix!")
if ( !(se=='bootstrap') & (CItype=='percentile') ) {
 CItype='pse'
 message ('Percentile Confidence intervals are avaible only for bootstrap: pse confidence intervals are constructed.')
}


### 2018-08-14, GZ

if (exfactors > factors) {
  stop (paste(factors, "factors are fewer than ",exfactors," exogenous factors."))
}

enfactors = factors - exfactors

if (enfactors==0) stop("Because all factors are exogenous, the model is a factor analysis model. You should estimate it using efa()") # 2022-07-11, GZ


### 2018-08-14, GZ


###------------------------------------------------------
# Start of 2020-05-28, GZ


positiveloadings = rep(0,2) # A possible fake variable to facilitate column reflection

if (useorder) {
  if (is.null(MTarget)) {
    stop ("The Order Matrix (MTarget) is not specified when ordering is request.")
  } else {

    p = nrow(MTarget)
    m = ncol(MTarget)
    MTarget.p = MTarget
    MTarget.p[is.na(MTarget)] <- - 9

    positiveloadings = rep(0,m)

    for (j in 1:m) {
      for (i in 1:p) {

        if (MTarget.p[i,j]==9) {
          positiveloadings[j] = i
          break
        }
      }
    }

    MTarget[is.na(MTarget)] <- 9


    if (is.null(MWeight)) {

      MWeight = MTarget

      MWeight[MTarget != 9] <- 1 # 1 corresponds to small (zero) loadings
      MWeight[MTarget == 9] <- 0 # 0 corresponds to large loadings


    } #  if (is.null(MWeight))


  }

} # if (useorder)




# making MWeight matrix optional

if ((rotation=='geomin') & (is.null(geomin.delta))) geomin.delta = 0.01


if (rotation=='target') {

  if (is.null(MTarget)) {
    stop ("MTarget is not specified for target rotation")
  } else {

    MTarget[is.na(MTarget)] <- 9


    if (is.null(MWeight)) {

      MWeight = MTarget

      MWeight[MTarget != 9] <- 1 # 1 corresponds to small (zero) loadings
      MWeight[MTarget == 9] <- 0 # 0 corresponds to large loadings


    } #  if (is.null(MWeight))

  }

} # (rotation=='target')


if (rotation=='semtarget') {

  if ( is.null(MTarget) & (is.null(PhiTarget)) &  (is.null(BGTarget))) {
    stop ("MTarget or PhiTarget is not specified for semtarget rotation")
  } else {

    MTarget[is.na(MTarget)] <- 9
    PhiTarget[is.na(PhiTarget)] <- 9
    BGTarget[is.na(BGTarget)] <- 9

    if (is.null(MWeight)) {

      MWeight = MTarget

      MWeight[MTarget != 9] <- 1 # 1 corresponds to small (zero) loadings
      MWeight[MTarget == 9] <- 0 # 0 corresponds to large loadings


    } #  if (is.null(MWeight))

    if (is.null(PhiWeight)) {

      PhiWeight = PhiTarget

      PhiWeight[PhiTarget != 9] <- 1 # 1 corresponds to small (zero) factor correlations
      PhiWeight[PhiTarget == 9] <- 0 # 0 corresponds to large factor correlations

      diag(PhiWeight) = 0 # diagonal elements having 0 weights. Note that it is slower than assigning
                        # values with a for loop

    } #  if (is.null(PhiWeight))


    if (is.null(BGWeight)) {

      BGWeight = BGTarget

      BGWeight[BGTarget != 9] <- 1 # 1 corresponds to small (zero) factor correlations
      BGWeight[BGTarget == 9] <- 0 # 0 corresponds to large factor correlations

      enfactors = nrow(BGTarget)

      for (j in 1:enfactors) { # diagonal and lower diagonal elements of BGTarget having 0 weights.
        BGWeight[j:enfactors,j] = 0
      }


    } #  if (is.null(PhiWeight))



  }

} # if (is.null(MTarget))





# End of 2020-05-28, GZ


###------------------------------------------------------

if (!(is.null(x))) {

x = data.matrix(x) # 2018-12-23
p = ncol (x)
n = nrow (x)

if (n <p) stop ("The sample size is less than the number of manifest variables!")

if ((dist=='ordinal') & (is.null(acm))) {

  if ((! mtest) & ( se == 'none'))  {
    # polychor =   PolychoricRM(x, NCore=4, IAdjust=1)
    polychor = get.RGamma(x, gamma=FALSE)
    R0 = polychor$R
  } else {

    polychor = get.RGamma(x, gamma=TRUE)
    # polychor =  PolychoricRM(x, NCore=4, IAdjust=1, estimate.acm=TRUE)
    R0 = polychor$R
    u.r = polychor$GammaR
  }

} else { ## Not ordinal variables
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

#### reconciling two acm types
if (is.null(acm)) {
  acm.type=2
} else {
  acm.type=1
}

# if (dist=='ordinal') acm.type=1 # 2020-07-08, GZ
### end of reconciling two acm types


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


### 2017-08-08, adding test statistics

# if (A.lst$heywood==0) {

if (mtest) { # 2020-05-12, GZ


if (fm=='ml') {
  statistic = (n-1) * A.lst$f
  df = ((p - factors)**2 - p - factors ) /2         # 2019-01-19, GZ
  stat.stage1 = list(statistic=statistic, df = df)  # 2019-01-19, GZ

} else if (fm=='ols') {

  if (!(is.null(acm))) u.r = acm # 2020-05-12, GZ

  if (dist=='continuous') {
    if (is.null(acm)) u.r = AsyCovCorr(x)$asc  # 2020-05-12, GZ

  } else if (dist=='ts') {
    if (is.null(acm)) u.r = TSCovCorr(x)$asc # 2020-05-12, GZ

  } else if (dist=='ordinal') {
#    if (is.null(acm)) u.r = get.RGamma(x, gamma=TRUE)$GammaR # 2020-05-12, GZ

  } else if (dist=='normal') {
    if (is.null(acm)) u.r = EliU(R0) # 2020-05-12, GZ

  }
  stat.stage1 = Compute.stat(R0,u.r,A.lst$Unrotated,acm.type)
  statistic = stat.stage1$statistic

  statistic = ifelse(is.nan(statistic), NaN, statistic*(n-1))


} # ols

} # if (mtest) 2020-05-12, GZ


if (! mtest) stat.stage1 = list(statistic = 10, df = 10 ) # 2020-05-12, GZ, the code is added to accommodate the next line
                                                          # it does not do anything


if ((stat.stage1$df <= 0) | (! mtest)) { # 2020-05-12

  #### 2018-12-31
  if (! mtest) message('The EFA statistic is not requested.')
  if ( mtest) message('The EFA test statistic is invalid because the degrees of freedom are not positive.')

  statistic = p*(p-1)/2
  ModelF = Model.Fit(statistic,fm,p,factors,n,LConfid[2])

  ModelF$f.stat = 0
  ModelF$RMSEA = NaN
  ModelF$p.perfect = NaN
  ModelF$p.close = NaN
  ModelF$RMSEA.l = NaN
  ModelF$RMSEA.u = NaN
  ModelF$ECVI = NaN
  ModelF$ECVI.l = NaN
  ModelF$ECVI.u = NaN


} else{
  if (is.nan(statistic)) {
    message('The EFA test statistic is invalid because the estimate of the asymptotic covariance matrix of correlations has negative eigenvalues.')

    statistic = p*(p-1)/2
    ModelF = Model.Fit(statistic,fm,p,factors,n,LConfid[2])

    ModelF$f.stat = NaN
    ModelF$RMSEA = NaN
    ModelF$p.perfect = NaN
    ModelF$p.close = NaN
    ModelF$RMSEA.l = NaN
    ModelF$RMSEA.u = NaN
    ModelF$ECVI = NaN
    ModelF$ECVI.l = NaN
    ModelF$ECVI.u = NaN

  }else {
    ModelF = Model.Fit(statistic,fm,p,factors,n,LConfid[2])
  }
}  # if ((stat.stage1$df <= 0) | (! mtest)) { # 2020-05-12

#} else { # A.lst$heywood > 0, 2017-08-14

#  message('The EFA test statistic is invalid because a Heywood case occurs.')


# }  # A.lst$heywood > 0, 2017-08-14


### 2017-08-08


# factor rotation
# I need an envelope function to handle multiple rotation methods.

## I replaced 'xtarget' with 'semtarget'.

if (factors > 1) {

transformation = NULL
if (rotation=='semtarget') {  # 2018-08-14, GZ, replacing 'xtarget' with 'semtarget'
  if(rtype=='orthogonal') {
    rtype = 'oblique'
    message('semtarget requires oblique rotation.')
    }
rotation = 'target'
Rot.Controls <- Make.Rot.Args(rtype,rotation,normalize,p,factors,geomin.delta,MTarget, MWeight,BGTarget, BGWeight, PhiTarget, PhiWeight,wxt2)
Lambda.lst = do.call (Rot.Controls$fnames, append(list(A.lst$Unrotated),Rot.Controls$Rot.Args))
transformation = (t(A.lst$Unrotated) %*% A.lst$Unrotated) %*% solve(t(Lambda.lst$loadings) %*% A.lst$Unrotated)

rotation = 'semtarget'

} # rotation = 'semtarget'

Rot.Controls <- Make.Rot.Args(rtype,rotation,normalize,p,factors,geomin.delta,MTarget, MWeight,BGTarget, BGWeight, PhiTarget, PhiWeight, wxt2,transformation)
Lambda.lst = do.call (Rot.Controls$fnames, append(list(A.lst$Unrotated),Rot.Controls$Rot.Args))



# if ((rotation == 'target') | (rotation =='semtarget')) useorder = FALSE # 2018-08-31, Friday, GZ

if (useorder) {

M.in.temp = matrix(0,(p+factors),factors)
M.in.temp [1:p,1:factors] = Lambda.lst$loadings

if (rtype=='oblique') {
M.in.temp [(p+1):(p+factors),1:factors] = Lambda.lst$Phi
} else{
M.in.temp [(p+1):(p+factors),1:factors] = diag(factors)
}

M.out.temp = Align.Matrix (MTarget, M.in.temp) # I have removed another bug, 2018-08-31, GZ

Lambda.lst$loadings = M.out.temp[1:p,1:factors]
if (rtype=='oblique') Lambda.lst$Phi = M.out.temp[(p+1):(p+factors),1:factors]


if (sum(positiveloadings)>0) {

  for (j in 1:factors) {
    if (positiveloadings[j] > 0) {
      if (Lambda.lst$loadings[positiveloadings[j],j] < 0) {
        Lambda.lst$loadings[1:p,j] = Lambda.lst$loadings[1:p,j] * (-1)
        Lambda.lst$Phi[1:factors,j] = Lambda.lst$Phi[1:factors,j] * (-1)
        Lambda.lst$Phi[j,1:factors] = Lambda.lst$Phi[j,1:factors] * (-1)
      }
    }
  } # (j in 1:factors)

} # if (sum(positiveloadings)>0)



} # (useorder)


} # if (factors > 1)


########


# standard errors
# I need an envelope function to handle multiple standard error procedures.

## modify se if necessary
if ( ( ! (( dist=='normal') & (merror=='NO'))) & (se =='information') ) {
  se = 'sandwich'
  message('The fisher information SEs are valid only for normal data and no model error. Sandwich SEs are computed.' )
}


# if (  ( dist=='normal') & (se =='sandwich') ) {
#  se = 'information'
#  message('The fisher information SE estimates are computed for normal data.')
# }

if (  ( (dist=='ordinal') | (dist=='ts')) & (se =='jackknife') ) {
  se = 'bootstrap'
  message('The jackknife SE estimates for ordinal data and time series data have not been developped yet; bootstrap SE estimates are computed instead.')
}
###

Phi = diag(factors)


if (factors > 1) {
  if (rtype == 'oblique') Phi = Lambda.lst$Phi

  SE = Compute.se (x, R0, n, rotated=Lambda.lst$loadings, phi=Phi, dist, fm, rtype, rotation, normalize, geomin.delta, MTarget, MWeight,
                          BGTarget, BGWeight, PhiTarget, PhiWeight, se, confid, Ib, FE.Arg, Rot.Controls,acm.type,wxt2)
}


## Problem 1: SE estimates for methods other than information
## Problem 2: the rotated BG

## Outputs
acm.in = FALSE
if (! (is.null(acm))) acm.in = TRUE


if (! (is.null(MTarget))) MTarget[MTarget==9] = NA
if (! (is.null(BGTarget))) BGTarget[BGTarget==9] = NA
if (! (is.null(PhiTarget))) PhiTarget[PhiTarget==9] = NA

if (sum(positiveloadings)>0) {
  for (j in 1:factors) {
    if (positiveloadings[j]==0) next
    MTarget[positiveloadings[j],j] = 9
  }
}


details = list(manifest=p,factors=factors, enfactors = factors - exfactors, exfactors = exfactors,n.obs=n, dist=dist, acm.in = acm.in, fm=fm, rtype=rtype, rotation=rotation, normalize=normalize,
           geomin.delta=geomin.delta, wxt2=wxt2, MTarget=MTarget, MWeight=MWeight, BGTarget=BGTarget, BGWeight=BGWeight,
                          PhiWeight = PhiWeight, PhiTarget = PhiTarget, merror=merror, mtest = mtest, se=se, LConfid=LConfid, Ib=Ib)

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
  Phise =  SE$Phi.se[1:factors,1:factors]
  BG = SE$BG
  BG.se = SE$BG.se
  psi = SE$psi
  psi.se = SE$psi.se
  Phi.xi = SE$Phi.xi
  Phi.xi.se = SE$Phi.xi.se
}

if (factors==1) {
  rotatedse = matrix(rotatedse,nrow=p,ncol=factors)
  Phise = matrix(0, factors, factors)
}

#### 2017-08-08, compute confidence intervals

alpha = 1 - LConfid[1]
CI.lambda = CIs(rotated, rotatedse, alpha, type = 'lambda.oblq')
CI.Phi = CIs(Phi, Phise, alpha, type = 'Phi')

rotatedlow = CI.lambda$LowerLimit
rotatedupper = CI.lambda$UpperLimit
Philow = CI.Phi$LowerLimit
Phiupper = CI.Phi$UpperLimit

if (enfactors>0) {
CI.BG = CIs(BG, BG.se, alpha, type='lambda.oblq') # 2018-08-14, GZ
CI.psi = CIs(psi, psi.se, alpha, type = 'uv') # 2018-08-15, GZ

BGlow = CI.BG$LowerLimit
BGupper = CI.BG$UpperLimit

psilow = CI.psi$LowerLimit
psiupper = CI.psi$UpperLimit


dimnames(BG) = list(fnames[1:enfactors],fnames)
dimnames(BG.se) = list(fnames[1:enfactors],fnames)
dimnames(BGlow) = list(fnames[1:enfactors],fnames)
dimnames(BGupper) = list(fnames[1:enfactors],fnames)


names(psi) = fnames[1:enfactors]
names(psi.se) = fnames[1:enfactors]
names(psilow) = fnames[1:enfactors]
names(psiupper) = fnames[1:enfactors]

} else {

BG = NULL
BG.se = NULL
BGlow = NULL
BGupper = NULL
psi = NULL
psi.se = NULL
psilow = NULL
psiupper = NULL

} # (enfactors>0)

if (exfactors>1) {
CI.Phi.xi = CIs(Phi.xi, Phi.xi.se, alpha, type = 'Phi') # 2018-08-14, GZ
Phixilow = CI.Phi.xi$LowerLimit
Phixiupper = CI.Phi.xi$UpperLimit

} else {

Phi.xi = diag(1)
Phi.xi.se = matrix(0,1,1)
Phixilow =diag(1)
Phixiupper = diag(1)

} # (exfactors>1)

dimnames(Phi.xi) = list(fnames[(enfactors+1):factors],fnames[(enfactors+1):factors])
dimnames(Phi.xi.se) = list(fnames[(enfactors+1):factors],fnames[(enfactors+1):factors])
dimnames(Phixilow) = list(fnames[(enfactors+1):factors],fnames[(enfactors+1):factors])
dimnames(Phixiupper) = list(fnames[(enfactors+1):factors],fnames[(enfactors+1):factors])



Phat = unrotated %*% t(unrotated)
Phat = Phat - diag(diag(Phat)) + diag(p)
Residual = R0 - Phat


#### 2017-08-08


rownames(A.lst$compsi) = mnames

dimnames(unrotated) = list(mnames,fnames)
dimnames(rotated) = list(mnames,fnames)
dimnames(rotatedse) = list(mnames,fnames)
dimnames(rotatedlow) = list(mnames,fnames)
dimnames(rotatedupper) = list(mnames,fnames)




dimnames(Phi) = list(fnames,fnames)
dimnames(Phise) = list(fnames,fnames)
dimnames(Philow) = list(fnames,fnames)
dimnames(Phiupper) = list(fnames,fnames)




dimnames(R0) = list(mnames,mnames)
dimnames(Phat) = list(mnames,mnames)
dimnames(Residual) = list(mnames,mnames)


## 2018-08-14, Tuesday!

## The task is to modify efaout to ssemout.

ssemout = list(details = details, unrotated=unrotated,fdiscrepancy=fdiscrepancy,convergence=convergence,heywood=heywood, nq = (p*(factors+1) - factors*(factors-1)/2), compsi=A.lst$compsi, R0=R0, Phat=Phat, Residual=Residual,
              rotated=rotated,Phi=Phi, BG = BG, psi = psi, Phi.xi = Phi.xi, rotatedse=rotatedse, Phise= Phise, BGse = BG.se, psise = psi.se, Phi.xise = Phi.xi.se,
              ModelF = ModelF, rotatedlow=rotatedlow, rotatedupper=rotatedupper, Philow=Philow, Phiupper=Phiupper,
              BGlow = BGlow, BGupper=BGupper, psilow=psilow, psiupper = psiupper, Phixilow = Phixilow, Phixiupper = Phixiupper)

class(ssemout) <- "ssem" ## I need to think about this: Should I change it to ssem? 2018-08-14, GZ.

return(ssemout)

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

} # ssem

#'@export

print.ssem <- function(x, moredetail=FALSE, ...) {

if (moredetail) {
  cat("\nAnalysis Details: \n")
  cat("   Number of Manifest Variable:  ", x$details$manifest, "\n")
  cat("   Number of Exogenous Factors:  ",x$details$exfactors,"\n")
  cat("   Number of Endogenous Factors:  ",x$details$enfactors,"\n")
  cat("   Sample Size:  ",x$details$n.obs,"\n")
  cat("   Manifest Variable Distribution:  ",x$details$dist,"\n")
  cat("   Estimation Method:  ",x$details$fm,"\n")
  cat("   Rotation Type:  ",x$details$rtype,"\n")
  cat("   Rotation Criterion:  ",x$details$rotation,"\n")
  cat("   Rotation Standardization:  ",x$details$normalize,"\n")
  cat("   Standard Error:  ",x$details$se,"\n")


  cat("\nMeasures of Fit \n")
  cat(" Root Mean Square Error of Approximation, RMSEA \n")
  cat("   Point estimate:  ",round(x$ModelF$RMSEA,3),"\n")
  cat("  ",x$details$LConfid[2]*100,"% Confidence Intervals:  (",round(x$ModelF$RMSEA.l,3), ",",round(x$ModelF$RMSEA.u,3), ") \n")
  cat(" Test Statistic:  ",round(x$ModelF$f.stat,3),"\n")
  cat("   Degrees of Freedom:  ",x$ModelF$df,"\n")
  cat("   Effect numbers of Parameters:  ",x$nq,"\n")
  cat("   P value for perfect fit:  ",round(x$ModelF$p.perfect,3),"\n")
  cat("   P value for close fit:  ",round(x$ModelF$p.close,3),"\n")

  cat("\nEigenvalues, SMCs, Communalities, and Unique Variance: \n")
  print(round(x$compsi,3))


  cat("\nUnrotated Factor Loadings: \n")
  print(round(x$unrotated,3))

  cat("\nRotated Factor Loadings: \n")
  print(round(x$rotated,3))


  cat("\nLatent Regression Weights: \n")
  print(round(x$BG,3))

  cat("\nLatent Residual Variances: \n")
  print(round(x$psi,3))

  cat("\nExogenous Factor Correlations: \n")
  print(round(x$Phi.xi,3))



  #  cat("\nFactor Correlations: \n")
  #  print(round(x$Phi,3))

  cat("\nSE for rotated Factor Loadings: \n")
  print(round(x$rotatedse,3))

  #  cat("\nSE for Factor Correlations: \n")
  #  print(round(x$Phise,3))


  cat("\nSE for Latent Regression Weights: \n")
  print(round(x$BGse,3))

  cat("\nSE for Latent Residual Variances: \n")
  print(round(x$psise,3))

  cat("\nSE for Exogenous Factor Correlations: \n")
  print(round(x$Phi.xise,3))


  cat("\nLower Bounds of", x$details$LConfid[1]*100,"% CIs for Rotated Factor Loadings: \n")
  print(round(x$rotatedlow,3))

  cat("\nUpper Bounds of", x$details$LConfid[1]*100,"% CIs for Rotated Factor Loadings: \n")
  print(round(x$rotatedupper,3))

  cat("\nLower Bounds of", x$details$LConfid[1]*100,"% CIs for Latent Regression Weights: \n")
  print(round(x$BGlow,3))

  cat("\nUpper Bounds of", x$details$LConfid[1]*100,"% CIs for Latent Regression Weights: \n")
  print(round(x$BGupper,3))

  cat("\nLower Bounds of", x$details$LConfid[1]*100,"% CIs for Residual Variances: \n")
  print(round(x$psilow,3))

  cat("\nUpper Bounds of", x$details$LConfid[1]*100,"% CIs for Residual Variances: \n")
  print(round(x$psiupper,3))


  cat("\nLower Bounds of", x$details$LConfid[1]*100,"% CIs for Exogenous Factor Correlations: \n")
  print(round(x$Phixilow,3))

  cat("\nUpper Bounds of", x$details$LConfid[1]*100,"% CIs for Exogenous Correlations: \n")
  print(round(x$Phixiupper,3))

} else {

  cat("\nSummary of Analysis: \n")
  cat("   Estimation Method:  ",x$details$fm,"\n")
  cat("   Rotation Criterion:  ",x$details$rotation,"\n")
  cat("   Test Statistic:  ",round(x$ModelF$f.stat,3),"\n")
  cat("   Degrees of Freedom:  ",x$ModelF$df,"\n")
  cat("   Effect numbers of Parameters:  ",x$nq,"\n")
  cat("   P value for perfect fit:  ",round(x$ModelF$p.perfect,3),"\n")


  cat("\nRotated Factor Loadings: \n")
  print(round(x$rotated,3))

  cat("\nLatent Regression Weights: \n")
  print(round(x$BG,3))

  cat("\nLatent Residual Variances: \n")
  print(round(x$psi,3))

  cat("\nExogenous Factor Correlations: \n")
  print(round(x$Phi.xi,3))
}

} # print.ssem

# summary.ssem <- function(object, ...){

# res <- object

# class(res) <- "summary.ssem"

# return(res)

# } # summary.efa

# print.summary.ssem <- function(x, ...) {



# } # print.summary.ssem

