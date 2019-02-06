#'@importFrom graphics title

#'@importFrom stats pchisq

#'@importFrom stats rnorm

#'@importFrom stats uniroot

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

#'@importFrom plyr count

#'@importFrom graphics barplot

#'@export efa

#'@export efaMR

#'@export Align.Matrix

#'@method print efa

#'@export



efa <- function(x=NULL, factors=NULL, covmat=NULL, n.obs=NULL, dist='normal', fm='ols', rtype='oblique', rotation='CF-varimax', normalize=FALSE, maxit=1000, geomin.delta=NULL, MTarget=NULL, MWeight=NULL,
PhiWeight = NULL, PhiTarget = NULL, useorder=FALSE, se='sandwich', LConfid=c(0.95,0.90), CItype='pse', Ib=2000, mnames=NULL, fnames=NULL, merror='YES', wxt2 = 1e0) {
    
    
    ##--------------------------------------------------------------------------------------------------------
    Make.Rot.Args <- function(rtype,rotation,normalize,p,m, geomin.delta, MTarget, MWeight,
    PhiWeight, PhiTarget,wxt2=1e0,transformation=NULL) {
        
        
        if ((rotation=='geomin') & (is.null(geomin.delta))) geomin.delta = 0.01
        if ((rotation=='target') & ((is.null(MWeight)) | (is.null(MTarget)))) stop ("MWeight or MTarget is not specified for target rotation")
        if ((rotation=='xtarget') & ((is.null(MWeight)) | (is.null(MTarget)) | (is.null(PhiWeight)) | (is.null(PhiTarget)) )) stop ("MWeight or MTarget is not specified for xtarget rotation")
        
        
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
                
                
            } else if (rotation=='xtarget') {
                
                
                fnames = 'xpstQ'
                Rot.Args <- list(Tmat=transformation, normalize=normalize, eps=1e-6, maxit=maxit,
                method="pst",methodArgs = list(W = MWeight, Target = MTarget),PhiWeight = PhiWeight, PhiTarget = PhiTarget, wxt2 = wxt2)
                
                
            } else {
                stop (paste(rotation, ' has not been implemented yet.'))
            }
            
            
            
        } else {   ### orthogonal rotations
            
            
            if (rotation=='CF-varimax') {
                
                fnames = 'cfT'
                Rot.Args <- list(Tmat=diag(m),kappa = 1/p, normalize=normalize, eps=1e-6, maxit=maxit)
                
            } else if (rotation=='CF-quartimax') {
                
                fnames = 'cfT'
                Rot.Args <- list(Tmat=diag(m),kappa = 0, normalize=normalize, eps=1e-6, maxit=maxit)
                
                
            } else if (rotation=='CF-facparsim') {
                
                fnames = 'cfT'
                Rot.Args <- list(Tmat=diag(m),kappa = 1, normalize=normalize, eps=1e-6, maxit=maxit)
                
                
            } else if (rotation=='CF-equamax') {
                
                fnames = 'cfT'
                Rot.Args <- list(Tmat=diag(m),kappa = m/(2*p), normalize=normalize, eps=1e-6, maxit=maxit)
                
                
            } else if (rotation=='CF-parsimax') {
                
                fnames = 'cfT'
                Rot.Args <- list(Tmat=diag(m),kappa = (m-1)/(p+m-2), normalize=normalize, eps=1e-6, maxit=maxit)
                
                
            } else if (rotation=='geomin') {
                
                
                fnames = 'geominT'
                Rot.Args <- list(Tmat=diag(m),delta = geomin.delta, normalize=normalize, eps=1e-6, maxit=maxit)
                
                
            } else if (rotation=='target') {
                
                fnames = 'pstT'
                Rot.Args <- list(Tmat=diag(m), W = MWeight, Target=MTarget, normalize=normalize, eps=1e-6, maxit=maxit)
                
            } else {
                stop (paste(rotation, ' has not been implemented yet.'))
            }
            
            
        } # End of the orthogonal rotations
        
        
        list(fnames = fnames, Rot.Args = Rot.Args)
        
    } # Make.Rot.Args
    
    
    Compute.se <- function (x,R0, n, rotated, phi, dist, fm, rtype, rotation, normalize, geomin.delta, MTarget, MWeight,
    PhiWeight, PhiTarget, se, confid,Ib, FE.Arg, Rot.Controls,wxt2=1e0) {
        
        
        if (se=='information') {
            
            # information deals with only normal variables and correctly specified models
            
            if (rtype=='oblique') {
                
                analytic.se = oblq.se.augmt(Lambda = rotated, Phi = phi, Rsample=R0, N=n, extraction=fm,
                normalize=normalize, rotation=rotation, modelerror='NO', geomin.delta = geomin.delta,
                MTarget=MTarget, MWeight=MWeight, PhiWeight=PhiWeight, PhiTarget=PhiTarget,wxt2=wxt2)
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
                
            } else if (dist=='normal') {
                
                u.r = EliU(R0)
            } else {
                stop (paste(dist, ' is not recognized as a data type. Four types of data are allowed: continuous, ts, ordinal, and normal.'))
            }
            
            
            
            
            if (rtype=='oblique') {
                
                analytic.se = oblq.se.augmt(Lambda = rotated, Phi = phi, Rsample=R0, N=n, extraction=fm,
                normalize=normalize, rotation=rotation, modelerror=merror, geomin.delta = geomin.delta,
                MTarget=MTarget, MWeight=MWeight, PhiWeight=PhiWeight, PhiTarget=PhiTarget,u.r=u.r,wxt2=wxt2)
            } else { # orthogonal
                
                analytic.se = orth.se.augmt(Lambda = rotated, Rsample=R0, N=n, extraction=fm,
                normalize=normalize, rotation=rotation, modelerror=merror, geomin.delta = geomin.delta,
                MTarget=MTarget, MWeight=MWeight, u.r = u.r)
            } # orthogonal
            
            
        } else if (se == 'jackknife') {
            
            Jack.Arg <- list(bj='jackknife',Ib=n, rtype=rtype,dist=dist, Level.Confid=confid,FE.Arg=FE.Arg, fnames = Rot.Controls$fnames, Rotation.Arg=Rot.Controls$Rot.Args)
            Jack = BootJack(x,rotated,Jack.Arg)
            
        } else if (se =='bootstrap') {
            
            Boot.Arg <- list(bj='bootstrap',Ib=Ib, rtype=rtype,dist=dist, Level.Confid=confid,FE.Arg=FE.Arg, fnames = Rot.Controls$fnames, Rotation.Arg=Rot.Controls$Rot.Args)
            Boot = BootJack(x,rotated,Boot.Arg)
            
        } # bootstrap
        
        
        # output
        
        # analytic.se
        
    } # Compute.se
    
    
    #-------------------------------------------------------------------------------------------------------
    
    # check input arguments
    # external functions: get.RGamma, fa.extract
    
    
    confid = LConfid[1]
    
    if ( (is.null(x)) & (is.null(covmat)) ) stop ("Neither raw data nor the correlation matrix is provided!")
    if (! (dist == 'normal') & (is.null(x)) ) stop ("Raw data are required for non-normal distributions!")
    if ( (is.null(x)) & (is.null(n.obs)) ) stop ("The sample size is not provided for the correlation matrix!")
    if ( !(se=='bootstrap') & (CItype=='percentile') ) {
        CItype='pse'
        message ('Percentile Confidence intervals are avaible only for bootstrap: pse confidence intervals are constructed.')
    }
    ###------------------------------------------------------
    if (!(is.null(x))) {
        
        x = data.matrix(x)
        
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
    
    
    # if (A.lst$heywood==0) {
    
    if (fm=='ml') {
        statistic = (n-1) * A.lst$f
        df = ((p - factors)**2 - p - factors ) /2
        stat.stage1 = list(statistic=statistic, df = df)
        
    } else if (fm=='ols') {
        
        if (dist=='continuous') {
            u.r = AsyCovCorr(x)$asc
            
        } else if (dist=='ts') {
            u.r = TSCovCorr(x)$asc
            
        } else if (dist=='ordinal') {
            u.r = get.RGamma(x, gamma=TRUE)$GammaR
            
        } else if (dist=='normal') {
            u.r = EliU(R0)
            
        }
        stat.stage1 = Compute.stat(R0,u.r,A.lst$Unrotated)
        statistic = stat.stage1$statistic
        
        statistic = ifelse(is.nan(statistic), NaN, statistic*(n-1))
        
        
    } # ols
    
    if (stat.stage1$df <= 0) {
        
        
        message('The EFA test statistic is invalid because the degrees of freedom are not positive.')
        
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
    }
    
    
    
    
    if (factors > 1) {
        
        transformation = NULL
        if (rotation=='xtarget') {
            if(rtype=='orthogonal') {
                rtype = 'oblique'
                message('xtarget requires oblique rotation.')
            }
            rotation = 'target'
            Rot.Controls <- Make.Rot.Args(rtype,rotation,normalize,p,factors,geomin.delta,MTarget, MWeight,PhiWeight, PhiTarget,wxt2)
            Lambda.lst = do.call (Rot.Controls$fnames, append(list(A.lst$Unrotated),Rot.Controls$Rot.Args))
            transformation = (t(A.lst$Unrotated) %*% A.lst$Unrotated) %*% solve(t(Lambda.lst$loadings) %*% A.lst$Unrotated)
            
            rotation = 'xtarget'
            
        } # rotation = 'xtarget'
        
        Rot.Controls <- Make.Rot.Args(rtype,rotation,normalize,p,factors,geomin.delta,MTarget, MWeight,PhiWeight, PhiTarget,wxt2,transformation)
        Lambda.lst = do.call (Rot.Controls$fnames, append(list(A.lst$Unrotated),Rot.Controls$Rot.Args))
        
        
        
        if (useorder) {
            
            if ( ( is.null(MWeight) ) | ( is.null(MTarget) ) ) stop ("MWeight or MTarget is not specified when ordering is request.")
            
            M.in.temp = matrix(0,(p+factors),factors)
            M.in.temp [1:p,1:factors] = Lambda.lst$loadings
            
            if (rtype=='oblique') {
                M.in.temp [(p+1):(p+factors),1:factors] = Lambda.lst$Phi
            } else{
                M.in.temp [(p+1):(p+factors),1:factors] = diag(factors)
            }
            
            M.out.temp = Align.Matrix (MTarget, M.in.temp)
            
            Lambda.lst$loadings = M.out.temp[1:p,1:factors]
            if (rtype=='oblique') Lambda.lst$Phi = M.out.temp[(p+1):(p+factors),1:factors]
            
        } # (useorder)
        
        
        
        
        for (j in 1:factors) {
            if (sum(Lambda.lst$loadings[1:p,j])<0) {
                Lambda.lst$loadings[1:p,j] = Lambda.lst$loadings[1:p,j] * (-1)
                Lambda.lst$Phi[1:factors,j] = Lambda.lst$Phi[1:factors,j] * (-1)
                Lambda.lst$Phi[j,1:factors] = Lambda.lst$Phi[j,1:factors] * (-1)
            } # if
            
        } # (j in 1:factors)
        
        
        
    } # if (factors > 1)
    
    
    
    if ( ( ! (( dist=='normal') & (merror=='NO'))) & (se =='information') ) {
        se = 'sandwich'
        message('The fisher information SEs are valid only for normal data and no model error. Sandwich SEs are computed.' )
    }
    
    
    
    if (  ( (dist=='ordinal') | (dist=='ts')) & (se =='jackknife') ) {
        se = 'bootstrap'
        message('The jackknife SE estimates for ordinal data and time series data have not been developped yet; bootstrap SE estimates are computed instead.')
    }
    ###
    
    Phi = diag(factors)
    
    
    if (factors > 1) {
        if (rtype == 'oblique') Phi = Lambda.lst$Phi
        
        if (!(se =='none'))  SE = Compute.se (x, R0, n, rotated=Lambda.lst$loadings, phi=Phi, dist, fm, rtype, rotation, normalize, geomin.delta, MTarget, MWeight,
        PhiWeight, PhiTarget, se, confid, Ib, FE.Arg, Rot.Controls,wxt2)
    } else {
        if (!(se =='none')) SE = Compute.se (x, R0, n, rotated=A.lst$Unrotated, phi=Phi, dist, fm, rtype='orthogonal', rotation='unrotated', normalize, geomin.delta, MTarget, MWeight,
        PhiWeight, PhiTarget, se, confid, Ib, FE.Arg, Rot.Controls,wxt2)
    } # (factors > 1)
    
    
    
    ## Outputs
    details = list(manifest=p,factors=factors, n.obs=n, dist=dist, fm=fm, rtype=rtype, rotation=rotation, normalize=normalize,
    geomin.delta=geomin.delta, wxt2=wxt2, MTarget=MTarget, MWeight=MWeight,
    PhiWeight = PhiWeight, PhiTarget = PhiTarget, merror=merror, se=se, LConfid=LConfid, Ib=Ib)
    
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
    
    if ((factors==1)&(!(se=='none'))) {
        rotatedse = matrix(rotatedse,nrow=p,ncol=factors)
        Phise = matrix(0, factors, factors)
    }
    
    
    
    if (se=='none') {
        rotatedse=matrix(NA,p,factors)
        Phise= matrix(NA,factors,factors)
        rotatedlow=matrix(NA,p,factors)
        rotatedupper=matrix(NA,p,factors)
        Philow=matrix(NA,factors,factors)
        Phiupper=matrix(NA,factors,factors)
        
    } else {
        
        alpha = 1 - LConfid[1]
        
        if (rtype=='orthogonal') {
            Ltype = 'lambda.orth'
        } else {
            Ltype = 'lambda.oblq'
        }
        
        CI.lambda = CIs(rotated, rotatedse, alpha, type = Ltype)
        CI.Phi = CIs(Phi, Phise, alpha, type = 'Phi')
        
        rotatedlow = CI.lambda$LowerLimit
        rotatedupper = CI.lambda$UpperLimit
        
        Philow = CI.Phi$LowerLimit
        Phiupper = CI.Phi$UpperLimit
    } # (se=='none')
    
    
    Phat = unrotated %*% t(unrotated) 
    Phat = Phat - diag(diag(Phat)) + diag(p)
    Residual = R0 - Phat
    
    
    
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
    
    
    efaout = list(details = details, unrotated=unrotated,fdiscrepancy=fdiscrepancy,convergence=convergence,heywood=heywood, nq = (p*(factors+1) - factors*(factors-1)/2), compsi=A.lst$compsi,R0=R0, Phat=Phat, Residual=Residual,
    rotated=rotated,Phi=Phi, rotatedse=rotatedse, Phise= Phise, ModelF = ModelF, rotatedlow=rotatedlow, rotatedupper=rotatedupper, Philow=Philow, Phiupper=Phiupper)
    
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
    
    cat("\nSummary of Analysis: \n")
    cat("   Estimation Method:  ",x$details$fm,"\n")
    cat("   Rotation Type:  ",x$details$rtype,"\n")  
    cat("   Rotation Criterion:  ",x$details$rotation,"\n") 
    cat("   Test Statistic:  ",round(x$ModelF$f.stat,3),"\n")
    cat("   Degrees of Freedom:  ",x$ModelF$df,"\n")
    cat("   Effect numbers of Parameters:  ",x$nq,"\n")
    cat("   P value for perfect fit:  ",round(x$ModelF$p.perfect,3),"\n")
    
    
    cat("\nRotated Factor Loadings: \n")
    print(round(x$rotated,3))
    
    cat("\nFactor Correlations: \n")
    print(round(x$Phi,3))
    
    
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
    
    
    cat("\nMeasures of Fit \n")
    cat(" Root Mean Square Error of Approximation, RMSEA \n")
    cat("   Point estimate:  ",round(x$ModelF$RMSEA,3),"\n")
    cat("  ",x$details$LConfid[2]*100,"% Confidence Intervals:  (",round(x$ModelF$RMSEA.l,3), ",",round(x$ModelF$RMSEA.u,3), ") \n")
    cat("   Test Statistic:  ",round(x$ModelF$f.stat,3),"\n")
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
    
    cat("\nFactor Correlations: \n")
    print(round(x$Phi,3))
    
    cat("\nSE for rotated Factor Loadings: \n")
    print(round(x$rotatedse,3))
    
    cat("\nSE for Factor Correlations: \n")
    print(round(x$Phise,3))
    
    cat("\nLower Bounds of", x$details$LConfid[1]*100,"% CIs for Rotated Factor Loadings: \n")
    print(round(x$rotatedlow,3))
    
    cat("\nUpper Bounds of", x$details$LConfid[1]*100,"% CIs for Rotated Factor Loadings: \n")
    print(round(x$rotatedupper,3))
    
    cat("\nLower Bounds of", x$details$LConfid[1]*100,"% CIs for Factor Correlations: \n")
    print(round(x$Philow,3))
    
    cat("\nUpper Bounds of", x$details$LConfid[1]*100,"% CIs for Factor Correlations: \n")
    print(round(x$Phiupper,3))
    
    
    
} # print.summary.efa

