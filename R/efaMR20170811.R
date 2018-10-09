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

# efaMR: Exploratory factor analysis with multiple rotations
# new functions: CompareSolutions and MultRandRotation

# Last Updated: 8pm, August 9, 2017


# ===================================================================================
efaMR <- function(x=NULL, factors=NULL, covmat=NULL, n.obs=NULL, 
                  dist='normal', fm='ols', rtype='oblique', rotation = 'CF-varimax', 
                  input.A=NULL, additionalRC = NULL, 
                  nstart = 100, compare = 'First', 
                  normalize=FALSE, geomin.delta=NULL, 
                  MTarget=NULL, MWeight=NULL, PhiTarget = NULL, PhiWeight = NULL, 
                  useorder=FALSE, mnames=NULL, fnames=NULL) {
  
# 1. check input arguments -----------------------------------------------------------
  if ( (is.null(x)) & (is.null(covmat)) & (is.null(input.A)) ) stop ("Neither raw data, the correlation matrix, nor the unrotated factor loading matrix is provided!")
  if (is.null(input.A)){
    if (! (dist == 'normal') & (is.null(x)) ) stop ("Raw data are required for non-normal distributions!")
    if ( (is.null(x)) & (is.null(n.obs)) ) stop ("The sample size is not provided for the correlation matrix!")
  } else message ("Factor rotation is conducted with the input unrotated factor loading matrix; factor extraction is unnecessary.")
  
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
      if(min(abs(R0 - covmat)) > 0.0001) message ('covmat is different from the one computed from the raw data! The one computed from raw data is used for EFA!') 
    } # if (!(is.null(covmat)))
  } else {            # (!(is.null(x)))

    if(is.null(input.A)){
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
    }# if(is.null(input.A))
  } # is.null(x) 
  
  if(!is.null(input.A)){
    p <- dim(input.A)[1]
    factors <- dim(input.A)[2]
    n <- n.obs
  }
  
  if(is.null(input.A)){
      max.factors = floor(((2*p + 1) - sqrt(8*p+1))/2)
      
      ev = eigen(R0)$values
      
      if (is.null(factors)) {factors = length(which(ev > 1))
      } else {
        if (factors > max.factors) stop (paste(factors, "factors is too many for",p,"variables."))
      }
  }
# 2. manifest variable names and factor names -------------------------------------
  if (is.null(mnames)) {
    mnames = rep(" ", p)  
    for(i in 1:p) {
      mnames[i] = paste("MV",i,sep="")
    }
  } else{
    if(length(mnames) != p) stop("The number of MV names is different from the number of MVs!")
  }
  
  if (is.null(fnames)) {
    fnames = rep(" ", factors)  
    for(i in 1:factors) {
      fnames[i] = paste("F",i,sep="")
    }
  } else{
    if(length(fnames) != factors) stop("The number of factor names is different from the number of factors!")
  }
  
# 3. extract m factors -------------------------------------------------------
  if(is.null(input.A)){
    A.lst = fa.extract(R0, factors, extraction = fm)
    FE.Arg <- list(factors = factors, extraction = fm)
    unrotated <- A.lst$Unrotated
  } else { 
  # if unrotated factor loading matrix is provided instead of (in addition to) x
    unrotated <- input.A
  }

  
# 4. factor rotation with multiple random starts ------------------------------

# 4A. nstart!=1 --- mutiple random starts for one rotation criterion ----------
  if(nstart!=1){  
      if(factors == 1){
        MultipleSolutions <- 'No rotation has been done for one-factor models.'
        Comparisons <- NULL
      }  else  { # if(factors!=1)

        
        
        multiple <- MultRandRotation(unrotated, epsilon = geomin.delta, nstart = nstart, 
                                     plot = T, cex = .5, rotation = rotation, rtype = rtype, 
                                     normalize = normalize, MWeight = MWeight, MTarget = MTarget)

        NumberSolutions <- length(multiple$loadings)
        Solutions <- list()
        for(i in 1:NumberSolutions){
          Solutions[[i]] <- list(Lambda = multiple$loadings[[i]], Phi = multiple$Phi[[i]])
        }

  
        MultipleSolutions <- list(nstart = nstart,
                                  RotationCriterion    = multiple$RotationCriterion, 
                                  NumberSolutions      = multiple$N.Solutions,
                                  FrequenciesSolutions = multiple$Frequencies,
                                  Solutions            = Solutions)
      
        if(multiple$N.Solutions != 1){
          comparison <- CompareSolutions(multiple$loadings, compare)
          Comparisons <- list(MinimumCongruence = comparison$MinimumCongruence,
                              RawCongruences    = comparison$RawCongruences)
        } else { 
          Comparisons <- 'No comparison with only one solution!'  
        }
      } # if(factors!=1)
    
    } else { # if(nstart==1)

# 4B. nstart==1 --- no random start for multiple rotation criteria --------
      if(factors == 1){
        MultipleSolutions <- 'No rotation has been done for one-factor models.'
        Comparisons <- NULL
      }  else  { 
      # if(factors!=1) ---------------------------------------
        if(is.null(additionalRC)) stop('We need additional rotation criteria against which we compare!!')
        allRC <- c(rotation, additionalRC)
        nsolutions <- length(allRC)
        Lam3 <- Phis3 <- list()
        
        for(i in 1:nsolutions){ # rotation --------------------------------
          if (rtype=='oblique') {              # oblique rotation --------
            if (allRC[i]=='CF-varimax') {
              res3 <- cfQ(unrotated, kappa = 1/p, normalize = normalize, maxit=100000)   
            } else if (allRC[i]=='CF-quartimax') {
              res3 <- cfQ(unrotated, kappa = 0, normalize=normalize, maxit=100000)   
            } else if (allRC[i]=='geomin') {
              res3 <- geominQ(unrotated, delta = geomin.delta, normalize =normalize, maxit = 100000)
            } else if (allRC[i]=='target') {
              res3 <- pstQ(unrotated, W = MWeight, Target=MTarget, normalize=normalize, maxit=100000)   
            } else if (allRC[i]=='xtarget'){
              res3 <- xpstQ(unrotated, methodArgs = list(W = MWeight, Target = MTarget), 
                           PhiWeight = PhiWeight, PhiTarget = PhiTarget,
                           normalize=normalize, maxit=100000)
            } else {
              stop (paste(allRC[i], ' has not been implemented yet.'))
            }
          } else {                  # orthogonal rotations ---------
            if (allRC[i]=='CF-varimax') {
              res3 <- cfT(unrotated, kappa = 1/p, normalize=normalize, maxit=100000)   
            } else if (allRC[i]=='CF-quartimax') {
              res3 <- cfT(unrotated, kappa = 0, normalize=normalize, maxit=100000)   
            } else if (allRC[i]=='geomin') {
              res3 <- geominT(unrotated, delta=geomin.delta, maxit = 100000)
            } else if (allRC[i]=='target') {
              res3 <- pstT(unrotated, W = MWeight, Target=MTarget, normalize=normalize, maxit=100000)
            } else {
              stop (paste(allRC[i], ' has not been implemented yet.'))
            }
          } # end of rotatation ------------------------------------
          
          Lam3[[i]] <- round(res3$loadings,3)   # rotated factor loadings
          
          if(rtype=='oblique'){
            Phis3[[i]] <- round(res3$Phi,3)      # rotated factor corelations 
          } else{  
            Phis3[[i]] <- diag(factors) # identity phi with orthogonal rotation
          }
          
          # factor alignment -----------------------
          if(i != 1){
            aligned.lambdas2 <- Align.Matrix(Lam3[[1]], rbind(Lam3[[i]], Phis3[[i]]))
            Lam3[[i]] <- aligned.lambdas2[1:p,]
            Phis3[[i]] <- aligned.lambdas2[-c(1:p,p+factors+1),]
          } # the end of factor alignment ----------

        } #for(i in 1:nsolutions) -----------------------------------

        NumberSolutions <- length(Lam3)
        Solutions <- list()
        for(i in 1:NumberSolutions){
          Solutions[[i]] <- list(Lambda = Lam3[[i]], Phi = Phis3[[i]])
        }
        
        multiple <- list(loadings = Lam3, Phi = Phis3)
        MultipleSolutions <- list(nstart             = nstart,
                                  RotationCriterion  = allRC, 
                                  Solutions          = Solutions)
        

        comparison <- CompareSolutions(multiple$loadings, compare)
        Comparisons <- list(MinimumCongruence = comparison$MinimumCongruence,
                            RawCongruences    = comparison$RawCongruences)
      } # if(factors!=1)
  } # if(nstart==1) ----------------------------------------------------
  

# 5. organizing details -------------------------------------------------------
  details = list(manifest = p, factors = factors, n.obs = n, dist = dist, 
                 fm = fm, rtype = rtype, rotation = rotation, 
                 normalize = normalize, geomin.delta = geomin.delta, 
                 MTarget = MTarget, MWeight = MWeight)
  
# 6. organizing outputs ------------------------------------------------------- 
  if(is.null(input.A)){
    fdiscrepancy = A.lst$f
    convergence = A.lst$convergence
  } else {
    fdiscrepancy <- NULL
    convergence <- NULL
  }
  
  dimnames(unrotated) = list(mnames,fnames)

# 7. output ----------------------------------------------------------------  
  efaout = list(details = details, unrotated = round(unrotated,3), 
                fdiscrepancy = fdiscrepancy, convergence = convergence,
                MultipleSolutions = MultipleSolutions, Comparisons = Comparisons)
  
  class(efaout) <- "efaMRS" 
  
  return(efaout)
} # efaMRS ==================================================================
# ==========================================================================



# ==========================================================================
# print.efaMRS <- function(x, ...) {
#   cat("\nSummary of Analysis: \n")
#   cat("   Estimation Method:  ",x$details$fm,"\n")
#   cat("   Rotation Type:  ",x$details$rtype,"\n")
#   cat("   Rotation Criterion:  ",x$details$rotation,"\n")
#
#   if(x$MultipleSolutions$nstart!=1){
#     if(x$details$factors == 1){
#       cat("\nNo rotation has been done for one-factor models \n")
#     } else{
#       cat("\nMultiple Solutions: \n")
#       cat('   There were',x$MultipleSolutions$NumberSolutions, 'solutions with',x$MultipleSolutions$nstart,'random starts.\n')
#       cat('   Rotation criteria values are',x$MultipleSolutions$RotationCriterion,'\n')
#       cat('   Frequencies for each solution are',x$MultipleSolutions$FrequenciesSolutions,'\n')
#       cat('      *** More details are available through summary(object) and all information can be retrieved using \'$\' ...\n')
#     }
#   }
# } # print.efaMRS
#
# # ======================================================================
# summary.efaMRS <- function(object, ...){
#   res <- object
#   class(res) <- "summary.efaMRS"
#   return(res)
# } # summary.efaMRS
#
# # ========================================================================
# print.summary.efaMRS <- function(x, ...) {
#   cat("\nAnalysis Details: \n")
#   cat("   Number of Manifest Variable:  ", x$details$manifest, "\n")
#   cat("   Number of Factors:  ",x$details$factors,"\n")
#   cat("   Sample Size:  ",x$details$n.obs,"\n")
#   cat("   Manifest Variable Distribution:  ",x$details$dist,"\n")
#   cat("   Estimation Method:  ",x$details$fm,"\n")
#   cat("   Rotation Type:  ",x$details$rtype,"\n")
#   cat("   Rotation Criterion:  ",x$details$rotation,"\n")
#   cat("   Rotation Standardization:  ",x$details$normalize,"\n")
#
#   cat("\nUnrotated Factor Loadings: \n")
#   print(round(x$unrotated,3))
#
#   if(x$details$factors !=1){
#     print(x$MultipleSolutions)
#     print(x$Comparisons)
#   }
# } # print.summary.efaMRS
