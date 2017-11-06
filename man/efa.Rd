\name{efa}

\alias{efa}

\title{Exploratory Factor Analysis}


\description{
Performs exploratory factor analysis under a variety of conditions. In particular, it provides standard errors for rotated factor loadings and factor correlations for normal variables, nonnormal continuous variables, and Likert scale variables with and without model error. 
}

\usage{
 efa(x=NULL, factors=NULL, covmat=NULL, n.obs=NULL, dist='normal', 
     fm='ols', rtype='oblique', rotation='CF-varimax', 
     normalize=FALSE, geomin.delta=NULL, MTarget=NULL, MWeight=NULL,
     PhiWeight = NULL, PhiTarget = NULL, useorder=FALSE, 
     se='information', LConfid=c(0.95,0.90), CItype='pse', 
     Ib=2000, mnames=NULL, fnames=NULL, merror='YES')
}

\arguments{
  \item{x}{The raw data: an n-by-p matrix where n is number of participants and p is the number of manifest variables.}
  \item{factors}{The number of factors m: specified by a researcher; the default one is the Kaiser rule which is the number of eigenvalues of covmat larger than one.}
  \item{covmat}{A p-by-p manifest variable correlation matrix.}  
  \item{n.obs}{The number of participants used in calculating the correlation matrix. This is not required when the raw data (x) is provided.}
  \item{dist}{Manifest variable distributions: 'normal'(default), 'continuous', 'ordinal' and 'ts'. 'normal' stands for normal distribution.
              'continuous' stands for nonnormal continuous distributions.  'ordinal' stands for Likert scale variable. 
               'ts' stands for distributions for time-series data.}  
  \item{fm}{Factor extraction methods: 'ols' (default) and 'ml'} 
  \item{rtype}{Factor rotation types: 'oblique' (default) and 'orthogonal'. Factors are correlated in 'oblique' rotation, and they are uncorrelated in 'orthogonal' rotation.}
  \item{rotation}{ Factor rotation criteria: 'CF-varimax' (default), 'CF-quartimax', 'target', and 'geomin'. These rotation criteria can be used in both orthogonal and oblique rotation. In addition, a fifth rotation criteria option 'xtarget' is available for oblique rotation.}
  \item{normalize}{Row standardization in factor rotation: FALSE (default) and TRUE (Kaiser standardization).} 
  \item{geomin.delta}{The controlling parameter in Geomin rotation, 0.01 as the default value.}
  \item{MTarget}{The p-by-m target matrix for the factor loading matrix in target rotation and xtarget rotation.}
  \item{MWeight}{The p-by-m weight matrix for the factor loading matrix in target rotation and xtarget rotation.}
  \item{PhiWeight}{The m-by-m target matrix for the factor correlation matrix in xtarget rotation.}
  \item{PhiTarget}{The m-by-m weight matrix for the factor correlation matrix in xtarget rotation}
  \item{useorder}{Whether an order matrix is used for factor alignment: FALSE (default) and TRUE}
  \item{se}{Methods for estimating standard errors for rotated factor loadings and factor correlations, 'information', 'sandwich', 'bootstrap', and 'jackknife'. For normal variables and ml estimation, the default method is 'information'. For all other situations, the default method is 'sandwich'. In addition, the 'bootstrap' and 'jackknife' methods can be used for raw data.}
  \item{LConfid}{Confidence levels for model parameters (factor loadings and factor correlations) and RMSEA, respectively: c(.95, .90) as default.}
  \item{CItype}{Type of confidence intervals: 'pse' (default) or 'percentile'. CIs with 'pse' are based on point and standard error estimates; CIs with 'percentile' are based on bootstrap percentiles.}
  \item{Ib}{The Number of bootstrap samples when se='bootstrap': 2000 (default)}
  \item{mnames}{Names of p manifest variables: Null (default)} 
  \item{fnames}{Names of m factors: Null (default)}
  \item{merror}{ Model error: 'YES' (default) or 'NO'. In general, we expect our model is a parisonmious representation to the complex real world. Thus, some amount of model error is unavailable. When merror = 'NO', the efa model is assumed to fit perfectly in the population.}
}


\details{
The function \code{\link{efa}} conducts exploratory factor analysis (EFA) (Gorsuch, 1983) in a variety of conditions. Data can be normal variables, non-normal continuous variables, and Likert variables. Our implementation of EFA includes three major steps: factor extraction, factor rotation, and estimating standard errors for rotated factor loadings and factor correlations.

Factors can be extracted using two methods: maximum likelihood estimation (ml) and ordinary least squares (ols). These factor loading matrices are referred to as unrotated factor loading matrices. The ml unrotated factor loading matrix is obtained using \code{\link[stats]{factanal}}. The ols unrotated factor loading matrix is obtained using \code{\link[stats]{optim}} where the residual sum of squares is minimized. The starting values for communalities are squared multiple correlations (SMCs). The test statistic and model fit measures are provided.

Four rotation criteria (CF-varimax, CF-quartimax, geomin, and target) are available for both orthogonal rotation and oblique rotation (Browne, 2001). Additionally, a new rotation criteria, xtarget, can be specified for oblique rotation. The factor rotation methods are achieved by calling functions in the package GPArotation. CF-varimax and CF-quartimax are members of the Crawford-Fugersion family (Crawford, & Ferguson, 1970) whose kappa = 1/p and kappa = 0, respectively. They are equivalent to varimax and quartimax rotation in orthogonal rotation. The equivalence does not carry over to oblique rotation, however. Although varimax and quartimax often fail to give satisfactory results in oblique rotation, CF-varimax and CF-quartimax do give satisfactory results in many oblique rotation applications. CF-quartimax rotation is equivalent to direct oblimin rotation for oblique rotation. The target matrix in target rotation can either be a fully specified matrix or a partially specified matrix. Target rotation can be considered as a procedure which is located between EFA and CFA. In CFA, if a factor loading is specified to be zero, its value is fixed to be zero; if target rotation, if a factor loading is specified to be zero, it is made to zero as close as possible. In xtarget rotation, target values can be specified on both factor loadings and factor correlations.

Confidences for rotated factor loadings and correlation matrices are constructed using point estimates and their standard error estimates. Standard errors for rotated factor loadings and factor correlations are computed using a sandwich method (Ogasawara, 1998; Yuan, Marshall, & Bentler, 2002), which generalizes the augmented information method (Jennrich, 1974). The sandwich standard error are consistent estimates even when the data distribution is non-normal and model error exists in the population. Sandwich standard error estimates require a consistent estimate of the asymptotic covariance matrix of manifest variable correlations. Such estimates are described in Browne & Shapiro (1986) for non-normal continuous variables and in Yuan & Schuster (2013) for Likert variables. Estimation of the asymptotic covariance matrix of polychoric correlations is slow if the EFA model involves a large number of Likert variables.

When manifest variables are normally distributed (dist = 'normal') and model error does not exist (merror = 'NO'), the sandwich standard errors are equivalent to the usual standard error estimates, which come from the inverse of the information matrix. The information standard error estimates in EFA is available CEFA (Browne, Cudeck, Tateneni, & Mels, 2010) and SAS Proc Factor. Mplus (Muthen & Muthen, 2015) also implemented a version of sandwich standard errors for EFA, which are robust against non-normal distribution but not model error. Sandwich standard errors computed in \code{\link{efa}} tend to be larger than those computed in Mplus. Sandwich standard errors for non-normal distributions and with model error are equivalent to the infinitesimal jackknife standard errors described in Zhang, Preacher, & Jennrich (2012). Two computationally intensive standard error methods (se='bootstrap' and se='jackknife') are also implemented. More details on standard error estimation methods in EFA are documented in Zhang (2014).

}

\value{
An object of class efa, which includes:

\item{details}{summary information about the analysis such as number of manifest variables, number of factors, sample size, factor extraction method, factor rotation method, target values for target rotation and xtarget rotation, and levels for confidence intervals.}

\item{unrotated}{the unrotated factor loading matrix}

\item{fdiscrepancy}{discrepancy function value used in factor extraction}

\item{convergence}{whether the factor extraction stage converged successfully, successful convergence indicated by 0}

\item{heywood}{the number of heywood cases}

\item{R0}{the sample correlation matrix}

\item{Phat}{the model implied correlation matrix}

\item{Residual}{the residual correlation matrix}

\item{rotated}{the rotated factor loadings}

\item{Phi}{the rotated factor correlations}

\item{rotatedse}{the standard errors for rotated factor loadings}

\item{Phise}{the standard errors for rotated factor correlations}

\item{ModelF}{the test statistic and measures of model fit}

\item{rotatedlow}{the lower bound of confidence levels for factor loadings}

\item{rotatedupper}{the upper bound of confidence levels for factor loadings}

\item{Philow}{the lower bound of confidence levels for factor correlations}

\item{Phiupper}{the lower bound of confidence levels for factor correlations}

}


\references{
Browne, M. W. (2001). An overview of analytic rotation in exploratory factor analysis. Multivariate Behavioral Research, 36, 111-150.

Browne, M. W., Cudeck, R., Tateneni, K., & Mels, G. (2010). CEFA 3.04: Comprehensive Exploratory Factor Analysis. Retrieved from http://faculty.psy.ohio-state.edu/browne/.

Browne, M. W., & Shapiro, A. (1986). The asymptotic covariance matrix of sample correlation coefficients under general conditions. Linear Algebra and its applications, 82, 169-176.

Crawford, C. B., & Ferguson, G. A. (1970). A general rotation criterion and its use in orthogonal
rotation. Psychometrika, 35 , 321-332.

Engle, R. W., Tuholsjki, S.W., Laughlin, J.E., & Conway, A. R. A. (1999). Working memory, short-term memory, and general fluid intelligence: a latent-variable approach. Journal of Experimental Psychology: General, 309-331.

Gorsuch, R. L. (1983). Factor analysis (2nd ed.). Mahwah, NJ: Lawrence Erlbaum Associates. 

Jennrich, R. I. (1974). Simplified formula for standard errors in maximum-likelihood factor analysis. British Journal of Mathematical and Statistical Psychology, 27, 122-131.

Jennrich, R. I. (2002). A simple general method for oblique rotation. Psychometrika, 67, 7-19. 

Muthen, L. K., & Muthen, B. O. (1998-2015). Mplus user's guide (7th ed.). Los Angeles, CA:
Muthen & Muthen.

Ogasawara, H. (1998). Standard errors of several indices for unrotated and rotated factors. Economic Review, Otaru University of Commerce, 49(1), 21-69.

Yuan, K., Marshall, L. L., & Bentler, P. M. (2002). A unified approach to exploratory factor analysis with missing data, nonnormal data, and in the presence of outliers. Psychometrika , 67 , 95-122. 

Yuan, K.-H., & Schuster, C. (2013). Overview of statistical estimation methods. In T. D. Little (Ed.), The Oxford handbook of quantitative methods (pp. 361-387). New York, NY: Oxford University Press.

Zhang, G. (2014). Estimating standard errors in exploratory factor analysis. Multivariate Behavioral Research, 49, 339-353.

Zhang, G., Preacher, K. J., & Jennrich, R. I. (2012). The infinitesimal jackknife with exploratory
factor analysis. Psychometrika, 77 , 634-648.


}

\author{Guangjian Zhang, Ge Jiang, Minami Hattori, and Lauren Trichtinger}


\examples{ 
#Examples using the data sets included in the packages:
  
data("CPAI537")    # Chinese personality assessment inventory (N = 537)
   
#1a) normal, ml, oblique, CF-varimax, information, merror='NO'
res1 <- efa(x=CPAI537,factors=4, fm='ml')
res1
  
#1b) confidence intervals: normal, ml, oblique, CF-varimax, information, merror='NO'
res1$rotatedlow       # lower bound for 95 percent confidence intervals for factor loadings
res1$rotatedupper     # upper bound for 95 percent confidence intervals for factor loadings
res1$Philow           # lower bound for 95 percent confidence intervals for factor correlations
res1$Phiupper         # upper bound for 95 percent confidence intervals for factor correlations	
 
#2) continuous, ml, oblique, CF-quartimax, sandwich, merror='YES'
#efa(x=CPAI537, factors=4, dist='continuous',fm='ml',rotation='CF-quartimax', merror='YES')
 
#3) continuous, ols, orthogonal, geomin, sandwich, merror='Yes'
#efa(x=CPAI537, factors=4, dist='continuous',rtype= 'orthogonal',rotation='geomin', merror='YES')
 
#4) ordinal, ols, oblique, CF-varimax, sandwich, merror='Yes'
#data("BFI228")      # Big-five inventory (N = 228)
# For ordinal data, estimating SE with the sandwich method 
#   can take time with a dataset with 44 variables
#reduced2 <- BFI228[,1:17] # extracting 17 variables corresponding to the first 2 factors
#efa(x=reduced2, factors=2, dist='ordinal', merror='YES')

#5) continuous, ml, oblique, Cf-varimax, jackknife
#efa(x=CPAI537,factors=4, dist='continuous',fm='ml', merror='YES', se= 'jackknife')
	
#6) extracting the test statistic 
#res2 <-efa(x=CPAI537,factors=4)
#res2
#res2$ModelF$f.stat
	
# #7) extended target rotation, ml
# # The data come from Engle et al. (1999) on memory and intelligence.
# datcor <- matrix(c(1.00, 0.51, 0.47, 0.35, 0.37, 0.38, 0.28, 0.34,
#                    0.51, 1.00, 0.32, 0.35, 0.35, 0.31, 0.24, 0.28,
#                    0.47, 0.32, 1.00, 0.43, 0.31, 0.31, 0.29, 0.32,
#                    0.35, 0.35, 0.43, 1.00, 0.54, 0.44, 0.19, 0.27,
#                    0.37, 0.35, 0.31, 0.54, 1.00, 0.59, 0.05, 0.19,
#                    0.38, 0.31, 0.31, 0.44, 0.59, 1.00, 0.20, 0.21,
#                    0.28, 0.24, 0.29, 0.19, 0.05, 0.20, 1.00, 0.68,
#                    0.34, 0.28, 0.32, 0.27, 0.19, 0.21, 0.68, 1.00),
#                  ncol = 8)
# 
# # Prepare target and weight matrices for lambda -------
# MTarget1 <- matrix(c(9, 0, 0,
#                      9, 0, 0,
#                      9, 0, 0, # 0 corresponds to targets
#                      0, 9, 0,
#                      0, 9, 0,
#                      0, 9, 0,
#                      0, 0, 9,
#                      0, 0, 9), ncol = 3, byrow = TRUE)
# MWeight1 <- matrix(0, ncol = 3, nrow = 8)
# MWeight1[MTarget1 == 0] <- 1 # 1 corresponds to targets
# 
# # Prepare target and weight matrices for phi ---------
# PhiTarget1 <- matrix(c(1, 9, 9,
#                        9, 1, 0,
#                        9, 0, 1), ncol = 3)
# PhiWeight1 <- matrix(0, ncol = 3, nrow = 3)
# PhiWeight1[PhiTarget1 == 0] <- 1
# 
# # Conduct extended target rotation -------------------
# mod.xtarget <- efa(covmat = datcor, factors = 3, n.obs = 133,
#                    rotation ='xtarget', fm = 'ml', useorder = T,
#                    MTarget = MTarget1, MWeight = MWeight1,
#                    PhiTarget = PhiTarget1, PhiWeight = PhiWeight1)
# mod.xtarget
# 



}

\keyword{ exploratory factor analysis }
\keyword{ factor rotation }
\keyword{ standard error }
\keyword{ factor loadings }

