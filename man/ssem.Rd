\name{ssem}
\alias{ssem}
\title{Saturated Structural Equation Models}


\description{This function fits saturated Structural equation models (SSEM) under a variety of conditions. SSEM re-parameterizes the obliquely rotated factor correlation matrix such that factors can be either endogenous or exogenous.  In comparison, all factors are exogenous in exploratory factor analysis. Manifest variables can be normal variables, nonnormal variables, nonnormal continuous variable, Likert scale variables and time series. It also provides standard errors and confidence intervals for rotated factor loadings and structural parameters.

}

\usage{
ssem(x=NULL, factors=NULL, exfactors=1, covmat=NULL,
acm=NULL, n.obs=NULL, dist='normal', fm='ml', mtest = TRUE,
rotation='semtarget', normalize=FALSE, maxit=1000, geomin.delta=NULL,
MTarget=NULL, MWeight=NULL, BGWeight = NULL, BGTarget = NULL,
PhiWeight = NULL, PhiTarget = NULL, useorder=TRUE, se='sandwich',
LConfid=c(0.95,0.90), CItype='pse', Ib=2000, mnames=NULL, fnames=NULL,
merror='YES', wxt2 = 1e0)
}

\arguments{
  \item{x}{The raw data: an n-by-p matrix where n is number of participants and p is the number of manifest variables.}
  \item{factors}{The number of factors m: specified by a researcher; the default one is the Kaiser rule which is the number of eigenvalues of covmat larger than one.}
  \item{exfactors}{The number of exogenous factors: 1 (default)}
  \item{covmat}{A p-by-p manifest variable correlation matrix.}
  \item{acm}{A p(p-1)/2 by p(p-1)/2 asymptotic covariance matrix of correlations: specified by the researcher.  }
  \item{n.obs}{The number of participants used in calculating the correlation matrix. This is not required when the raw data (x) is provided.}
  \item{dist}{Manifest variable distributions: 'normal'(default), 'continuous', 'ordinal' and 'ts'. 'normal' stands for normal distribution.
              'continuous' stands for nonnormal continuous distributions.  'ordinal' stands for Likert scale variable.
               'ts' stands for distributions for time-series data.}
  \item{fm}{Factor extraction methods: 'ml' (default) and 'ols'}
   \item{mtest}{Whether the test statistic is computed: TRUE (default) and FALSE}
  \item{rotation}{ Factor rotation criteria: 'semtarget' (default),'CF-varimax', 'CF-quartimax', 'CF-equamax', 'CF-parsimax', 'CF-facparsim','target', and 'geomin'. These rotation criteria can be used in both orthogonal and oblique rotation. In addition, a fifth rotation criterion 'xtarget'(extended target) rotation is available for oblique rotation. The ssem target rotation allows targets to be specified on both factor loadings and factor structural parameters.}
  \item{normalize}{Row standardization in factor rotation: FALSE (default) and TRUE (Kaiser standardization).}
  \item{maxit}{Maximum number of iterations in factor rotation: 1000 (default)}
  \item{geomin.delta}{The controlling parameter in Geomin rotation, 0.01 as the default value.}
  \item{MTarget}{The p-by-m target matrix for the factor loading matrix in target rotation and semtarget rotation.}
  \item{MWeight}{The p-by-m weight matrix for the factor loading matrix in target rotation and semtarget rotation. Optional}
  \item{BGWeight}{The m1-by-m weight matrix for the [Beta | Gamma] matrix in semtarget rotation (see details) Optional}
  \item{BGTarget}{The m1-by-m target matrix for the [Beta | Gamma] matrix in semtarget rotation where m1 is the number of endogenous factors (see details)}
  \item{PhiWeight}{The m2-by-m2 target matrix for the exogenous factor correlation matrix in semtarget rotation.Optional}
  \item{PhiTarget}{The m2-by-m2 weight matrix for the exogenous factor correlation matrix in semtarget rotation}
  \item{useorder}{Whether an order matrix is used for factor alignment: TRUE (default) and FALSE}
  \item{se}{Methods for estimating standard errors for rotated factor loadings and factor correlations, 'sandwich' (default),'information',  'bootstrap', and 'jackknife'. The 'bootstrap' and 'jackknife' methods require raw data.}
  \item{LConfid}{Confidence levels for model parameters (rotated factor loadings and structural parameters) and RMSEA, respectively: c(.95, .90) as default.}
  \item{CItype}{Type of confidence intervals: 'pse' (default) or 'percentile'. CIs with 'pse' are based on point and standard error estimates; CIs with 'percentile' are based on bootstrap percentiles.}
  \item{Ib}{The Number of bootstrap samples when se='bootstrap': 2000 (default)}
  \item{mnames}{Names of p manifest variables: Null (default)}
  \item{fnames}{Names of m factors: Null (default)}
  \item{merror}{ Model error: 'YES' (default) or 'NO'. In general, we expect our model is a parsimonious representation to the complex real world. Thus, some amount of model error is unavailable. When merror = 'NO', the ssem model is assumed to fit perfectly in the population.}
  \item{wxt2}{The relative weight for structural parameters in 'semtarget' rotation: 1 (default)}
}


\details{
The function \code{\link{ssem}} conducts saturated structural equation modeling (ssem) in a variety of conditions. Data can be normal variables, non-normal continuous variables, and Likert variables. Our implementation of SSEM includes three major steps: factor extraction, factor rotation, and estimating standard errors for rotated factor loadings and factor correlations.

Factors can be extracted using two methods: maximum likelihood estimation (ml) and ordinary least squares (ols). These factor loading matrices are referred to as unrotated factor loading matrices. The ml unrotated factor loading matrix is obtained using \code{\link[stats]{factanal}}. The ols unrotated factor loading matrix is obtained using \code{\link[stats]{optim}} where the residual sum of squares is minimized. The starting values for communalities are squared multiple correlations (SMCs). The test statistic and model fit measures are provided.

Eight rotation criteria (semtarget, CF-varimax, CF-quartimax, CF-equamax, CF-parsimax, CF-facparsim, target, and geomin) are available for oblique rotation (Browne, 2001). Additionally, a new rotation criteria, ssemtarget, can be specified for oblique rotation. The factor rotation methods are achieved by calling functions in the package GPArotation. CF-varimax, CF-quartimax, CF-equamax, CF-parsimax, and CF-facparsim are members of the Crawford-Fugersion family (Crawford, & Ferguson, 1970) whose kappa = 1/p and kappa = 0, respectively.  The target matrix in target rotation can either be a fully specified matrix or a partially specified matrix. Target rotation can be considered as a procedure which is located between EFA and CFA. In CFA, if a factor loading is specified to be zero, its value is fixed to be zero; if target rotation, if a factor loading is specified to be zero, it is made to zero as close as possible. In xtarget rotation, target values can be specified on both factor loadings and factor correlations. In ssemtarget, target values can be specified for the [Beta | Gamma] matrix where Beta is the regression weights of the endogenous factors on itself and the Gamma is the regression weights of the endogenous factors on the exogenous factors.

Confidence intervals for rotated factor loadings and correlation matrices are constructed using point estimates and their standard error estimates. Standard errors for rotated factor loadings and factor correlations are computed using a sandwich method (Ogasawara, 1998; Yuan, Marshall, & Bentler, 2002), which generalizes the augmented information method (Jennrich, 1974). The sandwich standard error are consistent estimates even when the data distribution is non-normal and model error exists in the population. Sandwich standard error estimates require a consistent estimate of the asymptotic covariance matrix of manifest variable correlations. Such estimates are described in Browne & Shapiro (1986) for non-normal continuous variables and in Yuan & Schuster (2013) for Likert variables. Estimation of the asymptotic covariance matrix of polychoric correlations is slow if the EFA model involves a large number of Likert variables.

When manifest variables are normally distributed (dist = 'normal') and model error does not exist (merror = 'NO'), the sandwich standard errors are equivalent to the usual standard error estimates, which come from the inverse of the information matrix. The information standard error estimates in EFA is available CEFA (Browne, Cudeck, Tateneni, & Mels, 2010) and SAS Proc Factor. Mplus (Muthen & Muthen, 2015) also implemented a version of sandwich standard errors for EFA, which are robust against non-normal distribution but not model error. Sandwich standard errors computed in \code{\link{efa}} tend to be larger than those computed in Mplus. Sandwich standard errors for non-normal distributions and with model error are equivalent to the infinitesimal jackknife standard errors described in Zhang, Preacher, & Jennrich (2012). Two computationally intensive standard error methods (se='bootstrap' and se='jackknife') are also implemented. More details on standard error estimation methods in EFA are documented in Zhang (2014).

}

\value{
An object of class ssem, which includes:

\item{details}{summary information about the analysis such as number of manifest variables, number of factors, number of endogenous factors, number of exogenous factors, sample size, distribution, factor extraction method, factor rotation method, target values for target rotation, xtarget rotation and ssemtarget rotation, and levels for confidence intervals.}

\item{unrotated}{the unrotated factor loading matrix}

\item{fdiscrepancy}{discrepancy function value used in factor extraction}

\item{convergence}{whether the factor extraction stage converged successfully, successful convergence indicated by 0}

\item{heywood}{the number of heywood cases}

\item{nq}{the number of effective parameters}

\item{compsi}{contains eigenvalues, SMCs, communalities, and unique variances}

\item{R0}{the sample correlation matrix}

\item{Phat}{the model implied correlation matrix}

\item{Residual}{the residual correlation matrix}

\item{rotated}{the rotated factor loadings}

\item{Phi}{the rotated factor correlations}

\item{BG}{the [Beta | Gamma] latent regression coefficients}

\item{psi}{the endogenous residuals}

\item{Phi.xi}{the exogenous correlation}

\item{rotatedse}{the standard errors for rotated factor loadings}

\item{Phise}{the standard errors for rotated factor correlations}

\item{BGse}{the standard errors for the [Beta | Gamma] latent regression coefficients}

\item{psise}{the standard errors for the endogenous residuals}

\item{Phi.xise}{the standard errors for the exogenous correlation}

\item{ModelF}{the test statistic and measures of model fit}

\item{rotatedlow}{the lower bound of confidence levels for factor loadings}

\item{rotatedupper}{the upper bound of confidence levels for factor loadings}

\item{Philow}{the lower bound of confidence levels for factor correlations}

\item{Phiupper}{the lower bound of confidence levels for factor correlations}

\item{BGlower}{the lower bound of the [Beta | Gamma] latent regression coefficients}

\item{BGupper}{the upper bound of the [Beta | Gamma] latent regression coefficients}

\item{psilower}{the lower bound of the endogenous residuals}

\item{psiupper}{the upper bound of the endogenous residuals}

\item{Phixilower}{the lower bound of the exogenous correlation}

\item{Phixiupper}{the upper bound of the exogenous correlation}

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

\author{Guangjian Zhang, Minami Hattori, and Lauren Trichtinger}


\examples{
#cormat <- matrix(c(1, .865, .733, .511, .412, .647, -.462, -.533, -.544,
#                   .865, 1, .741, .485, .366, .595, -.406, -.474, -.505,
#                   .733, .741, 1, .316, .268, .497, -.303, -.372, -.44,
#                   .511, .485, .316, 1, .721, .731, -.521, -.531, -.621,
#                   .412, .366, .268, .721, 1, .599, -.455, -.425, -.455,
#                   .647, .595, .497, .731, .599, 1, -.417, -.47, -.521,
#                  -.462, -.406, -.303, -.521, -.455, -.417, 1, .747, .727,
#                   -.533, -.474, -.372, -.531, -.425, -.47, .747, 1, .772,
#                   -.544, -.505, -.44, -.621, -.455, -.521, .727, .772, 1),
#                 ncol = 9)


#p <- 9      # a number of manifest variables

#m <- 3      # a total number of factors

#m1 <- 2     # a number of endogenous variables
#N <- 138    # a sample size

#mvnames <- c("H1_likelihood", "H2_certainty", "H3_amount", "S1_sympathy",
#             "S2_pity", "S3_concern", "C1_controllable", "C2_responsible", "C3_fault")

#fnames <- c('H', 'S', 'C')
# Step 2: Preparing target and weight matrices =========================
# a 9 x 3 matrix for lambda; p = 9, m = 3

#MT <- matrix(0, p, m, dimnames = list(mvnames, fnames))

#MT[c(1:3,6),1] <- 9

#MT[4:6,2] <- 9

#MT[7:9,3] <- 9

#MW <- matrix(0, p, m, dimnames = list(mvnames, fnames))

#MW[MT == 0] <- 1

# a 2 x 3 matrix for [B|G]; m1 = 2, m = 3

# m1 = 2
#BGT <- matrix(0, m1, m, dimnames = list(fnames[1:m1], fnames))

#BGT[1,2] <- 9

#BGT[2,3] <- 9

#BGT[1,3] <- 9

#BGW <- matrix(0, m1, m, dimnames = list(fnames[1:m1], fnames))

#BGW[BGT == 0] <- 1

#BGW[,1] <- 0

#BGW[2,2] <- 0
# a 1 x 1 matrix for Phi.xi; m - m1 = 1 (only one exogenous factor)

#PhiT <- matrix(9, m - m1, m - m1)

#PhiW <- matrix(0, m - m1, m - m1)
#SSEMres <- ssem(covmat = cormat, factors = m, exfactors = m - m1,
#                dist = 'normal', n.obs = N, fm = 'ml', rotation = 'semtarget',
#                maxit = 10000,
#                MTarget = MT, MWeight = MW, BGTarget = BGT, BGWeight = BGW,
#                PhiTarget = PhiT, PhiWeight = PhiW,  useorder = TRUE, se = 'information',
#                mnames = mvnames, fnames = fnames)
#}

\keyword{saturated structural equation model}



