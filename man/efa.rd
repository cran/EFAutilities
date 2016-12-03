\name{efa}

\alias{efa}

\title{Exploratory Factor Analysis
}


\description{
Performs exploratory factor analysis under a variety of conditions. In particular, it provides standard errors for rotated factor loadings and factor correlations for normal variables, nonnormal continuous variables, and Likert scale variables with and without model error. 
}

\usage{
efa(x = NULL, factors = NULL, covmat = NULL, n.obs = NULL,
dist = "normal", fm = "ols", rtype = "oblique",
rotation = "CF-varimax", normalize = FALSE,
geomin.delta = NULL, MTarget = NULL, MWeight = NULL,
PhiWeight = NULL, PhiTarget = NULL, useorder = FALSE,
se = "information", Ib = 2000, 
mnames = NULL, fnames = NULL, merror = "YES")
}

\arguments{
  \item{x}{The raw data: an n-by-p matrix where n is number of participants and p is the number of manifest variables.}

  \item{covmat}{A p-by-p manifest variable correlation matrix.}
  
 \item{n.obs}{The number of participants used in calculating the correlation matrix. This is not required when the raw data (x) is provided.}

  \item{factors}{The number of factors m: specified by a researcher; the default one is the Kaiser rule which is the number of eigenvalues of covmat larger than one.}
  
  \item{fm}{Factor extraction methods: 'ols' (default) and 'ml'} 
  
  \item{dist}{Manifest variable distributions: 'normal'(default), 'continuous', 'ordinal' and 'ts'. 'normal' stands for normal distribution.
'continuous' stands for nonnormal continuous distributions.  'ordinal' stands for Likert scale variable. "ts" stands for distributions for time-series data.}
  
  \item{merror}{ Model error: 'YES' (default) or 'NO'. In general, we expect our model is a parisonmious representation to the complex real world. Thus, some amount of model error is unavailable. When merror = 'NO', the efa model is assumed to fit perfectly in the population.}

  \item{rtype}{Factor rotation types: 'oblique' (default) and 'orthogonal'. Factors are correlated in 'oblique' rotation, and they are uncorrelated in 'orthogonal' rotation.}

  \item{rotation}{ Factor rotation criteria: 'CF-varimax' (default), 'CF-quartimax', 'target', and 'geomin'. These rotation criteria can be used in both orthogonal and oblique rotation.}

  \item{normalize}{Row standardization in factor rotation: FALSE (default) and TRUE (Kaiser standardization).}
  
 \item{se}{Methods for estimating standard errors for rotated factor loadings and factor correlations, 'information', 'sandwich', 'bootstrap', and 'jackknife'. For normal variables and ml estimation, the default method is 'information'. For all other situations, the default method is 'sandwich'. In addition, the 'bootstrap' and 'jackknife' methods can be used for raw data. }
 
  \item{geomin.delta}{The controlling parameter in Geomin rotation, 0.01 as the default value.}

  \item{MTarget}{The p-by-m target matrix for the factor loading matrix in target rotation.}

  \item{MWeight}{The p-by-m weight matrix for the factor loading matrix in target rotation.}

  \item{PhiWeight}{The m-by-m target matrix for the factor correlation matrix in xtarget rotation}
  \item{PhiTarget}{The m-by-m weight matrix for the factor correlation matrix in xtarget rotation}
  \item{useorder}{Whether an order matrix is used for factor alignment: FALSE (default) and TRUE}
  \item{Ib}{The Number of bootstrap samples when se='bootstrap': 2000 (default)}
  \item{mnames}{Names of p manifest variables: Null (default)} 
  \item{fnames}{Names of m factors: Null (default)}
}


\details{
% Introduce the factor model. 

% Describe possible manifest variables data distirbuions.

The function \code{\link{efa}} conducts exploratory factor analysis (EFA) (Gorsuch, 1983) in a variety of conditions. Data can be normal variables, non-normal continuous variables, and Likert variables. Our implementation of EFA includes three major steps: factor extraction, factor rotation, and estimating standard errors for rotated factor loadings and factor correlations.

Factors can be extracted using two methods: maximum likelihood estimation (ml) and ordinary least squares (ols). These factor loading matrices are referred to as unrotated factor loading matrices. The ml unrotated factor loading matrix is obtained using \code{\link[stats]{factanal}}. The ols unrotated factor loading matrix is obtained using \code{\link[stats]{optim}} where the residual sum of squares is minimized. The starting values for communalities are squared mutiple correlations (SMCs).

Four rotation criteria (CF-varimax, CF-quartimax, geomin, and target) are available for both orthogonal rotation and oblique rotation (Browne, 2001). The factor rotation methods are achieved by calling functions in the package GPArotation. CF-varimax and CF-quartimax are members of the Crawford-Fugersion family (Crawford, & Ferguson, 1970) whose kappa = 1/p and kappa = 0, respectively. They are equivalent to varimax and quartimax rotation in orthogonal rotation. The equivalence does not carry over to oblique rotation, however. Althogh variamx and quartimax often fail to give satisfactory results in oblique rotation, CF-varimax and CF-quartimax do give satisfactory results in many oblique rotation applications. CF-quartimax rotation is equivalent to direct oblimin rotation for oblique rotation. The target matrix in target rotation can either be a fully specified matrix or a partially specified matrix. Target rotation can be considered as a procedure which is located between EFA and CFA. In CFA, if a factor loading is specified to be zero, its value is fixed to be zero; if target rotation, if a factor loading is specified to be zeor, it is made to zero as close as possible.



Standard errors for rotated factor loadings and factor correlations are computed using a sandwich method (Ogasawara, 1998; Yuan, Marshall, & Bentler, 2002), which generalizes the augmented information method (Jennrich, 1974). The sandwich standard error are consistent estimates even when the data distribution is non-normal and model error exists in the population. Sandwich standard error estimates require a consistent estimate of the asymptotic covariance matrix of manifest variable correlations. Such estimates are described in Browne & Shapiro (1986) for non-normal continuous variables and in Yuan & Schuster (2013) for Likert variables. Estimation of the asymptotic covariance matrix of polychoric correlations is slow if the EFA model involves a large number of Likert variables.


When manifest variables are normally distributed (dist = 'normal') and model error does not exist (merror = 'NO'), the sandwich standard errors are equivalent to the usual standard error estimates, which come from the inverse of the information matrix. The information standard error estimates in EFA is available CEFA (Browne, Cudeck, Tateneni, & Mels, 2010) and SAS Proc Factor. Mplus (Muthen & Muthen, 2015) also implemented a version of sandwich standard errors for EFA, which are robust against non-normal distrubiton but not model error. Sandwich standard errors computed in \code{\link{efa}} tend to be larger than those computed in Mplus. Sandwich standard errors for non-normal distributions and with model error are equivalent to the infinitesimal jackknife standard errors described in Zhang, Preacher, & Jennrich (2012). Two computationally intensive standard error methods (se='bootstrap' and se='jackknife') are also implmented. More details on standard error estimation methods in EFA are documented in Zhang (2014).


}


\references{
Browne, M. W. (2001). An overview of analytic rotation in exploratory factor analysis. Multivariate Behavioral Research, 36, 111-150.

Browne, M. W., Cudeck, R., Tateneni, K., & Mels, G. (2010). CEFA 3.04: Comprehensive Ex- ploratory Factor Analysis. Retrieved from http://faculty.psy.ohio-state.edu/browne/.

Browne, M. W., & Shapiro, A. (1986). The asymptotic covariance matrix of sample correlation coefficients under general conditions. Linear Algebra and its applications, 82, 169-176.

Crawford, C. B., & Ferguson, G. A. (1970). A general rotation criterion and its use in orthogonal
rotation. Psychometrika, 35 , 321-332.

Gorsuch, R. L. (1983). Factor analysis (2nd ed.). Mahwah, NJ: Lawrence Erlbaum Associates. 

Jennrich, R. I. (1974). Simplified formulae for standard errors in maximum-likelihood factor analysis. British Journal of Mathematical and Statistical Psychology, 27, 122-131.

Jennrich, R. I. (2002). A simple general method for oblique rotation. Psychometrika, 67, 7-19. 

Muthen, L. K., & Muthen, B. O. (1998-2015). Mplus user's guide (7th ed.). Los Angeles, CA:
Muthen & Muthen.

Ogasawara, H. (1998). Standard errors of several indices for unrotated and rotated factors. Economic Review, Otaru University of Commerce, 49(1), 21-69.

Yuan, K., Marshall, L. L., & Bentler, P. M. (2002). A unified approach to exploratory factor analysis with missing data, nonnormal data, and in the presence of outliers. Psychometrika , 67 , 95-122. 

Yuan, K.-H., & Schuster, C. (2013). Overview of statistical estimation methods. In T. D. Little (Ed.), The Oxford handbook of quantitative methods (pp. 361-387). New York, NY: Oxford University Press.

Zhang, G. (2014). Estimating standard errors in exploratory factor analysis. Multivariate Behavioral Research, 49, 339-353.

Zhang, G., Preacher, K. J., & Jennrich, R. I. (2012). The innitesimal jackknife with exploratory
factor analysis. Psychometrika, 77 , 634-648.

%% ~put references to the literature/web site here ~
}

\author{Guangjian Zhang, Ge Jiang, Minami Hattori, and Lauren Trichtinger
%%  ~~who you are~~
}



\examples{
 
  #Examples using the data sets included in the packages:
  
  data("CPAI537")    # Chinese personality assessment inventory (N = 537)
   
  #1) normal, ml, oblique, CF-varimax, information, merror='NO'
  efa(x=CPAI537,factors=4, fm='ml')
 
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

% # discuss with Guangjian

}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ exploratory factor analysis }
\keyword{ factor rotation }
\keyword{ standard error }
\keyword{ factor loadings }


