\name{efaMR}
\alias{efaMR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Exploratory Factor Analysis with Multiple Rotations

}
\description{The function compares EFA solutions from  multiple random starts or from multiple rotation criteria.

}
\usage{
efaMR(x=NULL, factors=NULL, covmat=NULL, n.obs=NULL, dist='normal',
      fm='ols', rtype='oblique', rotation = 'CF-varimax', input.A=NULL,
       additionalRC = NULL,nstart = 100, compare = 'First',
      normalize=FALSE, geomin.delta=NULL, MTarget=NULL, MWeight=NULL,
       PhiTarget = NULL, PhiWeight = NULL, useorder=FALSE,
       mnames=NULL, fnames=NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The raw data: an n-by-p matrix where n is number of participants and p is the number of manifest variables.}
  \item{factors}{The number of factors m: specified by a researcher; the default one is the Kaiser rule which is the number of eigenvalues of covmat larger than one.}

  \item{covmat}{A p-by-p manifest variable correlation matrix.}  
  \item{n.obs}{The number of participants used in calculating the correlation matrix. This is not required when the raw data (x) is provided.}
  
  \item{dist}{Manifest variable distributions: 'normal'(default), 'continuous', 'ordinal' and 'ts'. 'normal' stands for normal distribution.
'continuous' stands for nonnormal continuous distributions.  'ordinal' stands for Likert scale variable. "ts" stands for distributions for time-series data.}

  \item{fm}{Factor extraction methods: 'ols' (default) and 'ml'} 
  
  \item{rtype}{Factor rotation types: 'oblique' (default) and 'orthogonal'. Factors are correlated in 'oblique' rotation, and they are uncorrelated in 'orthogonal' rotation.}
  \item{rotation}{ Factor rotation criteria: 'CF-varimax' (default), 'CF-quartimax', 'target', and 'geomin'. These rotation criteria can be used in both orthogonal and oblique rotation.}
  \item{input.A}{A p-by-m unrotated factor loading matrix. It can replace x or covmat as input arguments. Only factor rotation will be conducted; factor extraction will not be conducted.}
  \item{additionalRC}{A string of factor extraction methods against which the main rotation is compared. Required only when nstart = 1. See details. }

  \item{nstart}{The number random orthogonal starts used, with 100 as the default value.  With nstart = 1, additionalRC is required.  See details.}
  \item{compare}{'First' (default) or 'All': The global solution is compared against all local solutions with 'First'; All solutions are compared with each other with 'All'.}

  \item{normalize}{Row standardization in factor rotation: FALSE (default) and TRUE (Kaiser standardization).}   
  \item{geomin.delta}{The controlling parameter in Geomin rotation, 0.01 as the default value.}

  \item{MTarget}{The p-by-m target matrix for the factor loading matrix in target rotation or xtarget rotation.}
  \item{MWeight}{The p-by-m weight matrix for the factor loading matrix in target rotation or xtarget rotation.}
  \item{PhiTarget}{The m-by-m target matrix for the factor correlation matrix in xtarget rotation.}
  \item{PhiWeight}{The m-by-m weight matrix for the factor correlation matrix in xtarget rotation.}

  \item{useorder}{Use order matrix: FALSE(default) and TRUE}
  \item{mnames}{Names of p manifest variables: Null (default)} 
  \item{fnames}{Names of m factors: Null (default)}
  
%%     ~~Describe \code{x} here~~
}

\details{
efaMR performs EFA with multiple rotation using random starts.  

Geomin rotation, in particular, is known to produce multiple local solutions; the use of random starts is advised (Hattori, Zhang, & Preacher, in press). 

The p-by-m unrotated factor loading matrix is post-multiplied by an m-by-m random orthogonal matrices before rotation. 

The number of random starts can be specified with the default value of nstart = 100. 
Bar plot that represents frequencies of each solution is provided. 
If multiple solutions are found, they are compared with each other using congruence coefficient. 


If nstart = 1, no random start is used. The solution is compared against solutions using additional rotation criterion provided by additionalRC. 
For example, with rotation = geomin, additionalRC = c('CF-varimax', 'CF-quartimax), the geomin solution is compared against those with CF-varimax and CF-quartimax. 


Estimation of standard errors and construction of confidence intervals are disabled with multiple efaMR. 
They are available with efa. 
 
%%  ~~ If necessary, more details than the description above ~~
}

\references{
Hattori, M., Zhang, G., & Preacher, K. J. (in press).  Multiple local solutions and geomin rotation. Multivariate Behavioral Research.}
\author{Minami Hattori, Guangjian Zhang}

\examples{
#data("CPAI537")    # Chinese personality assessment inventory (N = 537)

# # Example 1: Oblique geomin rotation with 10 random starts
# res1 <- efaMR(CPAI537, factors = 5, fm = 'ml', 
#               rtype = 'oblique', rotation = 'geomin',
#               geomin.delta = .01, nstart = 10)
# res1
# summary(res1)
# res1$MultipleSolutions
# res1$Comparisons

# In practice, we recommend nstart = 100 or more (Hattori, Zhang, Preacher, in press) 
# although it may take longer time.   


# Example 2: Oblique geomin rotation (no random starts)
#            compared against CF-varimax and CF-quartimax rotation solutions
# res2 <- efaMR(CPAI537, factors = 5, fm = 'ml', 
#               rtype = 'oblique', rotation = 'geomin',
#               additionalRC = c('CF-varimax', 'CF-quartimax'), 
#               geomin.delta = .01, nstart = 1)
# res2$MultipleSolutions
# res2$Comparisons


# Example 3: Obtaining multiple solutions from the unrotated factor loading matrix as input
# res3 <- efa(CPAI537, factors = 5, fm = 'ml', 
#             rtype = 'oblique', rotation = 'geomin')
# set.seed(2017)
# res3MR <- efaMR(input.A = res3$unrotated, rtype = 'oblique',
#                 rotation = 'geomin', geomin.delta = .01)
# res3MR$MultipleSolutions
# res3MR$Comparisons
}

\keyword{ exploratory factor analysis }
\keyword{ factor rotation }
\keyword{ factor loadings }
