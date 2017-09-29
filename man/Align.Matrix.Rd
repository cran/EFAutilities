\name{Align.Matrix}

\alias{Align.Matrix}

\title{Factor Alignment}

\description{The function is to align a factor solution according to an order matrix.  The output matrix is a (p+m+1) by m matrix, where the first p rows are factor loadings of the best match, the next m rows are factor correlations of the best match, and the last row contains the sums of squared deviations for the best match and the second best match. The difference between the best match and the second best match could be considered as a confidence on the success of the aligning procedure.}

\usage{
Align.Matrix(Order.Matrix, Input.Matrix, Weight.Matrix=NULL)
}

\arguments{
  \item{Order.Matrix}{A p by m matrix: p is the number of manifest variables and m is the number of latent factors}
  \item{Input.Matrix}{A (p+m) by m matrix, the first p rows are factor loadings, the last m rows are factor correlations}
  \item{Weight.Matrix}{A p by m matrix that assigns weight to the order matrix: NULL (default)}
}


\details{
Align.Matrix is an R function to reflect and interchange columns of Input.Matrix to match those of Order.Matrix. Because it considers all possible permutations of columns of Input.Matrix,the best match in terms of the smallest sum of squared deviations between these two matrices can always be found. It may be slow if there are too many factors. 
}



\author{Guangjian Zhang}

\examples{
#Order Matrix
A <- matrix(c(0.8,0.7,0,0,0,0,0.8,0.7),nrow=4,ncol=2)

#Input.Matrix
B <-matrix(c(0,0,-0.8,-0.7,1,-0.2,0.8,0.7,0,0,-0.2,1),nrow=6,ncol=2)

Align.Matrix(Order.Matrix=A, Input.Matrix=B)

}

\keyword{ exploratory factor analysis }
\keyword{ factor rotation }
\keyword{ factor loadings }








