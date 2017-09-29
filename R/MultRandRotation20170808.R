
# Updated: August 8, 2017
# nstart: 1
# rotation: 'CF-varimax', 'CF-quartimax', 'target', 'geomin'
# rtype: 'oblique', 'orthogonal'

MultRandRotation <- function(unrotated, epsilon = .01, nstart = 100, plot = T, cex = .5,
                             rotation = 'geomin', rtype = 'oblique', normalize=F, MWeight=NULL, MTarget=NULL){
  A <- unrotated                 # unrotated factor loading
  p <- dim(A)[1]; m <- dim(A)[2] # number of variables, number of factors
  K <- nstart                    # number of random starts
  Lam <- Phi <- ali <- list()    # factor loadings, correlations, aligned loadings
  f.values <- rep(NA, K)         # geomin rotation function values

  # Rotation with random starts =======================================
  for(i in 1:K){
    RS <- qr.Q(qr(matrix(rnorm(m*m),m))) # random (m by m) orthogonal matrix
    A_start <- A %*% RS                  # random start (p by m) factor loadings

    if (rtype=='oblique') {              # oblique rotation
      if (rotation=='CF-varimax') {
        res <- cfQ(A_start, kappa = 1/p, normalize=normalize, maxit=100000)   
      } else if (rotation=='CF-quartimax') {
        res <- cfQ(A_start, kappa = 0, normalize=normalize, maxit=100000)   
      } else if (rotation=='geomin') {
        res <- geominQ(A_start, delta=epsilon, normalize=normalize, maxit = 100000)
      } else if (rotation=='target') {
        res <- pstQ(A_start, W = MWeight, Target=MTarget, normalize=normalize, maxit=100000)   
      } else {
        stop (paste(rotation, ' has not been implemented yet.'))
      }
    } else {                             # orthogonal rotations
      if (rotation=='CF-varimax') {
        res <- cfT(A_start, kappa = 1/p, normalize=normalize, maxit=100000)   
      } else if (rotation=='CF-quartimax') {
        res <- cfT(A_start, kappa = 0, normalize=normalize, maxit=100000)   
      } else if (rotation=='geomin') {
        res <- geominT(A_start, delta=epsilon, normalize=normalize, maxit = 100000)
      } else if (rotation=='target') {
        res <- pstT(A_start, W = MWeight, Target=MTarget, normalize=normalize, maxit=100000)
      } else {
        stop (paste(rotation, ' has not been implemented yet.'))
      }
    }                                    # end of rotatation
    
    Lam[[i]] <- res$loadings             # rotated factor loadings
    Phi[[i]] <- res$Phi                  # rotated factor corelations
    f.values[i] <- sum(exp(apply(log(Lam[[i]]^2 + epsilon), 1, sum) / m))
  }
  
  # Determine the number of solutions =================================
  deviation <- matrix(0, K, K)
  for(i in 1:K){
    if(i < K)
    for(j in (i+1):K){
      deviation[j,i] <- Align.Matrix(Lam[[i]], rbind(Lam[[j]],diag(m)))[(p+m+1),1]
      deviation[i,j] <- deviation[j,i]
    }
  }
  diag(deviation) <- 1
  
  solution <- rep(0, K)
  for(j in 1:K){
    first.s <- which(solution==0)[1]
    solution[first.s] <- j
    for(i in 1:K){
      if(solution[i] == 0)
      if(round(deviation[first.s,i],4) == 0) solution[i] <- j
    }
  }
  
  # Check the number of solutions also using rotation criterion values -------
  solution2 <- rep(0, K)
  for(j in 1:K){
    first.s <- which(solution2==0)[1]
    solution2[first.s] <- j # first Solution j
    for(i in 1:K){
      if(solution2[i] == 0)
      if(round(f.values[first.s] - f.values[i], 3) == 0) solution2[i] <- j
    }
  }
  if(sum(solution2 == solution) < K) stop('Something wrong in counting solutions!')
  
  # Basic descriptive outputs ==========================================
  nsol <- max(solution)
  f.solutions <- rep(NA, nsol)
  for(j in 1:nsol) f.solutions[j] <- f.values[which(solution == j)[1]]
  count.sol <- count(solution)[,2]   # here we need library(plyr)
  count.sol <- count.sol[order(f.solutions)]
  f.solutions <- round(f.solutions[order(f.solutions)],3)
  
  # Align solutions =================================================
  LAMBDAS <- PHIS <- list()
  global.solution.number <- which(round(f.values,3) == f.solutions[1])[1]
  global.lambda <- Lam[[global.solution.number]]
  global.phi <- Phi[[global.solution.number]]
  LAMBDAS[[1]] <- round(global.lambda,3)
  PHIS[[1]] <- round(global.phi,3)
  for(j in 2:nsol){
    solution.number <- which(round(f.values,3) == f.solutions[j])[1]
    lambdas <- Lam[[solution.number]]; phis <- Phi[[solution.number]]
    aligned.lambdas2 <- Align.Matrix(global.lambda, rbind(lambdas,phis))
    LAMBDAS[[j]] <- round(aligned.lambdas2[1:p,], 3)
    PHIS[[j]] <- round(aligned.lambdas2[-c(1:p,p+m+1),], 3)
  }
  
  # Make a bar plot =======================================================
  if(plot == T){  
    TAB <- as.data.frame(cbind(f.solutions, count.sol))
    if(nsol > 8){    # Case A: 9 or more solutions --------------------
      global <- TAB[1,]; locals <- TAB[-1,]        
      n.locals <- dim(locals)[1]
      seven <- locals[order(locals[,2]),][-c(1:(n.locals-7)),]
      seven.ordered <- seven[order(seven[,1]),]
      other.locals <- c(NA, sum(locals[order(locals[,2]),][c(1:(n.locals-7)),2]))
      TAB2 <- rbind(global, seven.ordered, other.locals)
      p <- barplot(TAB2[,2]/K, ylab = 'Proportions', 
                   xlab = "Rotation Criterion Values",
                   names.arg = c(as.vector(TAB2[-9,1]), 'Others'), 
                   cex.names = cex, col = c('grey50',rep(NA,7),'gray85'))
      title(main=paste('Local Solutions (',K,' random starts, ',
                       nsol,' solutions)',sep=''), cex.main = .8)
    } else {    # Case B: 8 or less solutions: No 'others' ------------
      TAB2 <- TAB; n.bars <- nsol
      p <- barplot(TAB2[,2]/K, ylab = 'Proportions', 
                   xlab = "Rotation Criterion Values",
                   names.arg = c(as.vector(TAB2[,1])), 
                   cex.names = cex, col = c('grey50',rep(NA,n.bars-1)))
      title(main = paste('Local Solutions (',K,' random starts, ', 
                         nsol,' solutions)',sep=''), cex.main = .8)
    } 
  } # if(plot == T)
  
  # Output ==========================================================
  list(RotationCriterion = f.solutions, 
       N.Solutions = nsol,
       Frequencies = count.sol,
       loadings = LAMBDAS,
       Phi = PHIS)
} # MultRandRotation() -------------------------------------------------


