pi_st <-
function(s, t, rho, tau)
{	
  tau[[1]] <- c(-Inf, tau[[1]], Inf); tau[[2]] <- c(-Inf, tau[[2]], Inf)
  p = pmvnorm(lower = c(tau[[1]][s],tau[[2]][t]), upper = c(tau[[1]][s+1],tau[[2]][t+1]), mean = c(0,0), sigma = matrix(c(1, rho, rho, 1), 2))    
  return(p)		
}
