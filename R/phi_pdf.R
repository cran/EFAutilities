phi_pdf <-
function(s, t, rho, tau)
{
  n1 = length(tau[[1]])+1; n2 = length(tau[[2]])+1; 
  if (s==0|s==n1|t==0|t==n2) p = 0 else {
    p = dmvnorm(c(tau[[1]][s], tau[[2]][t]), mean = c(0,0), sigma = matrix(c(1, rho, rho, 1), 2))}
  return(p)		
}

