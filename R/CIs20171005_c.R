# This function constructs confidence intervals using point and SE estimates. 
# It accepts factor loadings, factor correlations, communalities, and unique variances

# type: 'lambda', 'Phi', 'h'(communalities), 'uv'(unique variances)

######### lines to include in efa.R
# alpha <- 1 - LConfid[1]
# type <- 'lambda'
# point <- rotated
# se <- rotatedse

# type <- 'Phi'
# point <- Phi
# se <- Phise


CIs <- function(point, se, alpha, type = type){
  
  crit <- qnorm(1- alpha/2 ) 
  
  if(any(se<0)) stop ('Wait, there is a negative standard error!')
  if(sum(dim(point)==dim(se))!=2) stop (paste('Oops,',type,"matrix and its SEs are of different size!"))
  
  
  if(type=='lambda'){
    
    what <- 'Factor Loadings'
    
    ll <- point - crit * se
    ul <- point + crit * se
    
  } else if(type=='Phi'){
    
    what <- 'Factor Correlations'
    
    if(any(point - t(point) > .0001)) stop(paste('It may not be a correlation matrix!'))
    if(any(se - t(se) > .0001)) stop(paste('SE for phi is not symmetric...?'))
    
    p <- dim(point)[1]
    
    phi <- point[lower.tri(point)]
    phi.se <- se[lower.tri(se)]
    z.phi <- log((1 + phi)/(1-phi))/2
    z.phise <- phi.se/(1-phi^2)

    z.phi.cl.low <- z.phi - crit * z.phise
    z.phi.cl.high <- z.phi + crit * z.phise
    
    phi.cl.low <- (exp(2*z.phi.cl.low)-1)/(exp(2*z.phi.cl.low)+1)
    phi.cl.high <- (exp(2*z.phi.cl.high)-1)/(exp(2*z.phi.cl.high)+1)
        
    ll <- matrix(0, ncol = p, nrow = p)
    ll[lower.tri(ll)] <- phi.cl.low
    ll <- ll + t(ll)+diag(p)
    ul <- matrix(0, ncol = p, nrow = p)
    ul[lower.tri(ul)] <- phi.cl.high
    ul <- ul + t(ul)+diag(p)
    
  } else if(type=='h'){
    
    what <- 'Communalities'
    
    h <- point
    z.h <- -log(1/h-1)
    z.hse <- se/(h*(1-h))
    
    z.h.cl.low <- z.h - crit * z.hse
    z.h.cl.high <- z.h + crit * z.hse
    
    ll <- 1/(1+exp(-z.h.cl.low))
    ul <- 1/(1+exp(-z.h.cl.high))
    
  } else if(type=='uv'){
    
    what <- 'Unique Variances'
    
    uvar <- point
    z.uvar <- log(uvar)
    z.uvarse <- se / uvar

    z.uvar.cl.low <- z.uvar - crit * z.uvarse
    z.uvar.cl.high <- z.uvar + crit * z.uvarse
    
    ll <- exp(z.uvar.cl.low)
    ul <- exp(z.uvar.cl.high)
    
  } else stop('Hmm, specification is not accepted!')
  
  out <- list(#What = what,
              #Point.Estimate = point,
              #SE.Estimate = se,
              LowerLimit = ll,
              UpperLimit = ul)#, 
              #Alpha = alpha)
  return(out)
}

