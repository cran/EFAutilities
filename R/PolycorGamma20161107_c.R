## ==============================================================================
## Estimating polychoric correlations and its asymptotic Covariance matrix
## ==============================================================================
## Nov. 3rd: Speed up the program by rewriting for loops with matrix operation
## Oct. 28th: Muted three lines of testing codes 
## Aug. 26th: Compile all functions into 'get.RGamma'
## Aug. 23th: Fisher scoring method is used to compute polychoric correation (faster)
## Aug. 18th: Update the functions check.data because of the naming error (X, Y, etc.)
## ==============================================================================
library("mvtnorm")
get.RGamma = function(dat, gamma=FALSE)
{
    #library("mvtnorm")
  # test if data set contains continuous/constant variables
  cont = apply(dat, 2, function(x18) sum(abs(x18-round(x18))))
  if (any(cont!=0)){
    stop("data set contains continuous variables: ", paste0("V",which(cont!=0),collapse=','))
  }
  Vari = apply(dat, 2, var)
  if (any(Vari==0)){
    stop("data set contains variables with zero variance: ", paste0("V",which(Vari==0),collapse=','))
  }
  #cat("All clear!",'\n')
  
  
  ## =============================================
  ## A list of relevant functions
  ## =============================================
  ## bivariate normal cdf over certain intervals / pi_st
  pi_st = function(s, t, rho, tau)
  {	
    tau[tau<=-1000] = -Inf;  tau[tau>=1000] = Inf; 
    p = pmvnorm(lower = c(tau[1,s],tau[2,t]), upper = c(tau[1,s+1],tau[2,t+1]), mean = c(0,0), sigma = matrix(c(1, rho, rho, 1), 2))    
    return(p)		
  }
  pi_st2 = function(pars)
  {	
    pars[pars<=-1000] = -Inf;  pars[pars>=1000] = Inf; 
    p = pmvnorm(lower = c(pars[1],pars[2]), upper = c(pars[3],pars[4]), mean = c(0,0), sigma = matrix(c(1, pars[5], pars[5], 1), 2))    
    return(p)		
  }
  ## bivariate normal pdf 
  phi_pdf = function(x1, x2, rho)
  {
    z1 = x1^2 + x2^2 - 2*rho*x1*x2
    z2 = exp(-z1/2/(1-rho*rho))
    p = z2/2/pi/sqrt(1-rho*rho)
    return(p)		
  }
  ## 1st order derivative  ## d(pi_st)/d(p_jk)
  pi_to_p = function(s, t, rho, tau)
  {
    x1 = tau[1,s+1]; x2 = tau[2,t+1]; x3 = tau[1,s]; x4 = tau[2,t]; 
    pi_st_to_p_jk = phi_pdf(x1, x2, rho) - phi_pdf(x1, x4, rho) - phi_pdf(x3, x2, rho) + phi_pdf(x3, x4, rho)
    return(pi_st_to_p_jk)
  }
  
  ## transfer GammaRho to GammaR
  get.covR = function(GammaRho, p)
  {
    cloc = diag(0, p, p); cloc[lower.tri(cloc)] = 1:ncol(GammaRho); cloc = cloc + t(cloc)
    ccloc = c(cloc)    #correlation location
    
    all.ind = as.matrix(expand.grid(ccloc, ccloc))
    zero.ind = which(rowSums(all.ind==0)==0)
    value = rep(0, nrow(all.ind))
    value[zero.ind] = GammaRho[all.ind]
    
    GammaR = matrix(value, p^2, p^2)
    
    return(GammaR)
  }
  ## merge the estimated thresholds allowing for unequal length
  mg_filled <- function(ls){
    max_length <- max(unlist(lapply(ls,length)))
    nm_filled <- sapply(ls, function(x) c(x, rep(NA, max_length - length(x))))
    return(t(nm_filled))
  }
  ## compute g6 to be used in estimation equations
  get.g6 = function(dat, R, n, p, catt, hi){
    dm = apply(rbind(1:p, dat), 2, function(x3) match(x3, catt[[x3[1]]]))[-1,]+1
    
    rho.ls = R[lower.tri(R)]
    bran = combn(1:p,2)
    ele = Reduce(rbind, lapply(1:n, function(i) t(combn(dm[i,],2))))
    
    ab.1 = cbind(hi[cbind(bran[1,],ele[,1])], hi[cbind(bran[2,],ele[,2])])
    ab.2 = cbind(hi[cbind(bran[1,],ele[,1]-1)], hi[cbind(bran[2,],ele[,2])])
    ab.3 = cbind(hi[cbind(bran[1,],ele[,1])], hi[cbind(bran[2,],ele[,2]-1)])
    ab.4 = cbind(hi[cbind(bran[1,],ele[,1]-1)], hi[cbind(bran[2,],ele[,2]-1)])
    
    num = phi_pdf(ab.1[,1], ab.1[,2], rho.ls) - phi_pdf(ab.2[,1], ab.2[,2], rho.ls) - phi_pdf(ab.3[,1], ab.3[,2], rho.ls) + phi_pdf(ab.4[,1], ab.4[,2], rho.ls)
    
    den.inp = cbind(ab.4, ab.1, rho.ls)
    den = apply(den.inp, 1, pi_st2)
    g6 = matrix(num/den, nrow = n, byrow = T)
    return(g6=g6)
  }

  ## Compute polychoric correlation matrix
  fisher.polychor = function(dat, hat.tau){
    #cat('Compute Polychoric Correlation...', '\n')
    n = nrow(dat); p = ncol(dat); 
    ncat = unname(apply(dat, 2, function(x18) length(unique(x18)))) 
    
    ep = 1E-5 # Fisher-scoring criterion
    itr = 50 ##Max number of iteration
    idiv = 0 ##idiv is the nonconvergence due to bad fit
    
    R = diag(1, p)
    
    scorr = cor(dat) #use pearson's correlations as initial values
    
    for (j in 1 : (p-1)){
      for (k in (j+1) : p){
        
        p_jk = scorr[j,k]
        tau = hat.tau[c(j,k),]
        ncount = table(dat[,j], dat[,k])
        if (min(ncount)==0){ncount = ncount + 1/ncat[j]/ncat[k]} #empty cells in the contingency table
        E_l_dotdot_p_jk = l_dot_p_jk = pii_st = matrix(NA, ncat[j], ncat[k])   
        
        for (i in 1 : itr){
          st = expand.grid(1:ncat[j], 1:ncat[k])
          pi_num = matrix(pi_to_p(st[,1], st[,2], p_jk, tau), ncat[j], ncat[k])
          
          for (s in 1 : ncat[j]){
            for (t in 1 : ncat[k]){
              
              pii_st[s,t] = pi_st(s, t, p_jk, tau)
              l_dot_p_jk[s,t] = ncount[s,t] * pi_num[s,t] / pii_st[s,t]
              E_l_dotdot_p_jk[s,t] = (pi_num[s,t]^2) / pii_st[s,t]		
            }
          }   
          
          l_dot = sum(l_dot_p_jk)
          E_l_dotdot = -n * sum(E_l_dotdot_p_jk)
          
          dt = l_dot / E_l_dotdot
          p_jk = p_jk - dt   
          
          if (abs(dt)<ep) break
          if (abs(p_jk)>1)  break  
        }
        if (i==itr|abs(p_jk)>1) cat('Cannot estimate polychoric correlation for V',j,' and V',k,' using Fisher Scoring', '\n',sep = '')
        if (i==itr|abs(p_jk)>1) break
        R[j,k] = R[k,j] = p_jk
      }
    }
    return(R)
  }  
  
  
  n = nrow(dat); p = ncol(dat); 
  dat = as.matrix(dat, n, p)
  ncat = unname(apply(dat, 2, function(x18) length(unique(x18))))
  catt = lapply(1:p, function(x18) sort(unique(dat[,x18])))
  # estimated thresholds
  htau.null = lapply(1:p, function(x18) head(qnorm( cumsum(table(dat[,x18])/n) ), -1))
  hat.tau = mg_filled(lapply(1:p, function(x18) c(-1E+30, head(qnorm( cumsum(table(dat[,x18])/n) ), -1), 1E+30)))
  
  ## compute polychoric correlations
  R = fisher.polychor(dat, hat.tau)
  
  ## No output of covariance matrix
  if (!gamma){return(list(cutoffs=htau.null, R=R))}
  
  ## ================================
  ###  Covariance matrix of R  ###
  ## ================================
  if (gamma){
    #cat('Gamma Matrix Estimation...', '\n')
    
    g3 = NULL
    for(t in 1 : p) { 
      g3 = cbind(g3, sapply(1:(length(catt[[t]])-1), function(x18) ifelse(dat[,t]<=catt[[t]][x18], 1, 0)) - matrix(rep(pnorm(htau.null[[t]]), n), n, byrow=T))
    }
    
    g6 = get.g6(dat, R, n, p, catt, hat.tau)
    G = cbind(g3, g6)
    
    
    B = matrix( rowMeans(apply(G, 1, function(x18) as.matrix(x18)%*%t(x18))), ncol(G), byrow=T)
    
    A33 = -diag(dnorm(unlist(htau.null)), ncol(g3))
    
    A63 = matrix(0, ncol(g6), ncol(g3))
    A66 = matrix(0, ncol(g6), ncol(g6))
    loc = list(); for (t in 1:p){loc[[t]]=seq((c(0,cumsum(ncat-1))[t]+1),cumsum(ncat-1)[t])}
    
    ini = 0
    for (j in 1 : (p-1)){
      for (k in (j+1) : p){
        ini = ini + 1
        tau = hat.tau[c(j,k),]
        rho = R[j,k]
        store63_js_k = matrix(NA, ncat[j]-1, ncat[k]); store63_j_kt = matrix(NA, ncat[k]-1, ncat[j]) 
        store66 = matrix(NA, ncat[j], ncat[k]) 
        
        st = expand.grid(1:ncat[j], 1:ncat[k])
        pi_num = matrix(pi_to_p(st[,1], st[,2], rho, tau), ncat[j], ncat[k])
        pi_den = matrix(NA, ncat[j], ncat[k])
          
        for (s in 1 : ncat[j]){
          for (t in 1 : ncat[k]){
            num5 = pi_num[s,t]
            pi_den[s,t] = den5 = pi_st(s, t, rho, tau)
            store66[s,t] = num5^2/den5 
          }
        }
        
        for (s in 1 : (ncat[j]-1)){
          for (t in 1 : ncat[k]){
            num1 = pi_num[s,t]
            den1 = pi_den[s,t]
            num2 = pi_num[s+1,t]
            den2 = pi_den[s+1,t]
            z1 = (tau[2,t+1]-rho*tau[1,s+1])/sqrt(1-rho^2)
            z2 = (tau[2,t]-rho*tau[1,s+1])/sqrt(1-rho^2)
            mul = (pnorm(z1)-pnorm(z2))*dnorm(tau[1,s+1])
            store63_js_k[s,t] = (num1/den1 - num2/den2)*mul
          }
        }
        
        for (t in 1 : (ncat[k]-1)){
          for (s in 1 : ncat[j]){
            num3 = pi_num[s,t]
            den3 = pi_den[s,t]
            num4 = pi_num[s,t+1]
            den4 = pi_den[s,t+1]
            z3 = (tau[1,s+1]-rho*tau[2,t+1])/sqrt(1-rho^2)
            z4 = (tau[1,s]-rho*tau[2,t+1])/sqrt(1-rho^2)
            mul = (pnorm(z3)-pnorm(z4))*dnorm(tau[2,t+1])
            store63_j_kt[t,s] = (num3/den3 - num4/den4)*mul
          }
        }    

        A63[ini,loc[[j]]] = -rowSums(store63_js_k)
        A63[ini,loc[[k]]] = -rowSums(store63_j_kt)
        A66[ini,ini] = -sum(store66)
      }
    }
    
    ## ===========================================
    ###   Sandwich Covariance A^(-1)BA'^(-1)   ###
    ## ===========================================
    ind = cumsum(c(ncol(g3), ncol(g6)))
    A = matrix(0, ncol(G), ncol(G))
    A[1:ind[1], 1:ind[1]] = A33
    A[(ind[1]+1):ind[2], 1:ind[1]] = A63
    A[(ind[1]+1):ind[2], (ind[1]+1):ind[2]] = A66
    A.inv = solve(A)
    
    Gamma = A.inv%*%B%*%t(A.inv)
    GammaRho = Gamma[(ind[1]+1):ind[2], (ind[1]+1):ind[2]]
    GammaR = get.covR(GammaRho, p)
    return(list(cutoffs=htau.null, R=R, GammaR=GammaR))
  } 
}
