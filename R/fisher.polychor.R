fisher.polychor <-
function(dat, hat.tau){
  cat('Compute Polychoric Correlation...', '\n')
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
      tau = hat.tau[c(j,k)]
      ncount = table(dat[,j], dat[,k])
      if (min(ncount)==0){ncount = ncount + 1/ncat[j]/ncat[k]} #empty cells in the contingency table
      E_l_dotdot_p_jk = l_dot_p_jk = pii_st = matrix(NA, ncat[j], ncat[k])   
      
      for (i in 1 : itr){
        for (s in 1 : ncat[j]){
          for (t in 1 : ncat[k]){
            
            pii_st[s,t] = pi_st(s, t, p_jk, tau)
            l_dot_p_jk[s,t] = ncount[s,t] * pi_to_p(s, t, p_jk, tau) / pii_st[s,t]
            E_l_dotdot_p_jk[s,t] = (pi_to_p(s, t, p_jk, tau)^2) / pii_st[s,t]		
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

