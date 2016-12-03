pi_to_p <-
function(s, t, rho, tau)
{
  pi_st_to_p_jk = phi_pdf(s, t, rho, tau) - phi_pdf(s, t-1, rho, tau) - phi_pdf(s-1, t, rho, tau) + phi_pdf(s-1, t-1, rho, tau)
  return(pi_st_to_p_jk)
}

