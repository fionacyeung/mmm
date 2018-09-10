logitProb_latent_opp_set = function(beta, GammaW, GammaM, Xd, Zd, gw, gm, n, symmetric, asc){
  
  Xdim = dim(Xd)
  Zdim = dim(Zd)
  nXu = Xdim[1]
  nZu = Xdim[2]

  num_var_w = Xdim[3]
  num_var_m = Zdim[3]
  
  # assumes 0 for outside option
  if (asc==1) { # same asc for all types of male (Menzel's)
    if (symmetric) {
      beta_w = beta[1:num_var_w]
      beta_m = beta_w
      asc_w = rep(beta[num_var_w+1], nZu)
      asc_m = rep(beta[num_var_w+1], nXu)
    } else {
      beta_w = beta[1:num_var_w]
      beta_m = beta[(num_var_w+1)+1:num_var_m]
      asc_w = rep(beta[num_var_w+1], nZu)
      asc_m = rep(beta[num_var_w+1+num_var_m+1], nXu)
    }
  } else {  # assume a separate asc for each type of male if asc != 1
    if (symmetric) {
      beta_w = beta[1:num_var_w]
      beta_m = beta_w
      asc_w = beta[num_var_w+(1:nZu)]
      asc_m = asc_w   # assumes nXu==nZu
    } else {
      beta_w = beta[1:num_var_w]
      beta_m = beta[(num_var_w+nZu)+1:num_var_m]
      asc_w = beta[num_var_w+(1:nZu)]
      asc_m = beta[(num_var_w+nZu+num_var_m)+(1:nXu)]
    }
  }
  
  CPW = matrix(0, nrow=nXu, ncol=nZu+1) # women decision makers indexed by row
  CPM = matrix(0, nrow=nZu, ncol=nXu+1) # men decision makers indexed by row
  
  # single option
  for (i in 1:nXu) {
    CPW[i,(1+nZu)] = 1/(1+GammaW[i])
  }
  for (i in 1:nZu) {
    CPM[i,(1+nXu)] = 1/(1+GammaM[i])
  }
  
  # married women
  for (j in 1:nXu) {
    for (k in 1:nZu) {
      Ustar = asc_w[k] # alternative-specific constant
      for (i in 1:num_var_w) {
        Ustar = Ustar + beta_w[i]*Xd[j,k,i]
      }
      CPW[j,k] = 1/sqrt(n) * exp(Ustar) / (1.0 + GammaW[j])
    }
  }
  
  # married men
  for (j in 1:nZu) {
    for (k in 1:nXu) {
      Vstar = asc_m[k]
      for (i in 1:num_var_m) {
        Vstar = Vstar + beta_m[i]*Zd[j,k,i] # alternative-specific constant
      }
      CPM[j,k] = 1/sqrt(n) * exp(Vstar) / (1.0 + GammaM[j])
    }
  }
  
  return(list(CPW = CPW, CPM = CPM))

}


