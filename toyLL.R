#---------------------------------------------------------------------
#  toySWIFT-LL
#  (Version 2.0, October 11, 2022)
#  (c) Ralf Engbert & Maximilian M. Rabe, Universität Potsdam
#---------------------------------------------------------------------#---------------------------------------------------------------------
toyLL <- function(lfreq,fixpos,fixdur,nu,r,mt,iota,beta) {
  
  # fixed parameters
  eta = -2
  kappa = 0
  
  # model variables
  sigma = 1/(1+2*nu+nu^2)   # normalization of processing span
  NW = length(lfreq)     # number of words
  amax = 1 - beta*lfreq  # word frequency dependent maximum activation
  a = rep(0,NW)  # word activation
  s = rep(0,NW)  # word saliency
  p = rep(0,NW)  # selection probability
  shape = 9    # shape parameter of the gamma distribution
  rate = shape/mt   # rate parameter
  LLtime = LLspat = 0   # log-likelihood components
  
  # simulation loop
  Nfix = length(fixdur)
  tfix = fixdur[1]     # fixation duration
  k = fixpos[1]        # fixated word
  for ( j in 1:(Nfix-1) ) {
    
    # 1. Update processing rates
    lambda = rep(0, NW)
    if(k-1 >= 1) lambda[k-1] = nu*sigma
    lambda[k] = sigma
    if(k+1 <= NW) lambda[k+1] = nu*sigma
    if(k+2 <= NW) lambda[k+2] = nu^2*sigma
    
    # 2. Evolve activations
    a = a + r*lambda*tfix/1000
    idx = a>amax
    a[idx] = amax[idx]

    # Compute word saliencies
    s = amax*sin(pi*a/amax) + 10^eta

    # Compute probability for target selection 
    p = s/sum(s)      
    
    # 3. Spatial loglik
    k = fixpos[j+1]   # fixation position observed in the data
    LLspat = LLspat + log(p[k])
    
    # 4. Temporal loglik
    tfix = fixdur[j+1]   # fixation duration observed in the data
    if ( k>1 )  leftact = sum(s[1:(k-1)])
    else leftact = 0
    rate2 = rate*(1+iota*a[k])/(1+kappa*leftact)
    LLtime = LLtime + dgamma(tfix,shape,rate2,log=TRUE)
  }
  loglik = sum(LLtime,LLspat)
  return(loglik)
}