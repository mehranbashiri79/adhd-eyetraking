#---------------------------------------------------------------------
#  toygen.R
#  Function for simulation of toySWIFT (generative model)
#  (Version 3.0, January 3, 2023)
#  (c) Ralf Engbert & Maximilian M. Rabe, Universität Potsdam
#---------------------------------------------------------------------
toygen <- function(nu=0.3,r=10,mt=200,iota=0.5,eta=-3,beta=0.6,kappa=0,gamma=1,lfreq,MODE=2) {
  # activation, saliency, probability
  NW = length(lfreq)     # number of words
  amax = 1 - beta*lfreq  # word frequency dependent maximum activation
  a = rep(0,NW)  # word activation
  s = rep(0,NW)  # word saliency
  p = rep(0,NW)  # selection probability
  # fixation duration 
  shape = 9                  # gamma density shape
  rate = shape/mt            # gamma density rate
  sigma = 1/(1+2*nu+nu^2)    # processing span normalization constant
  # other variables
  time = 0    # time
  k = 1       # fixated word
  traj = c()  # store trajectory
  
  # simulation loop
  while ( length(which(a<1))>0 & k<NW ) {
    
    # 1. Generate fixation duration
    if ( k>1 )  leftact = prod(1+kappa*(s[1:(k-1)]-10^eta))
    else leftact = 1
    tfix = rgamma(1,shape,(rate*(1+iota*a[k])/leftact))
    
    # 2. Update processing rates
    lambda = rep(0, NW)
    if(k-1 >= 1) lambda[k-1] = nu*sigma
    lambda[k] = sigma
    if(k+1 <= NW) lambda[k+1] = nu*sigma
    if(k+2 <= NW) lambda[k+2] = nu^2*sigma
    
    # 3. Evolve activations
    # MODE=1: millisecond computation
    if ( MODE==1 ) {
      for ( t in 1:tfix ) {
        time = time + 1
        a = a + r*lambda/1000
        idx = a>amax
        a[idx] = amax[idx]
        # Compute word saliencies
        s = amax*sin(pi*a/amax) + 10^eta
        traj = rbind(traj,c(time,k,tfix,s))
      }
    }  else  {
      # MODE=0/2: computation per fixation
      
      # 2. Update processing rates
      lambda = rep(0, NW)
      if (k-1 >= 1) lambda[k-1] = nu*sigma
      lambda[k] = sigma
      if (k+1 <= NW) lambda[k+1] = nu*sigma
      if (k+2 <= NW) lambda[k+2] = nu^2*sigma
      
      # 3. Evolve activations
      time = time + tfix
      a = a + r*lambda*tfix/1000
      idx = a>amax
      a[idx] = amax[idx]
      
      # Compute word saliencies 
      s = amax*sin(pi*a/amax) + 10^eta
      
      # store trajectory
      traj = rbind(traj,c(time,k,tfix,s))
    }
    
    # Compute probability for target selection
    p = s^gamma/sum(s^gamma)      
    
    # 6. Select saccade target
    k = sample.int(n=NW, size=1, prob=p) 
    
    # 4. Update fixation position
    if ( MODE==0 )  traj = rbind(traj,c(time,k,tfix,s))
    
  }
  if ( k>1 )  leftact = prod(1+kappa*(s[1:(k-1)]-10^eta))
  else leftact = 1
  tfix = rgamma(1,shape,(rate*(1+iota*a[k])/leftact))
  if ( MODE==1 ) traj = rbind(traj,c(time,k,tfix,rep(0,NW)))
  else  traj = rbind(traj,c(time,k,tfix,s))
  
  return(traj)
}