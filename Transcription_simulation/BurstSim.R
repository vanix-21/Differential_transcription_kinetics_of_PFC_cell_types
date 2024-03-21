
library(stats)
##function for simulation
#alpha: parameter representing burst frequency
#beta: parameter representing burst size
#alphamax: maximal value of alpha in the whole dataset
#betamax: maximal value of beta in the whole dataset
BurstSim<- function(alpha, beta, alphamax, betamax){
  tmax=1 #Maxtime 
  kon=alpha/alphamax #rate of switching on
  onprob=1-exp(-kon*tmax) #on state probability
  offprob=1-onprob #off state probability
  koff=(-log(1-offprob, base = exp(1)))/tmax  #rate of switching off
  ks=(beta/betamax)*koff #RNA synthesis rate
  
  rna0=0 #initial rna
  curr=rna0 #Current rna
  currt=0 #Current time 
  df<- data.frame(tvals=currt, rnavals=curr)#Recording of time and rna
  state=0 #initial state (0:off; 1:on)
  
  while(currt<tmax){
    onrate=kon
    offrate=koff
    total_rate=onrate+offrate
    
    deltat=rexp(1, rate=1/total_rate) #Find next event 
    currt=currt+(deltat/1000) #Move forward in time 
    
    #switch on
    randon = runif(1) #random number from unifrom distribution
    if(state==0 & randon<onrate) {
      state=1
      if (beta>10) curr=5 else curr=3
    }
    
    #random peaks at on state
    rand=runif(1)
    if (state == 1 & rand<0.5 & curr<=7 & curr>5) curr=curr-1 
    if (state == 1 & rand>=0.5 & curr<7 & curr>=3) curr=curr+1 
    
    #switch off
    rands= runif(1, max = ks) # rands depends on RNA synthesis rate and beta
    if (beta>10) rands= runif(1, max = ks*10) 
    if (beta>100) rands= runif(1, max = ks*100)
    
    if (rands<ks) randoff = runif(1, max = offrate*10) else randoff = runif(1, max = offrate*100)
    if(state==1 & randoff<offrate ){
      curr = 0
      state=0
    } #the larger the beta there is less probability to switch off
    
    #Recordthestate   
    df<- rbind(df, c(currt, curr))
  }
  return(df)
}




