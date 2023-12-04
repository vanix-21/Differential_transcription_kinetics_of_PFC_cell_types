library(ggplot2)
library(BPSC)
library(openssl)

##plot theme
theme_gen<- theme(axis.title.x = element_text(face="bold", size=28),
                  axis.text.x  = element_text(face= "bold", color = "black", size=26),
                  axis.title.y = element_text(face="bold", size=28),
                  axis.text.y  = element_text(face= "bold", color = "black", size=26),
                  axis.line = element_line(linewidth=1.5),
                  axis.ticks = element_line(linewidth = 0.8),
                  strip.text = element_text(face = "bold", size = rel(1.8)),
                  legend.title= element_blank(),
                  legend.text = element_text(face= "bold", color = "black", size=26))


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

betalabelt= c('beta  =')##label for plotting
alphalabelt=c('alpha =')##label for plotting

### SET ####
GOI_data_sim <- read.delim("GOI_data_sim")### alpha and beta parameter data of selected genes
set.seed(228)
betamax=999 ## maximal value of beta in the whole dataset
alphamax=2 ## maximal value of alpha in the whole dataset
####

##simulation and plotting
j=1
for (j in 1:nrow(GOI_data_sim)){
  gene=GOI_data_sim[j, 'Gene.name'] ##gene name
  print(gene)
  
  ###Pyr data
  alpha=GOI_data_sim[j, "alpha_pyr"] 
  beta=GOI_data_sim[j, "beta_pyr"]
  dfpyr<- BurstSim(alpha, beta, alphamax , betamax)  ##simulation
  
  alphalabeln=c(round(alpha, digits = 4)) ##label for plot
  betalabeln=c(round(beta, digits = 4)) ##label for plot
  
  ggplot(dfpyr, aes(x=tvals, y=rnavals))+
    geom_line(linewidth=1)+
    theme_classic() +
    xlab(bquote("Time"))+
    ylab(bquote("Transcription\nrate"))+
    scale_y_continuous(limits=c(0,10))+
    annotate("text", x = -Inf, y = Inf,  label = alphalabelt, hjust = -.2, vjust =2.5 , size= 8, fontface = "bold")+
    annotate("text", x = -Inf, y = Inf,  label = alphalabeln, hjust = -2.5, vjust =2.5 , size= 8, fontface = "bold")+
    annotate("text", x = -Inf, y = Inf,  label = betalabelt, hjust = -.2, vjust =4 , size= 8, fontface = "bold")+
    annotate("text", x = -Inf, y = Inf,  label = betalabeln, hjust = -2, vjust =4, size= 8, fontface = "bold")+
    annotate("text", x = -Inf, y = Inf,  label = gene, hjust = -.2, vjust =1 , size= 8, fontface = "bold")+
    theme_gen
  
  name=paste('pyr', gene, '.tiff', sep='_') #file name
  ggsave(name, units="mm", width=300, height=100, dpi=600, compression = 'lzw')  ##save to wd
  
  ##FS data
  alpha=GOI_data_sim[j, "alpha_fs"]
  beta=GOI_data_sim[j, "beta_fs"]
  dffs<- BurstSim(alpha, beta, alphamax , betamax) ##simulation

  alphalabeln=c(round(alpha, digits = 4))##label for plot
  betalabeln=c(round(beta, digits = 4))##label for plot
  
  ggplot(dffs, aes(x=tvals, y=rnavals))+
    geom_line(linewidth=1)+
    theme_classic() +
    xlab(bquote("Time"))+
    ylab(bquote("Transcription\nrate"))+
    scale_y_continuous(limits=c(0,10))+
    annotate("text", x = -Inf, y = Inf,  label = alphalabelt, hjust = -.2, vjust =2.5 , size= 8, fontface = "bold")+
    annotate("text", x = -Inf, y = Inf,  label = alphalabeln, hjust = -2.5, vjust =2.5 , size= 8, fontface = "bold")+
    annotate("text", x = -Inf, y = Inf,  label = betalabelt, hjust = -.2, vjust =4 , size= 8, fontface = "bold")+
    annotate("text", x = -Inf, y = Inf,  label = betalabeln, hjust = -2, vjust =4 , size= 8, fontface = "bold")+
    annotate("text", x = -Inf, y = Inf,  label = gene, hjust = -.2, vjust =1 , size= 8, fontface = "bold")+
    theme_gen
  
  name=paste('fs', gene, '.tiff', sep='_') #file name
  ggsave(name, units="mm", width=300, height=100, dpi=600, compression = 'lzw') #save to wd
}




