
library(ggplot2)
library(BPSC)
library(openssl)

source("BurstSim.R")
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



### SET ####

GOI_data_sim <- read.delim("GOI_data_sim")### alpha and beta parameter data of selected genes
set.seed(228)
betamax=999 ## maximal value of beta in the whole dataset
alphamax=2 ## maximal value of alpha in the whole dataset
####
betalabelt= c('beta  =')##label for plotting
alphalabelt=c('alpha =')##label for plotting

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




