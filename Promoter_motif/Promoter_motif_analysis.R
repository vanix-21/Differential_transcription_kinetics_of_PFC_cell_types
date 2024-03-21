library(ggplot2)
library(forcats)
library(dplyr)
library(tidyr)
library(multcomp)

theme_gen<- theme(axis.title.x = element_text(face="bold", size=26),
                  axis.text.x  = element_text(face= "bold", color = "black", size=24),
                  axis.title.y = element_text(face="bold", size=26),
                  axis.text.y  = element_text(face= "bold", color = "black", size=24),
                  axis.line = element_line(linewidth=1.5),
                  axis.ticks = element_line(linewidth = 0.8),
                  strip.text = element_text(face = "bold", size = rel(1.8)),
                  legend.title= element_blank(),
                  legend.text = element_text(face= "bold", color = "black", size=12))


df_20percent_both<- read.csv("burst_kinetics_copynumber.csv")


ccaat<- read.csv("CCAATbox.csv", sep=";")
tata<- read.csv("TATA.csv", sep=";")
gc<- read.csv("GCbox.csv", sep=";")

#adding promoter information to genes
merged_promoters_bk<-data.frame()
merged_promoters_bk<- merge(df_20percent_both, 
                            tata[tata$alternative_promoters==1,8:9], by= 'Gene.name', all.x = TRUE)
merged_promoters_bk<- merge(merged_promoters_bk, 
                            ccaat[ccaat$alternative_promoters==1,8:9], by= 'Gene.name', all.x = TRUE)
merged_promoters_bk<- merge(merged_promoters_bk, 
                            gc[gc$alternative_promoters==1,8:9], by= 'Gene.name', all.x = TRUE)


##Groups based on a combination of promoters
merged_promoters_bk$Group<- paste(merged_promoters_bk$CCAAT.box , merged_promoters_bk$Gcbox,
                                  merged_promoters_bk $TATA)

merged_promoters_bk[merged_promoters_bk$Group=='NA NA NA', 'Group']='No motif'
merged_promoters_bk[merged_promoters_bk$Group=='NA Gcbox TATA', 'Group']='GC-box and TATA'
merged_promoters_bk[merged_promoters_bk$Group=='CCAAT-box Gcbox NA', 'Group']='GC-box and CCAAT-box'
merged_promoters_bk[merged_promoters_bk$Group=='NA Gcbox NA', 'Group']='Only GC-box'
merged_promoters_bk[merged_promoters_bk$Group=='CCAAT-box NA NA', 'Group']='Only CCAAT-box'
merged_promoters_bk[merged_promoters_bk$Group=='NA NA TATA', 'Group']='Only TATA-box'
merged_promoters_bk[merged_promoters_bk$Group=='CCAAT-box Gcbox TATA', 'Group']='GC-box and CCAAT-box and TATA'
merged_promoters_bk[merged_promoters_bk$Group=='CCAAT-box NA TATA', 'Group']='CCAAT-box and TATA'

library(RColorBrewer) 

#
brewer.pal(n=12,"Set2")#"#FFFF99"

promoter_col<- c("No motif"="white"       , 
                  "GC-box and CCAAT-box" = "#A6CEE3" ,     
                  "Only GC-box"=   "#1F78B4" ,  
                  "Only TATA-box" = "#E31A1C", 
                  "GC-box and CCAAT-box and TATA" = "#FB9A99",
                  "GC-box and TATA" = "#CAB2D6",
                  "Only CCAAT-box"= "#33A02C"  ,
                  "CCAAT-box and TATA"=  "#FFFF99"   )

promoter_col2<- c("No motif"="white"       , 
                  "Gcbox"=   "#A6CEE3"  ,  
                  "TATA" = "#CAB2D6" , 
                  "CCAAT-box"=  "#FFFF99" )


promoters_pivot_longer<- pivot_longer(merged_promoters_bk, cols = c('TATA', 'Gcbox', 'CCAAT.box') , values_to = 'promoter')
promoters_pivot_longer[promoters_pivot_longer$Group=='No motif' & is.na(promoters_pivot_longer$Group)==FALSE , 'promoter']<- 'No motif'
promoters_pivot_longer$promoter<- factor(promoters_pivot_longer$promoter, levels = c('No motif',
                                                                                     'Gcbox', 
                                                                                     'CCAAT-box',
                                                                                     'TATA'))

merged_promoters_bk<- merge(merged_promoters_bk, df_20percent_both[,c('Gene.name', 'VAR')], by='Gene.name', all.x=TRUE)
merged_promoters_bk<- unique(merged_promoters_bk)
promoters_pivot_longer<-promoters_pivot_longer[is.na(promoters_pivot_longer$promoter)==FALSE, ] 

##beta 
ggplot(promoters_pivot_longer, aes(y=log2(beta), x=promoter, fill=promoter))+
  geom_boxplot( linewidth=1, width=0.7)+
  theme_minimal()+
  theme_gen+
  scale_fill_manual(values=promoter_col2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+    
  ylab(bquote('Burst size (beta) (log2)'))+
  xlab(bquote(''))

anova(lm(promoters_pivot_longer$beta ~promoters_pivot_longer$promoter))
TukeyHSD(aov(beta ~ promoter, data=promoters_pivot_longer))


##alpha
ggplot(promoters_pivot_longer, aes(y=log2(alpha), x=promoter, fill=promoter))+
  geom_boxplot( linewidth=1, width=0.7)+
  theme_minimal()+
  theme_gen+
  scale_fill_manual(values=promoter_col2)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+    
  ylab(bquote('Burst frequency (alpha) (log2)'))+
  xlab(bquote(''))


anova(lm(promoters_pivot_longer$alpha ~promoters_pivot_longer$promoter))
TukeyHSD(aov(alpha ~ promoter, data=promoters_pivot_longer))



merged_promoters_bk[is.na(merged_promoters_bk$Group )== FALSE &
                      merged_promoters_bk$Group!='CCAAT-box and TATA' & merged_promoters_bk$Group!='GC-box and CCAAT-box and TATA',] %>%
  mutate(Group = fct_reorder(Group, alpha, .na_rm = TRUE)) %>%
  ggplot(aes(x=Group, y=log2(alpha), fill=Group))+
  geom_boxplot( linewidth=1, width=0.7)+
  theme_minimal()+
  theme_gen+
  scale_fill_manual(values=promoter_col)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(face = "bold", size = rel(1.5)))+    
  ylab(bquote('Burst frequency (alpha) (log2)'))+
  xlab(bquote(''))


merged_promoters_bk[is.na(merged_promoters_bk$Group )== FALSE &
                      merged_promoters_bk$Group!='CCAAT-box and TATA' & merged_promoters_bk$Group!='GC-box and CCAAT-box and TATA',] %>%
  mutate(Group = fct_reorder(Group, beta, .na_rm = TRUE)) %>%
  ggplot(aes(x=Group, y=log2(beta), fill=Group))+
  geom_boxplot( linewidth=1, width=0.7)+
  theme_minimal()+
  theme_gen+
  scale_fill_manual(values=promoter_col)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(face = "bold", size = rel(1.5)))+    
  ylab(bquote('Burst size (beta) (log2)'))+
  xlab(bquote(''))


anova(lm(merged_promoters_bk[is.na(merged_promoters_bk$Group )== FALSE &
                               merged_promoters_bk$Group!='CCAAT-box and TATA' & 
                               merged_promoters_bk!='GC-box and CCAAT-box and TATA','alpha']~
           merged_promoters_bk[is.na(merged_promoters_bk$Group) == FALSE &
                                 merged_promoters_bk$Group!='CCAAT-box and TATA' & 
                                 merged_promoters_bk!='GC-box and CCAAT-box and TATA','Group'])) ##p= 0.002511


anova(lm(merged_promoters_bk[is.na(merged_promoters_bk$Group )== FALSE &
                               merged_promoters_bk$Group!='CCAAT-box and TATA' & 
                               merged_promoters_bk!='GC-box and CCAAT-box and TATA','beta']~
           merged_promoters_bk[is.na(merged_promoters_bk$Group) == FALSE &
                                 merged_promoters_bk$Group!='CCAAT-box and TATA' & 
                                 merged_promoters_bk!='GC-box and CCAAT-box and TATA','Group'])) ##0.3905





promoter_oprof <- read.csv("promoter_oprof.csv", sep=";")
subset<- c('AlphalargerinFS', 'BetalargerinPyr')
motif<- c('TATA', 'CCAAT', 'GC')
range<- list(c(-79.5:20.5), c(-199.5:20.5), c(-139.5:50.5))

j=1
i=1
for (i in 1:2){
  for (j in 1:3){
    
    spline_int <- as.data.frame(spline(promoter_oprof[promoter_oprof$genes==subset[i] & promoter_oprof$Motif==motif[j], 'Position'],
                                       promoter_oprof[promoter_oprof$genes==subset[i] & promoter_oprof$Motif==motif[j], 'frequency.percent.']))
    
    spline_int <- rbind(spline_int[,1:2], as.data.frame(spline(promoter_oprof[promoter_oprof$genes=='All' & promoter_oprof$Motif==motif[j], 'Position'],
                                                               promoter_oprof[promoter_oprof$genes=='All' & promoter_oprof$Motif==motif[j], 'frequency.percent.'])))
    
    spline_int[1:273, 'group']<- subset[i]
    spline_int[274:546, 'group']<- c('All')
    spline_int$group<- as.factor(spline_int$group)
    
    if (i==1 ) col= c("AlphalargerinFS"="#66c2a5")
    if (i==2) col= c("BetalargerinPyr"= "#fc8d62")
    
    if (j==1 ) col =c(col, "All" = "#E78AC3")
    if (j==2 ) col =c(col, "All" =   "#FFD92F")
    if (j==3 ) col =c(col, "All" = "#8DA0CB")
    
    
    ggplot(promoter_oprof[promoter_oprof$genes==subset[i] & promoter_oprof$Motif==motif[j],])+
      geom_point( aes(x=Position, y=frequency.percent.), alpha=0)+
      geom_line(data=spline_int, aes(x = x, y = y, colour=group), linewidth=2)+
      scale_color_manual(values = col )+
      scale_y_continuous(limits = c(0,65))+
      scale_x_continuous(limits = c(-200,100), breaks = c(-200, -100, 0, 100))+
      theme_minimal()+
      xlab(bquote('Position relative to TSS (bp)'))+
      ylab(bquote('Promoter motif frequency (% of genes)'))+
      #facet_grid(~genes)+
      theme_gen

    name=paste(motif[j], subset[i], '.tiff', sep='_')
    ggsave(name, units="mm", width=200, height=150, dpi=600, compression = 'lzw')
  }
}




for (j in 1:3){
  
  spline_int <- as.data.frame(spline(promoter_oprof[ promoter_oprof$Motif==motif[j], 'Position'],
                                     promoter_oprof[ promoter_oprof$Motif==motif[j], 'frequency.percent.']))
  
  ggplot(promoter_oprof[promoter_oprof$Motif==motif[j],])+
    geom_point( aes(x=Position, y=frequency.percent.), alpha=0)+
    geom_line(data=spline_int, aes(x = x, y = y), colour='grey20', linewidth=1)+
    scale_y_continuous(limits = c(0,65))+
    scale_x_continuous(limits = c(-200,50), breaks = c(-200, -150, -100, -50 ,0, 50))+
    theme_minimal()+
    theme_gen+
    theme(axis.title.x = element_text(face="bold", size=10),
          axis.text.x  = element_text(face= "bold", color = "black", size=8),
          axis.title.y = element_text(face="bold", size=10),
          axis.text.y  = element_text(face= "bold", color = "black", size=8),
          legend.title= element_blank(),
          legend.text = element_blank())+
    xlab(bquote('Position relative to TSS (bp)'))+
    ylab(bquote('Promoter motif frequency (% of genes)'))
  
  name=paste(motif[j], 'All', '.tiff', sep='_')
  ggsave(name, units="mm", width=100, height=100, dpi=600, compression = 'lzw')
  
}

par<- c('alpha', 'beta')
for (i in 1:2){
  for (j in 1:3){
    prom_ks <-ks.test(promoter_oprof[promoter_oprof$genes==subset[i] & promoter_oprof$Motif==motif[j] & promoter_oprof$Position %in% range[[j]],'frequency.percent.'],
                      promoter_oprof[promoter_oprof$genes=='All' & promoter_oprof$Motif==motif[j]& promoter_oprof$Position %in% range[[j]] ,'frequency.percent.'])
    print(c(subset[i], motif[j], prom_ks$p.value ))
  }
}



promoter_col3<- c("GC"=   "#A6CEE3"  ,  
                  "TATA" = "#CAB2D6" , 
                  "CCAAT"=  "#FFFF99" )

promoter_col4<- c("GC"=    "#8DA0CB" ,  
                  "TATA" = "#E78AC3" , 
                  "CCAAT"=   "#FFD92F" )

promoter_integrate <- read.csv("promoter_integrate.csv", sep=";")
promoter_integrate$Genes<-factor(promoter_integrate$Genes, levels = c('All', 'Alpha larger in FS', 'Beta larger in Pyr'))
promoter_integrate$Promoter<-factor(promoter_integrate$Promoter, levels = c('GC', 'CCAAT', 'TATA'))


ggplot(promoter_integrate,aes(x=Genes, y=Area, fill=Promoter, col=Promoter) )+
  geom_col(width=0.6, linewidth=2, position = position_dodge())+
  theme_minimal()+
  scale_fill_manual(values=promoter_col3)+
  scale_colour_manual(values=promoter_col4)+
  scale_y_continuous(breaks = c(1000,2000,3000,4000,5000))+
  theme_gen+
  theme(axis.text.x=element_text(size=22, angle=45, hjust = 1),
        axis.title.x = element_blank())+
  ylab(bquote('AUC of promoter motif frequency'))+
  facet_grid(~Promoter)


