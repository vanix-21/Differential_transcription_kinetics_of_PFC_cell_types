
library(dplyr)
library(matrixStats)
library(readxl)
library(ggplot2)
library(GGally)


df_20percent_both<- read.csv("burst_kinetics_copynumber.csv")
df_20percent_both<- df_20percent_both[,-1]



celltypecol2<- c("Fast-spiking"="black", "Pyramidal" = "grey40")
celltypecol3<- c("Fast-spiking"="white", "Pyramidal" = "grey40")

theme_gen<- theme(axis.title.x = element_text(face="bold", size=28),
                  axis.text.x  = element_text(face= "bold", color = "black", size=26),
                  axis.title.y = element_text(face="bold", size=28),
                  axis.text.y  = element_text(face= "bold", color = "black", size=26),
                  axis.line = element_line(linewidth=1.5),
                  axis.ticks = element_line(linewidth = 0.8),
                  strip.text = element_text(face = "bold", size = rel(1.8)),
                  legend.title= element_blank(),
                  legend.text = element_text(face= "bold", color = "black", size=26))


##correlation of copy number distribition measures and burst kinetics
#burst frequncy and size log2 vs mean expression log2 
ggplot(df_20percent_both, 
       aes(x=log2(MEAN), y=log2(beta), shape=CellType, color = CellType))+
  geom_point(size=2, stroke=2, alpha=0.8)+
  geom_smooth(method=lm, colour= "grey60")+
  theme_minimal()+
  scale_color_manual(values= celltypecol2)+
  scale_shape_manual(values= c(1,2))+
  xlab(bquote('Average of copy number (log2)'))+
  ylab(bquote('Burst size (beta) (log2)'))+
  theme_gen+
  facet_grid(~CellType)


ggplot(df_20percent_both,
       aes(x=log2(MEAN), y=log2(alpha), shape=CellType, color = CellType))+
  geom_point(size=2, stroke=2, alpha=0.8)+
  geom_smooth(method=lm, colour= "grey60")+
  theme_minimal()+
  scale_color_manual(values= celltypecol2)+
  scale_shape_manual(values= c(1,2))+
  xlab(bquote('Average of copy number (log2)'))+
  ylab(bquote('Burst frequency (alpha) (log2)'))+
  theme_gen+
  facet_wrap(~CellType)


##spearman correlation of parameters and mean expression
cor.test(log2(df_20percent_both[df_20percent_both$CellType=='Pyramidal', 'beta']), 
         log2(df_20percent_both[df_20percent_both$CellType=='Pyramidal', 'MEAN']), method='spearman')
cor.test(log2(df_20percent_both[df_20percent_both$CellType=='Fast-spiking', 'beta']), 
         log2(df_20percent_both[df_20percent_both$CellType=='Fast-spiking', 'MEAN']), method='spearman')

cor.test(log2(df_20percent_both[df_20percent_both$CellType=='Pyramidal', 'alpha']), 
         log2(df_20percent_both[df_20percent_both$CellType=='Pyramidal', 'MEAN']), method='spearman')
cor.test(log2(df_20percent_both[df_20percent_both$CellType=='Fast-spiking', 'alpha']), 
         log2(df_20percent_both[df_20percent_both$CellType=='Fast-spiking', 'MEAN']), method='spearman')

##box plot for  on state
df_20percent_both<- df_20percent_both %>% mutate(bin_on_freq = ntile(on_state_freq, n=10))

ggplot(df_20percent_both, aes(as.factor(bin_on_freq), log2(alpha), fill=CellType))+
  geom_violin(linewidth=1.2, color='black', trim = TRUE)+
  scale_fill_manual(values= celltypecol3)+
  geom_boxplot(width=0.1,size=1,outlier.alpha = 0.5, color='black', 
               position = position_dodge(width=0.87))+
  theme_minimal()+
  theme_gen+
  ylab("Burst frequency (alpha) (log2)")+
  xlab("Percentage of on state cells")+
  facet_wrap(~CellType)

ggplot(df_20percent_both, aes(as.factor(bin_on_freq), log2(beta), fill=CellType))+
  geom_violin(linewidth=1.2, color='black', trim = TRUE)+
  geom_boxplot(width= 0.1,size=1,outlier.alpha = 0.5, color='black', 
               position = position_dodge(width=0.87))+
  scale_fill_manual(values= celltypecol3)+
  theme_minimal()+
  theme_gen+
  ylab("Burst size (beta) (log2)")+
  xlab("Percentage of on state cells")+
  facet_wrap(~CellType)



#correlation of copy number distribition measures 
df_ggpairs<- cbind(log2(df_20percent_both[, c(2:3,11,12,14)]), 
                   df_20percent_both[, c(7,10)])


ggpairs(df_ggpairs[df_ggpairs$CellType=='Pyramidal',],
        columns = 3:6,
        upper = list(continuous = wrap("cor", method = "spearman")),
        diag = list(continuous = "blankDiag"),
        columnLabels = c( 'Median copy number (log2)','Average copy number (log2)', 
                          'Copy number variance  (log2)',
                          '% ON state cells'),
        ggplot2::aes(shape=CellType, color=CellType))+
  theme_classic()+
  geom_smooth(aes(), method=lm, color= 'black')+
  scale_color_manual(values= celltypecol2)+
  scale_shape_manual(values= c(2,1))

ggpairs(df_ggpairs[df_ggpairs$CellType=='Fast-spiking',],
        columns = 3:6,
        upper = list(continuous = wrap("cor", method = "spearman")),
        diag = list(continuous = "blankDiag"),
        columnLabels = c( 'Median copy number (log2)','Average copy number (log2)', 
                          'Copy number variance  (log2)',
                          '% ON state cells'),
        ggplot2::aes(shape=CellType, color=CellType))+
  theme_classic()+
  geom_smooth(aes(), method=lm, color= 'black')+
  scale_color_manual(values= celltypecol2)+
  scale_shape_manual(values= c(2,1))

