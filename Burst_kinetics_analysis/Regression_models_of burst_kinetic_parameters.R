
library(dplyr)
library(matrixStats)
library(readxl)
library(ggplot2)


df_20percent_both<- read.csv("burst_kinetics_copynumber.csv")
df_20percent_both<- df_20percent_both[,-1]

theme_gen<- theme(axis.title.x = element_text(face="bold", size=28),
                  axis.text.x  = element_text(face= "bold", color = "black", size=26),
                  axis.title.y = element_text(face="bold", size=28),
                  axis.text.y  = element_text(face= "bold", color = "black", size=26),
                  axis.line = element_line(linewidth=1.5),
                  axis.ticks = element_line(linewidth = 0.8),
                  strip.text = element_text(face = "bold", size = rel(1.8)),
                  legend.title= element_blank(),
                  legend.text = element_text(face= "bold", color = "black", size=26))


celltypecol2<- c("Fast-spiking"="black", "Pyramidal" = "grey40")
celltypecol3<- c("Fast-spiking"="white", "Pyramidal" = "grey40")

##regression model of alpha
lm_alpha<- lm(log2(alpha) ~ log2(MEAN)* on_state_freq * log2(MEDIAN), data=df_20percent_both[df_20percent_both$MEDIAN>2,])
df_20percent_both[df_20percent_both$MEDIAN>2,'alpha_fit']<- lm_alpha$fitted.values

alpha_model_corr<- cor.test(df_20percent_both[df_20percent_both$MEDIAN>2,'alpha_fit'],
                            df_20percent_both[df_20percent_both$MEDIAN>2,'alpha'], method = 'spearman')


ggplot(df_20percent_both[df_20percent_both$MEDIAN>2 & df_20percent_both$on_state>8
                         #(df_20percent_both$alpha > 0.5 & df_20percent_both$alpha_fit < -4) == FALSE
                         ,], 
       aes(x=log2(alpha), y=alpha_fit, shape=CellType, color=CellType))+
  geom_point( size=1, stroke=1.2, alpha=0.8)+
  geom_smooth(method="glm", color="#8DA0CB")+
  theme_minimal()+
  scale_color_manual(values= celltypecol2)+
  scale_shape_manual(values= c(1,2))+
  xlab(bquote('Burst frequency (alpha) (log2)'))+
  ylab(bquote('Burst frequency (model) '))+
  facet_grid(~CellType)+
  theme_gen

##regression model of beta

lm_beta<- lm(log2(beta) ~ log2(alpha) + on_state_freq* log2(MEAN)*log2(MEDIAN), data=df_20percent_both[df_20percent_both$MEDIAN>2,])

df_20percent_both[df_20percent_both$MEDIAN>2, 'beta_fit']<- lm_beta$fitted.values

beta_model_corr<- cor.test(df_20percent_both[df_20percent_both$MEDIAN>2,'beta_fit'],
                           df_20percent_both[df_20percent_both$MEDIAN>2,'beta'], method = 'spearman')

ggplot(df_20percent_both[df_20percent_both$MEDIAN>2,], 
       aes(x=log2(beta), y=beta_fit, shape=CellType, color=CellType))+
  geom_point( size=0.7, stroke=1.2, alpha=0.8)+
  geom_smooth(method="glm", color="#8DA0CB")+
  theme_minimal()+
  scale_color_manual(values= celltypecol2)+
  scale_shape_manual(values= c(1,2))+
  xlab(bquote('Burst size (beta) (log2)'))+
  ylab(bquote('Burst size (model) '))+
  facet_grid(~CellType)+
  theme_gen


df_ggpairs<- cbind(log2(df_20percent_both[, c(2:3,11,12,14)]), 
                   df_20percent_both[, c(7,10)])


corfs_mtx<- cor( df_ggpairs[df_ggpairs$CellType=='Fast-spiking', c(1:6)], method = 'spearman')
corpyr_mtx<- cor(df_ggpairs[df_ggpairs$CellType=='Pyramidal', c(1:6)], method = 'spearman')


label<- c('Burst frequency', 'Burst size','Median copy number', 
          'Average copy number', 'Copy number variance', '% ON state cells') 
library(gplots)
library(RColorBrewer)
library(GGally)
colfunc<- colorRampPalette(c("#fecc5c", "white","#74a9cf"))
#dev.off()

heatmap.2(corpyr_mtx, col = colfunc(128),
          trace = "none", density.info = "none",
          scale = "none",
          labRow = label,
          labCol = label,
          cexRow = 1,
          cexCol = 1,
          dendrogram ='none',
          Rowv = FALSE,
          Colv = FALSE,
          margins = c(5,5),
          key = TRUE,
          keysize = 0.0015, 
          cellnote=round(corpyr_mtx,2),
          notecol="grey20",
          lhei = c(10,10),
          lwid = c(10, 10)          
)

heatmap.2(corfs_mtx, col = colfunc(128),
          trace = "none", density.info = "none",
          scale = "none",
          labRow = label,
          labCol = label,
          cexRow = 1,
          cexCol = 1,
          dendrogram ='none',
          Rowv = FALSE,
          Colv = FALSE,
          margins = c(5,5),
          key = TRUE,
          keysize = 0.0015, 
          cellnote=round(corfs_mtx,2),
          notecol="grey20",
          lhei = c(10,10),
          lwid = c(10, 10)          
)
