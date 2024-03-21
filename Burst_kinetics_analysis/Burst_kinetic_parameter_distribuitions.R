
library(readxl)
library(ggplot2)


df_20percent_both<- read.csv("burst_kinetics_copynumber.csv")
df_20percent_both<- df_20percent_both[,-1]



celltypecol2<- c("Fast-spiking"="black", "Pyramidal" = "grey40")


theme_gen<- theme(axis.title.x = element_text(face="bold", size=28),
                  axis.text.x  = element_text(face= "bold", color = "black", size=26),
                  axis.title.y = element_text(face="bold", size=28),
                  axis.text.y  = element_text(face= "bold", color = "black", size=26),
                  axis.line = element_line(linewidth=1.5),
                  axis.ticks = element_line(linewidth = 0.8),
                  strip.text = element_text(face = "bold", size = rel(1.8)),
                  legend.title= element_blank(),
                  legend.text = element_text(face= "bold", color = "black", size=26))

#scatter plot alpha and beta
ggplot(df_20percent_both, aes(x=log2(alpha), y=log2(beta), col= CellType, shape= CellType))+
  geom_point(size=2, stroke=2, alpha=0.8)+
  theme_minimal()+
  scale_colour_manual(values= celltypecol2)+
  scale_shape_manual(values= c(1,2))+
  xlab(bquote('Burst frequency (alpha) (log2)'))+
  ylab(bquote('Burst size (beta) (log2)'))+
  theme_gen+
  facet_wrap(~CellType)

#overlapped
ggplot(df_20percent_both, aes(x=log2(alpha), y=log2(beta), color= CellType, shape= CellType))+
  geom_point( size=2, stroke=2, alpha=0.8)+
  theme_minimal()+
  scale_color_manual(values= celltypecol2)+
  scale_shape_manual(values= c(1,2))+
  xlab(bquote('Burst frequency (alpha) (log2)'))+
  ylab(bquote('Burst size (beta) (log2)'))+
  theme_gen

###Density plot of burst kinetics ditribution
ggplot(df_20percent_both, aes(x=log2(beta), col=CellType))+
  geom_density(alpha=0.5, linewidth=2 )+
  scale_colour_manual(values= celltypecol2)+
  ylab(bquote('Density (n/N)'))+
  xlab(bquote('Burst size (beta) (log2)'))+
  theme_minimal() +
  theme_gen

ggplot(df_20percent_both, aes(x=log2(alpha), col=CellType))+
  geom_density(alpha=0.5, linewidth=2 )+
  scale_colour_manual(values= celltypecol2)+
  ylab(bquote('Density (n/N)'))+
  xlab(bquote('Burst frequency (alpha) (log2)'))+
  theme_minimal() +
  theme_gen

##wilcoxon rank sum test for burst kinetic parameter differences between cell types
wilcox.test(df_20percent_both[df_20percent_both$CellType=='Fast-spiking', "alpha"], 
            df_20percent_both[df_20percent_both$CellType=='Pyramidal', "alpha"])

wilcox.test(df_20percent_both[df_20percent_both$CellType=='Fast-spiking', "beta"], 
            df_20percent_both[df_20percent_both$CellType=='Pyramidal', "beta"])


