


library(readxl)
library(ggplot2)


df_20percent_both<- read.csv("burst_kinetics_copynumber.csv")
df_20percent_both<- df_20percent_both[,-1]

pivot_df<- read.csv("burst_kinetics_copynumber_wider.csv")

theme_gen<- theme(axis.title.x = element_text(face="bold", size=28),
                  axis.text.x  = element_text(face= "bold", color = "black", size=26),
                  axis.title.y = element_text(face="bold", size=28),
                  axis.text.y  = element_text(face= "bold", color = "black", size=26),
                  axis.line = element_line(linewidth=1.5),
                  axis.ticks = element_line(linewidth = 0.8),
                  strip.text = element_text(face = "bold", size = rel(1.8)),
                  legend.title= element_blank(),
                  legend.text = element_text(face= "bold", color = "black", size=26))


###reshapeing data frames
pyr<- df_20percent_both[df_20percent_both$CellType == "Pyramidal",]
fs<- df_20percent_both[df_20percent_both$CellType == "Fast-spiking",]

delta_fspyr<- data.frame(Gene.name=fs$Gene.name, delta_alpha=pyr$alpha-fs$alpha, delta_beta=pyr$beta-fs$beta,
                         ratio_alpha=fs$alpha/pyr$alpha, ratio_beta = fs$beta/pyr$beta, ratio_MCpval=fs$MCpval/pyr$MCpval,
                         MEAN_fs= fs$MEAN, MEAN_pyr=pyr$MEAN, MEAN_ratio= fs$MEAN/pyr$MEAN,
                         MODE_fs= fs$MODE, MODE_pyr= pyr$MOD, MODE_ratio= fs$MODE/pyr$MODE,
                         on_state_freq_fs= fs$on_state_freq, on_state_freq_pyr=pyr$on_state_freq)

##standard deviation of alpha and beta parameters
sd_alpha <- sd(delta_fspyr$ratio_alpha)
sd_beta<- sd(delta_fspyr$ratio_beta)


ggplot(pivot_df, aes(x=log2(alpha_pyr), y=log2(alpha_fs)))+
  geom_point( size=2, stroke=2, alpha=0.7, color='grey20')+
  lapply(c(-log2(3*sd_alpha), log2(3*sd_alpha) ), function(o)
    geom_abline(slope=1, intercept=o, color='grey20', linewidth=1))+
  xlab(bquote('Burst frequency (alpha) Pyr (log2)'))+
  ylab(bquote('Burst frequency (alpha) FS (log2)'))+
  theme_classic()+
  theme_gen

ggplot(pivot_df, aes(x=log2(beta_pyr), y=log2(beta_fs)))+
  geom_point( size=2, stroke=2, alpha=0.7, color='grey20')+
  lapply(c(-log2(3*sd_beta), log2(3*sd_beta) ), function(o)
    geom_abline(slope=1, intercept=o, color='grey20', linewidth=1))+
  xlab(bquote('Burst size (beta) Pyr (log2)'))+
  ylab(bquote('Burst size (beta) FS (log2)'))+
  theme_classic()+
  theme_gen


fscol<- c("Larger in FS"="#66c2a5", "All" = "grey")
pyrcol<- c("All"="grey", "Larger in Pyr" = "#fc8d62")



ggplot(pivot_df, aes(x=log2(MEAN_pyr), y=log2(MEAN_fs)))+
  geom_point( size=2, stroke=2, alpha=0.7, color='grey20')+
  geom_point(data=pivot_df[pivot_df$direction_alpha == "Larger in FS" ,], 
             color='#66c2a5',  alpha=0.9,
             size=3, stroke=2)+
  geom_point(data=pivot_df[pivot_df$direction_alpha== 'Larger in Pyr',], 
             color='#fc8d62',  alpha=0.9,
             size=3, stroke=2)+
  lapply(c(log2(0.1), log2(10) ), function(o)
    geom_abline(slope=1, intercept=o, color='grey20', linewidth=1))+
  xlab(bquote('Average copy number Pyr (log2)'))+
  ylab(bquote('Average copy number FS (log2)'))+
  theme_classic()+
  theme_gen


ggplot(pivot_df, aes(x=log2(MEAN_pyr), y=log2(MEAN_fs)))+
  geom_point( size=2, stroke=2, alpha=0.7, color='grey20')+
  geom_point(data=pivot_df[pivot_df$direction_beta == "Larger in FS" ,], 
             color='#66c2a5',  alpha=0.9,
             size=3, stroke=2)+
  geom_point(data=pivot_df[pivot_df$direction_beta == 'Larger in Pyr',], 
             color='#fc8d62',  alpha=0.9,
             size=3, stroke=2)+
  lapply(c(log2(0.1), log2(10) ), function(o)
    geom_abline(slope=1, intercept=o, color='grey20', linewidth=1))+
  xlab(bquote('Average copy number Pyr (log2)'))+
  ylab(bquote('Average copy number FS (log2)'))+
  theme_classic()+
  theme_gen
