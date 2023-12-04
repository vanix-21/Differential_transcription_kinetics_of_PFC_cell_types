
library(goseq)
library(forcats)
library(dplyr)
library(tidyr)

#import dATA
df_20percent_both <- read.csv("filtered20percent_sets.csv", row.names=1)

ccaat<- read.csv("CCAATbox.csv", sep=";")
tata<- read.csv("TATA.csv", sep=";")
gc<- read.csv("GCbox.csv", sep=";")

#merge data.frames
merged_promoters_bk<- merge(df_20percent_both, 
                            tata[tata$alternative_promoters==1,8:9], by= 'Gene.name', all.x = TRUE)
merged_promoters_bk<- merge(merged_promoters_bk, 
                            ccaat[ccaat$alternative_promoters==1,8:9], by= 'Gene.name', all.x = TRUE)
merged_promoters_bk<- merge(merged_promoters_bk, 
                            gc[gc$alternative_promoters==1,8:9], by= 'Gene.name', all.x = TRUE)

merged_promoters_bk$Group<- paste(merged_promoters_bk$CCAAT.box , merged_promoters_bk$Gcbox,
                                  merged_promoters_bk $TATA)
merged_promoters_bk[merged_promoters_bk$Group=='NA NA NA', 'Group']='No motif'

#longer data form
promoters_pivot_longer<- pivot_longer(merged_promoters_bk, cols = c('TATA', 'Gcbox', 'CCAAT.box') , values_to = 'promoter')
promoters_pivot_longer[promoters_pivot_longer$Group=='No motif', 'promoter']<- 'No motif'
promoters_pivot_longer$promoter<- factor(promoters_pivot_longer$promoter, levels = c('No motif',
                                                                                     'Gcbox', 
                                                                                     'CCAAT-box',
                                                                                     'TATA'))
#filter out repeated data
promoters_pivot_longer<-promoters_pivot_longer[is.na(promoters_pivot_longer$promoter)==FALSE, ] 

#cut data into to data.frames
promoters_longer_fs<- as.data.frame(promoters_pivot_longer[(promoters_pivot_longer$CellType == 'Fast-spiking'),])
promoters_longer_pyr<- as.data.frame(promoters_pivot_longer[(promoters_pivot_longer$CellType == 'Pyramidal'),])


##anova and tukey statistics for parameter-promoter dependence
anova(lm(promoters_pivot_longer$beta ~promoters_pivot_longer$promoter))
TukeyHSD(aov(beta ~ promoter, data=promoters_pivot_longer))

anova(lm(promoters_pivot_longer$alpha ~promoters_pivot_longer$promoter))
TukeyHSD(aov(alpha ~ promoter, data=promoters_pivot_longer))

#promoter frequency data
promoter_oprof <- read.csv("promoter_oprof.csv", sep=";")

subset<- c('AlphalargerinFS', 'BetalargerinPyr')
motif<- c('TATA', 'CCAAT', 'GC')
range<- list(c(-79.5:20.5), c(-199.5:20.5), c(-139.5:50.5))
prom_ks <-c()


#kolmogorov smirnoff test for distribution differences between all genes and subsets
i=1
j=1
for (i in 1:length(subset)){
  for (j in 1:length(motif)){
    prom_ks <-ks.test(promoter_oprof[promoter_oprof$genes==subset[i] & promoter_oprof$Motif==motif[j] & promoter_oprof$Position %in% range[[j]],'frequency.percent.'],
                      promoter_oprof[promoter_oprof$genes=='All' & promoter_oprof$Motif==motif[j]& promoter_oprof$Position %in% range[[j]] ,'frequency.percent.'])
    print(c(subset[i], motif[j], prom_ks$p.value ))

  }
}


