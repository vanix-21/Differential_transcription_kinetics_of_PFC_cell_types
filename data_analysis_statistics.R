library(dplyr)
library(matrixStats)
library(readxl)
library(ggplot2)
library(ggpubr)

#import tables
PyrBurstKinetics_5758_Genes <- read_excel("C:/Users/Vanda/singlecell/PredictedBurstKinetics_5758_Genes.xlsx", 
                                          sheet = "Pyramidal cells")
FSBurstKinetics_5758_Genes <- read_excel("C:/Users/Vanda/singlecell/PredictedBurstKinetics_5758_Genes.xlsx", 
                                         sheet = "Fast-spiking cells")

FS <- read_excel("C:/Users/Vanda/singlecell/RawData_scRNA-seq.xlsx", 
                 sheet = "FSC")
PYR <- read_excel("C:/Users/Vanda/singlecell/RawData_scRNA-seq.xlsx", 
                  sheet = "Pyramidal")

#setcelltypes
PyrBurstKinetics_5758_Genes$CellType<-c("Pyramidal") 
FSBurstKinetics_5758_Genes$CellType<-c("Fast-spiking") 

df_two_both<- rbind(PyrBurstKinetics_5758_Genes, FSBurstKinetics_5758_Genes)
genes_two<- unique(df_two_both$Gene.name)

#Function for filtering data base on the number of state cells
  #df: dataframe used
  #n: the minimum number of cells that are needed to express each gene
  #colname: column name for the number of on state cells
OnStateFilter<- function(df, n){
  #filter genes that are expressed in at least n cells
  df_onfilt<- df[df$on_state >= n,]
  #filter genes that are expressed in at least n cells and expressed in both cell types
  both<-intersect(df_onfilt[df_onfilt$CellType == "Pyramidal", "Gene.name"], 
                  df_onfilt[df_onfilt$CellType == "Fast-spiking", "Gene.name"])
  
  df_onfilt_both<- rbind(df_onfilt[df_onfilt$Gene.name %in% both & 
                                     df_onfilt$CellType == "Pyramidal" ,],
                         df_onfilt[df_onfilt$Gene.name %in%  both &
                                     df_onfilt$CellType == "Fast-spiking",])
}



df_20percent_both<- OnStateFilter(as.data.frame(df_two_both),  df_two_both$total_cells*0.2) #20% expression filter
small_pars<- df_20percent_both[df_20percent_both$alpha< 1e-03,] ##genes with small lpha parameter
df_20percent_both<- df_20percent_both[!(df_20percent_both$Gene.name %in% small_pars$Gene.name),]#filering out too small parameters
genes_20percent<- unique(df_20percent_both$Gene.name) ##gene list


#df_20percent_both<- read.csv("C:/Users/Vanda/singlecell/filtered20percent_sets.csv")
#df_20percent_both<- df_20percent_both[,2:15]
#scatter plot alpha and beta

##Calculation of mean copy numbers
  #df1: copy number data of celltype 1
  #df2: copy number data of celltype 2  
  #genelist: list of genes 

MeanCopy<- function(df1, df2, genelist){
  df1_set<- as.data.frame(df1[df1$Gene.name %in% genelist,])
  df2_set<- as.data.frame(df2[df1$Gene.name %in% genelist,])
  
  df1_mtx<- as.matrix(df1_set[,-1])
  rownames(df1_mtx)<- df1_set$Gene.name
  df2_mtx<- as.matrix(df2_set[,-1])
  rownames(df2_mtx)<- df2_set$Gene.name
  
  df1_mtx[df1_mtx == 0]<- NA
  Mean1<- rowMeans(df1_mtx, na.rm = TRUE)
  Var1<- rowVars(df1_mtx, na.rm = TRUE)
  df2_mtx[df2_mtx == 0]<- NA
  Mean2<- rowMeans(df2_mtx, na.rm = TRUE)
  Var2<- rowVars(df2_mtx, na.rm = TRUE)
  
  MEAN<- c(Mean1, Mean2)
}


##Calculation of copy number variance
  #df1: copy number data of celltype 1
  #df2: copy number data of celltype 2  
  #genelist: list of genes 

Variance<- function(df1, df2, genelist){
  df1_set<- as.data.frame(df1[df1$Gene.name %in% genelist,])
  df2_set<- as.data.frame(df2[df1$Gene.name %in% genelist,])
  
  df1_mtx<- as.matrix(df1_set[,-1])
  rownames(df1_mtx)<- df1_set$Gene.name
  df2_mtx<- as.matrix(df2_set[,-1])
  rownames(df2_mtx)<- df2_set$Gene.name
  
  df1_mtx[df1_mtx == 0]<- NA
  Var1<- rowVars(df1_mtx, na.rm = TRUE)
  df2_mtx[df2_mtx == 0]<- NA
  Var2<- rowVars(df2_mtx, na.rm = TRUE)
  
  VAR<- c(Var1, Var2)
}

# Create mode() function to calculate highest density copy number for each gene
mode <- function(x) {
  
  x = x[!is.na(x)]
  plot<- density(x)
  
  peak=which.max(plot$y)
  res=plot$x[peak]
}

##Calculation of highest density copy number
  #df1: copy number data of celltype 1
  #df2: copy number data of celltype 2  
  #genelist: list of genes 

ModeDataFrame<- function(df1, df2, genelist){
  df1_set<- as.data.frame(df1[df1$Gene.name %in% genelist,])
  df2_set<- as.data.frame(df2[df1$Gene.name %in% genelist,])
  
  df1_mtx<- as.matrix(df1_set[,-1])
  rownames(df1_mtx)<- df1_set$Gene.name
  df2_mtx<- as.matrix(df2_set[,-1])
  rownames(df2_mtx)<- df2_set$Gene.name
  
  df1_mtx[df1_mtx == 0]<- NA
  df2_mtx[df2_mtx == 0]<- NA
  
  Mode1<- apply(df1_mtx, 1, mode)
  Mode2<- apply(df2_mtx, 1, mode)
  
  MODE<- c(Mode1, Mode2)
}


df_20percent_both$MEAN<- MeanCopy(PYR,  FS, genes_20percent)
df_20percent_both$MODE<- ModeDataFrame(PYR,  FS, genes_20percent)
df_20percent_both$VAR<- Variance(PYR,  FS, genes_20percent)
#write.csv(df_20percent_both, file = "filtered20percent_sets.csv")

##correlation of mean and highest density copy number
cor.test(log2(df_20percent_both[df_20percent_both$CellType=='Pyramidal', 'MEAN']), 
         log2(df_20percent_both[df_20percent_both$CellType=='Pyramidal', 'MODE']), method='spearman')
cor.test(log2(df_20percent_both[df_20percent_both$CellType=='Fast-spiking', 'MEAN']), 
         log2(df_20percent_both[df_20percent_both$CellType=='Fast-spiking', 'MODE']), method='spearman')


#correlation of mean copy number with parameters
cor.test(log2(df_20percent_both[df_20percent_both$CellType=='Pyramidal', 'beta']), 
         log2(df_20percent_both[df_20percent_both$CellType=='Pyramidal', 'MEAN']), method='spearman')
cor.test(log2(df_20percent_both[df_20percent_both$CellType=='Fast-spiking', 'beta']), 
         log2(df_20percent_both[df_20percent_both$CellType=='Fast-spiking', 'MEAN']), method='spearman')

cor.test(log2(df_20percent_both[df_20percent_both$CellType=='Pyramidal', 'alpha']), 
         log2(df_20percent_both[df_20percent_both$CellType=='Pyramidal', 'MEAN']), method='spearman')
cor.test(log2(df_20percent_both[df_20percent_both$CellType=='Fast-spiking', 'alpha']), 
         log2(df_20percent_both[df_20percent_both$CellType=='Fast-spiking', 'MEAN']), method='spearman')



##data cutting based on on state frequencies
df_20percent_both<- df_20percent_both %>% mutate(bin_on_freq = ntile(on_state_freq, n=10))


## wilcoxon-test for parameter distributions between cell types
pyr<- df_20percent_both[df_20percent_both$CellType == "Pyramidal",]
fs<- df_20percent_both[df_20percent_both$CellType == "Fast-spiking",]

wilcox.test(fs[, "alpha"], pyr[, "alpha"])
wilcox.test(fs[, "beta"], pyr[, "beta"])



##wider data frame
pivot_df<-cbind(pyr[,c(1:14)], fs[,c(2:14)])
colnames(pivot_df)<- c("Gene.name",  "alpha_pyr" ,"beta_pyr"  , "lambda1_pyr"  ,"lambda2_pyr"  ,
                       "MCpval_pyr"    ,"on_state_freq_pyr","on_state_pyr",  "total_cells_pyr" ,
                       "CellType_pyr","MEAN_pyr", "MODE_pyr", 'VAR_pyr', 'bin_on_freq_pyr',
                       "alpha_fs",  "beta_fs" ,       "lambda1_fs", "lambda2_fs"  , 
                       "MCpval_fs", "on_state_freq_fs",  "on_state_fs" , "total_cells_fs",  
                       "CellType_fs","MEAN_fs", "MODE_fs", 'VAR_fs','bin_on_freq_fs'
)


#ratio fs/pyr
pivot_df$MODE_ratio<- pivot_df$MODE_fs/pivot_df$MODE_pyr                  
pivot_df$alpha_ratio<- pivot_df$alpha_fs/pivot_df$alpha_pyr
pivot_df$beta_ratio<- pivot_df$beta_fs/pivot_df$beta_pyr
pivot_df$MEAN_ratio<- pivot_df$MEAN_fs/pivot_df$MEAN_pyr 



delta_fspyr<- data.frame(Gene.name=fs$Gene.name, delta_alpha=pyr$alpha-fs$alpha, delta_beta=pyr$beta-fs$beta,
                         ratio_alpha=fs$alpha/pyr$alpha, ratio_beta = fs$beta/pyr$beta, ratio_MCpval=fs$MCpval/pyr$MCpval,
                         MEAN_fs= fs$MEAN, MEAN_pyr=pyr$MEAN, MEAN_ratio= fs$MEAN/pyr$MEAN,
                         MODE_fs= fs$MODE, MODE_pyr= pyr$MOD, MODE_ratio= fs$MODE/pyr$MODE,
                         on_state_freq_fs= fs$on_state_freq, on_state_freq_pyr=pyr$on_state_freq)

#standard dev of parameter ratios
sd_alpha <- sd(delta_fspyr$ratio_alpha)
sd_beta<- sd(delta_fspyr$ratio_beta)



#Function for identifying genes of interest
  # par: either "alpha" or "beta"
  # pyr: data.frame of pyr burst kinetic data 
  # fs : data.frame of fs burst kinetic data 
  # delta_fspyr : data.frame of burst kinetic differences
GoiRatio<- function(par, pyr, fs, delta_fspyr){
  
  goi_all<- data.frame(no=c(1:2630, 1:2630), Gene.name= c(pyr$Gene.name, fs$Gene.name),
                       ratio_alpha=c(delta_fspyr$ratio_alpha, delta_fspyr$ratio_alpha),
                       ratio_beta=c(delta_fspyr$ratio_beta, delta_fspyr$ratio_beta),
                       MEAN=  c(delta_fspyr$MEAN_pyr, delta_fspyr$MEAN_fs),
                       MODE= c(delta_fspyr$MODE_pyr, delta_fspyr$MODE_fs),
                       CellType=c(pyr$CellType, fs$CellType),
                       on_state_freq= c(delta_fspyr$on_state_freq_pyr, delta_fspyr$on_state_freq_fs))
  if (par=='alpha'){
    goi_up<- goi_all[goi_all$ratio_alpha > sd_alpha*3,]
    
    goi_down<-goi_all[goi_all$ratio_alpha < 1/(sd_alpha*3),]
  }
  
  if (par=='beta'){
    goi_up<- goi_all[goi_all$ratio_beta > sd_beta*3, ]
    
    goi_down<- goi_all[goi_all$ratio_beta < 1/(sd_beta*3),]
  }
  
  goi_down$direction<- "Larger in Pyr"
  goi_up$direction<- "Larger in FS"
  
  goi<- rbind(goi_down, goi_up)
  return(goi)
}


#genes of interest list based on parameter differences
goi_alpha<- GoiRatio("alpha", pyr, fs, delta_fspyr)
goi_beta<- GoiRatio("beta", pyr, fs, delta_fspyr)


##adding labels for differential burst kinetics
pivot_df$direction_alpha<- c('')
pivot_df$direction_beta<- c('')

i=1
for (i in c(1:2630)){
  if (pivot_df[i,'alpha_ratio'] > sd_alpha*3){
    pivot_df[i,'direction_alpha'] <- 'Larger in FS'
  }
  if (pivot_df[i,'alpha_ratio']< 1/(sd_alpha*3)){
    pivot_df[i,'direction_alpha'] <- 'Larger in Pyr'
  }
}

i=1
for (i in  c(1:2630)){
  if (pivot_df[i,'beta_ratio'] > sd_beta*3){
    pivot_df[i,'direction_beta'] <- 'Larger in FS'
  }
  if (pivot_df[i,'beta_ratio']< 1/(sd_beta*3)){
    pivot_df[i,'direction_beta'] <- 'Larger in Pyr'
  }
}

#write.csv(pivot_df, file = 'filtered_20percent_wider.csv')