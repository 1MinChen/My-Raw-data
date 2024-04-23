if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
}
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
source('mg_base.R')
my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+   
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                
          plot.title = element_text(hjust = 0.5),text = element_text(family = 'Times',size = size),
          axis.text.x = element_text(angle = angle, hjust = hjust),
          panel.grid = element_blank()) 
  return(p)
}


#GPL18109######
# ###fasta
# gpl <- data.table::fread("GPL18109_family/GPL18109-tbl-1.txt",data.table = F)
# dim(gpl)
# head(gpl)
# gpl[1:15,]
# ## 
# gpl <- gpl[,c(1,6)]
# colnames(gpl)=c('ID','SEQUENCE')
# ## 
# library(dplyr)
# gpl <- gpl %>% 
#   filter(nchar(SEQUENCE)!=0)  
# head(gpl)
# dim(gpl)
# ## 
# gp <- paste0('>',gpl$ID,'\n', gpl$SEQUENCE)
# 
# ## 
# write.table(gp,'GPL18109_family/GPL.fasta', quote = F, row.names = F, col.names = F)


##
## 
probe2ID <- data.table::fread("seqmap_results.txt",data.table = F)
## 
library(tidyr)
library(dplyr)
colnames(probe2ID)
# "trans_id"     "trans_coord"  "target_seq"   "probe_id"     "probe_seq"    "num_mismatch" "strand"  
# str_split_fixed(probe2ID$trans_id,'\\|',10)[,5]
probe2ID <- probe2ID %>%
  select(probe_id,trans_id) %>% 
  separate(trans_id,into = c("Ensembl",
                             "drop1","drop2","drop3",
                             "trans_Symble","gene_Symble","drop4","trans_biotype"),sep = "\\|") %>% 
  select(probe_id,Ensembl,trans_Symble,gene_Symble,trans_biotype)


head(probe2ID)
table(probe2ID$trans_biotype)
length(probe2ID$probe_id)
length(unique(probe2ID$probe_id))



#mRNA和lncRNA
table(probe2ID$trans_biotype)
mrna_symbol=unique(probe2ID$gene_Symble[probe2ID$trans_biotype=='protein_coding'])
lncrna_symbol=unique(probe2ID$gene_Symble[probe2ID$trans_biotype=='lncRNA'])
length(mrna_symbol)
length(lncrna_symbol)


#GSE53625####
load('origin_datas/GEO/GSE53625.RData')
GSE53625.pheno=pData(GSE53625)
head(GSE53625.pheno)
GSE53625.pheno=data.frame(Samples=GSE53625.pheno$geo_accession,
                          tissue=str_split_fixed(GSE53625.pheno$`tissue:ch1`,' ',2)[,1],
                          OS.time=GSE53625.pheno$`survival time(months):ch1`,
                          OS=GSE53625.pheno$`death at fu:ch1`)
rownames(GSE53625.pheno)=GSE53625.pheno$Samples
head(GSE53625.pheno)
table(GSE53625.pheno$tissue)
GSE53625.pheno$cohort=rep(c('GSE53622','GSE53624'),c(120,238))
table(GSE53625.pheno$cohort)
# 
# GSE53625.exp=exprs(GSE53625)
# GSE53625.exp[1:5,1:5]
# rownames(GSE53625.exp)
# dim(GSE53625.exp)
# 
# #
# head(probe2ID)
# GSE53625.exp=exp_probe2symbol_v2(datExpr = GSE53625.exp,anno = probe2ID[,c(1,4)],method = 'mean')
# save(GSE53625.exp,file = 'origin_datas/GEO/GSE53625.exp.RData')
load('origin_datas/GEO/GSE53625.exp.RData')
GSE53625.exp[1:5,1:5]
range(GSE53625.exp)
dim(GSE53625.exp)
boxplot(GSE53625.exp)

# ###
# library(ggbiplot)
# before.pca <- prcomp(t(GSE53625.exp[,GSE53625.pheno$Samples]), scale=T)
# pca_before <- ggbiplot(before.pca, scale=1, groups = GSE53625.pheno$cohort,
#                        ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
#   scale_color_manual(values = ggsci::pal_lancet('lanonc')(9)) + 
#   # xlim(-5, 5) + ylim(-5,5) +
#   theme_light() +
#   theme(legend.direction = 'horizontal', legend.position = 'top') +
#   xlab('PCA1') + ylab('PCA2')
# pca_before
# 
# library(limma)
# library(sva)
# GSE53625.exp1 = ComBat(GSE53625.exp[,GSE53625.pheno$Samples],batch =  GSE53625.pheno$cohort)
# after.pca <- prcomp(t(GSE53625.exp1[,GSE53625.pheno$Samples]), scale=T)
# pca_after <- ggbiplot(after.pca, scale=1, groups = GSE53625.pheno$cohort,
#                        ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
#   scale_color_manual(values = ggsci::pal_lancet('lanonc')(9)) + 
#   # xlim(-5, 5) + ylim(-5,5) +
#   theme_light() +
#   theme(legend.direction = 'horizontal', legend.position = 'top') +
#   xlab('PCA1') + ylab('PCA2')
# pca_after
# boxplot(GSE53625.exp1)

# save(GSE53625.exp,file='GSE53625.exp.RData')

#GSE43732[miRNA]####
load('origin_datas/GEO/GSE43732.RData')
GSE43732.pheno=pData(GSE43732)
head(GSE43732.pheno)
table(GSE43732.pheno$`tissue type:ch1`)
GSE43732.pheno=data.frame(Samples=GSE43732.pheno$geo_accession,tissue=GSE43732.pheno$`tissue type:ch1`)

GSE43732.df=exprs(GSE43732)
dim(GSE43732.df)
GSE43732.df[1:5,1:5]
dim(GSE43732.df)
boxplot(GSE43732.df)



#01.#####
dir.create('results/01.DEmRNA')
##1.1 ####
pcd.genesets=readxl::read_excel('origin_datas/PCD.geneSets.PMID36341760.xlsx')
pcd.genesets=data.frame(pcd.genesets)
pcd.genesets=pcd.genesets[,-1]
head(pcd.genesets)
pcd.genesets.list=list()
pcd.genesets.df=c()
for(i in colnames(pcd.genesets)){
  pcd.genesets.list[[i]]=as.character(na.omit(pcd.genesets[,i]))
  pcd.genesets.df=rbind(pcd.genesets.df,data.frame(PCD=i,Symbol=pcd.genesets[,i],check.names = F,stringsAsFactors = F))
}
pcd.genes=unique(pcd.genesets.df$Symbol)
length(pcd.genes)
#1255

##1.2mRNA####
GSE53625.limma=mg_limma_DEG(exp = GSE53625.exp[rownames(GSE53625.exp)%in%mrna_symbol,GSE53625.pheno$Samples],
                            group = GSE53625.pheno$tissue,
                            ulab = 'cancer',dlab = 'normal')
GSE53625.limma$Summary
GSE53625.DEmRNA=GSE53625.limma$DEG[GSE53625.limma$DEG$adj.P.Val<0.05 & abs(GSE53625.limma$DEG$logFC)>1,]
head(GSE53625.DEmRNA)
write.csv(GSE53625.DEmRNA,'results/01.DEmRNA/GSE53625.DEmRNA.csv')

my_volcano=function(dat,p_cutoff=0.05,fc_cutoff=1,col=c("red","blue","black"),
                    ylab='-log10 (adj.PVal)',xlab='log2 (FoldChange)',leg.pos='right'){
  degs_dat=dat$DEG
  degs_dat$type=factor(ifelse(degs_dat$adj.P.Val<p_cutoff & abs(degs_dat$logFC) > fc_cutoff, 
                              ifelse(degs_dat$logFC> fc_cutoff ,'Up','Down'),'No Signif'),levels=c('Up','Down','No Signif'))
  p=ggplot(degs_dat,aes(x=logFC,y=-log10(adj.P.Val),color=type))+
    geom_point()+
    scale_color_manual(values=col)+#
    # geom_text_repel(
    #   data = tcga.diff$DEG[tcga.diff$DEG$adj.P.Val<p_fit & abs(tcga.diff$DEG$logFC)>fc_fit,],
    #   #aes(label = Gene),
    #   size = 3,
    #   segment.color = "black", show.legend = FALSE )+#
    theme_bw()+#
    theme(
      legend.title = element_blank(),#
      legend.position = leg.pos,
      text = element_text(family = 'Times')
    )+
    ylab(ylab)+#
    xlab(xlab)+#
    geom_vline(xintercept=c(-fc_cutoff,fc_cutoff),lty=3,col="black",lwd=0.5) +#|FoldChange|>2
    geom_hline(yintercept = -log10(p_cutoff),lty=3,col="black",lwd=0.5)#padj<0.05
  return(p)
}
my_volcano(dat = GSE53625.limma,p_cutoff = 0.05,fc_cutoff = 1)+ggtitle('GSE53625 DEmRNA')
ggsave('results/01.DEmRNA/degs_volcano.pdf',height = 6,width = 7.5)



##1.3 PCD####
PCD.degs=intersect(pcd.genes,rownames(GSE53625.DEmRNA))
length(PCD.degs)
#149
PCD.degs.FC=GSE53625.DEmRNA[PCD.degs,]
head(PCD.degs.FC)
PCD.degs.FC$type=ifelse(PCD.degs.FC$logFC>0,'up','down')
table(PCD.degs.FC$type)
write.csv(PCD.degs.FC[,-7],'results/PCD.degs.csv')

# 
test <- PCD.degs.FC[
  with(PCD.degs.FC, order(type, logFC)),
]
nrow(test)
test=test[c(1:25,125:149),]
pdf('results/degs_heatmap2.pdf',height = 6,width = 8,onefile = F)
Heatmap(as.matrix(t(scale(t(GSE53625.exp[rownames(test),GSE53625.pheno$Samples[order(GSE53625.pheno$tissue)]]))))
        , name = "Expr", row_split = test$type
        , column_split = GSE53625.pheno$tissue[order(GSE53625.pheno$tissue)]
        , row_title_gp = gpar(fill =ggsci::pal_nejm()(10)[4:5])
        , column_title_gp = gpar(fill =ggsci::pal_nejm()(10))
        , cluster_rows = T, cluster_columns = F
        , cluster_row_slices = F, cluster_column_slices=T
        , show_row_dend = F, show_column_dend = F
        , show_row_names = T, show_column_names = F
        , col = circlize::colorRamp2(c(-5, 0, 5), c('#3B4992FF', 'white', '#EE0000FF')))
dev.off()

##1.4 ######
PCD.degs.enrich=mg_clusterProfiler(genes =PCD.degs)
p1=enrichplot::dotplot(PCD.degs.enrich$KEGG)+theme(text = element_text(family = 'Times',size = 16))+ggtitle('KEGG')
p2=enrichplot::dotplot(PCD.degs.enrich$GO_BP)+theme(text = element_text(family = 'Times',size = 16))+ggtitle('Biological Process')
p3=enrichplot::dotplot(PCD.degs.enrich$GO_CC)+theme(text = element_text(family = 'Times',size = 16))+ggtitle('Cell Component')
p4=enrichplot::dotplot(PCD.degs.enrich$GO_MF)+theme(text = element_text(family = 'Times',size = 16))+ggtitle('Molecular Function')
p=mg_merge_plot(p1,p2,p3,p4,ncol=2,nrow=2,labels = LETTERS[3:6])
ggsave('results/01.DEmRNA/DE_PCDs_enrichment.pdf',p,height = 12,width = 15)

write.xlsx(list(KEGG=PCD.degs.enrich$KEGG@result,GO_BP=PCD.degs.enrich$GO_BP@result,
                GO_CC=PCD.degs.enrich$GO_CC@result,GO_MF=PCD.degs.enrich$GO_MF@result),
            'results/01.DEmRNA/DE_PCDs_enrichment.xlsx',overwrite = T)


#02.####
ML_dat=t(GSE53625.exp[PCD.degs,])
ML_dat=cbind.data.frame(type=GSE53625.pheno$tissue,ML_dat[GSE53625.pheno$Samples,])
ML_dat$type=as.factor(ML_dat$type)
head(ML_dat)
mat=as.matrix(ML_dat[,PCD.degs])
group=as.factor(ML_dat$type)

##2.1 （SVM）####
dir.create('results/02.SVM')
library(caret)
set.seed(123)
rfeControl = rfeControl(functions = caretFuncs, method ="cv", number= 10, verbose = FALSE)
subsetSizes = 1:50
rf1 = rfe(x = mat,y = group,sizes=subsetSizes,rfeControl = rfeControl, method ="svmLinear")
result_svm = rf1$result
write_tsv(result_svm,"results/02.SVM/SVM_feature_number.txt")
write_tsv(as.data.frame(rf1$optVariables),"results/02.SVM/SVM_selected_features.txt")

pdf('results/02.SVM/SVM_feature_select.pdf',height = 5,width = 7,onefile = F)
plot(result_svm$Variables,result_svm$Accuracy,type="l",xlim=c(1,30),col="blue",
     xlab="Number of features",ylab="10x CV accuracy")
points(9, 0.9943651, cex = 2, pch = 1, col ="red")
text(10,0.99,"9 - 0.9943651",col="red",cex=1)
dev.off()


##2.2 #####
dir.create('results/03.randomForest')
library(randomForest)

# 1.
n = ncol(as.matrix(mat))     # n=ncol(train_data_rf)
rate=1     # 

for(i in 1:(n-1)){
  set.seed(5)
  rf_train = randomForest(mat,group,mtry=i,ntree=1000)
  rate[i] = mean(rf_train$err.rate)   # 
  print(rf_train)    
}

rate # 

pdf("results/03.randomForest/RandomForest_mtry.pdf",width = 6,height = 6)
plot(rate)
points(which(rate==min(rate)),min(rate),cex = 2, pch = 3, col ="red")
mtry = 24
# points(mtry,rate[mtry],cex = 1, pch = 2, col ="red")
dev.off()


# 2.
set.seed(5)
rf_train = randomForest(mat,group,mtry=mtry,ntree=500)
write_tsv(as.data.frame(rf_train$err.rate),"results/03.randomForest/RandomForest_ntree.txt")

pdf("results/03.randomForest/RandomForest_ntree.pdf",width = 6,height = 6)
plot(rf_train)    #
dev.off()

ntree = 50
# 3.
set.seed(5)
rf_train = randomForest(mat,group,mtry=mtry,ntree=ntree,importance=TRUE,proximity=TRUE)    
rf_train
##
importance  = rf_train$importance
head(importance)
write.table(importance, 'results/03.randomForest/importance.txt', sep = '\t', col.names = NA, quote = FALSE)

# 
pdf("results/03.randomForest/RandomForest_features_top30.pdf",width = 10,height = 7)
varImpPlot(rf_train, n.var = min(30, nrow(rf_train$importance)), main = 'Top 30 - variable importance')
dev.off()



##2.3 lass###############
dir.create('results/04.LASSO')
#lasso
library(glmnet)
set.seed(321)
fit1=glmnet(x = mat,y = group,family = "binomial",nlambda=100, alpha=1) 
cv.fit<-cv.glmnet(x = mat,y = group,family = "binomial",nlambda=100, alpha=1)
sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
cv.fit$lambda.min
#0.0002895798

pdf('results/04.LASSO/LASSO.pdf',height = 6,width = 12,onefile = F)
par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()

coefficient = coef(cv.fit,s=cv.fit$lambda.min)
Active.Index = which(as.numeric(coefficient)!=0)
active.coefficients = as.numeric(coefficient)[Active.Index]
sig_gene_multi_cox = rownames(coefficient)[Active.Index]
sig_gene_multi_cox = sig_gene_multi_cox[-1]
active.coefficients=active.coefficients[-1]
write.table(data.frame(gene=sig_gene_multi_cox,coef=active.coefficients), 'results/04.LASSO/LASSO_result.txt', 
            sep = '\t', col.names = NA, quote = FALSE)


#03.hub #####
dir.create('results/05.Diagnostic Signatures')
RF = read_tsv("results/03.randomForest/importance.txt")
RF = RF %>% dplyr::arrange(desc(MeanDecreaseAccuracy))%>% head(n=30) %>% dplyr::select(`...1`) %>% unlist() 
write.table(as.character(RF),'results/03.randomForest/importance_TOP30.txt',quote = F,sep = '\t')

SVM=read_tsv('results/02.SVM/SVM_selected_features.txt')
SVM=SVM$`rf1$optVariables`

hub.genes=Reduce(intersect,list(sig_gene_multi_cox,SVM,RF))
hub.genes
library(eulerr)
v=list(LASSO=sig_gene_multi_cox,SVM=SVM,randomForest=RF)
venn.plot=plot(venn(v),labels = list(col = "gray20", font = 2), 
               edges = list(col="gray60", lex=1),
               fills = list(fill = c("#297CA0", "#E9EA77"), alpha = 0.6),
               quantities = list(cex=.8, col='gray20'))
venn.plot
ggsave('results/05.Diagnostic Signatures/venn.pdf',venn.plot,height = 5,width = 5)


#
library(pROC)
array.roc <- list()
for (i in 1:length(hub.genes)){
  roc1 <- roc(ML_dat$type, ML_dat[,hub.genes[i]])
  array.roc[[i]]=roc1
  names(array.roc)[i] <- paste0(hub.genes[i],' AUC=',round(roc1$auc[1],2))
}
ggroc(array.roc)+
  geom_segment(aes(x = 1, y = 0, xend =0, yend = 1), color="darkgrey", linetype="dashed")+
  theme_bw()+theme(panel.grid = element_blank(),legend.position = c(0.85,0.35))+
  ggtitle(label = 'GSE53625',subtitle = 'ROC curve of gene expression of Hub genes')+
  scale_colour_discrete(name="Gene",labels =names(array.roc))

library(e1071)
paste0(hub.genes,collapse = '+')
mylog_array=svm(type ~BOK+INHBA+LRRK2+HSP90AA1+HSPB8,data = ML_dat)
summary(mylog_array)
mylog_array<-predict(mylog_array)
mylog_array <- as.ordered(mylog_array)
modelroc_array <- roc(ML_dat$type,mylog_array)
modelroc_array
plot(modelroc_array, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.5, 0.2), grid.col=c("green", "red"), 
     max.auc.polygon=T, auc.polygon.col="skyblue", print.thres=TRUE)

pdf('results/05.Diagnostic Signatures/GSE53625_roc.pdf',height =6,width = 9)
par(mfrow=c(2,3))
plot(modelroc_array, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.5, 0.2), grid.col=c("green", "red"), 
     max.auc.polygon=T, auc.polygon.col="skyblue", print.thres=TRUE)
for (i in hub.genes){
  roc1 <- plot.roc(ML_dat$type, 
                   as.numeric(ML_dat[,i]),
                   percent=TRUE, col="2")
  title(main = paste0(i," AUC=",round(roc1$auc[1],2),'%'))
}
dev.off()


####hubgene
mg_PlotMutiBoxplot(data = t(GSE53625.exp[hub.genes,GSE53625.pheno$Samples]),group = GSE53625.pheno$tissue,
                   test_method = 'wilcox.test',add = 'boxplot',legend.pos = 'top')

data_ex = cbind.data.frame(t(GSE53625.exp1[hub.genes,GSE53625.pheno$Samples]),tissue=GSE53625.pheno$tissue)
library(ggpubr)
p=list()
for(i in 1:length(hub.genes)){
  dt<-data_ex[,c("tissue",hub.genes[i])]
  p[[i]]<-ggboxplot(dt,x='tissue',y=colnames(dt)[2],fill = "tissue", 
               palette = c("#00AFBB",  "#FC4E07"),main=hub.genes[i],ylab="Expression",xlab="")+
    theme(plot.title = element_text(hjust = 0.5),legend.position = "none",text = element_text(family = 'Times'))+
    stat_compare_means(method = "wilcox.test",label.x=1.5,label="p.signif")
  assign(paste("p", i, sep=""), p)
}
pdf('results/05.Diagnostic Signatures/hubgene_expr.pdf',height = 5,width = 15)
mg_merge_plot(p,ncol=5,labels = 'A')
dev.off()

#04.#####
dir.create('results/06.network')
# ##4.1 
# GSE53625.limma.lnc=mg_limma_DEG(exp = GSE53625.exp[rownames(GSE53625.exp)%in%lncrna_symbol,GSE53625.pheno$Samples],
#                                 group = GSE53625.pheno$tissue,
#                                 ulab = 'cancer',dlab = 'normal')
# GSE53625.limma.lnc$Summary
# GSE53625.DElncRNA=GSE53625.limma.lnc$DEG[GSE53625.limma.lnc$DEG$adj.P.Val<0.05 & abs(GSE53625.limma.lnc$DEG$logFC)>1,]
# head(GSE53625.DElncRNA)
# my_volcano(dat = GSE53625.limma.lnc,p_cutoff = 0.05,fc_cutoff = 1)

##4.1 mRNA-miRNA####
#
GSE43732.limma=mg_limma_DEG(exp = GSE43732.df[,GSE43732.pheno$Samples],group = GSE43732.pheno$tissue,
                            ulab = 'Cancer tissue',dlab = 'Adjacent normal tissue')
GSE43732.limma$Summary
GSE43732.DEmiRNA=GSE43732.limma$DEG[GSE43732.limma$DEG$adj.P.Val<0.05 & abs(GSE43732.limma$DEG$logFC)>1,]
GSE43732.DEmiRNA=GSE43732.DEmiRNA %>% drop_na(logFC)
head(GSE43732.DEmiRNA)
dim(GSE43732.DEmiRNA)
my_volcano(dat = GSE43732.limma,p_cutoff = 0.05,fc_cutoff = 1)+ggtitle('GSE43732 DEmiRNA')

## mRNA-miRNA
mRNA_miRNA =read.delim("results/06.network/Conserved_Family_Info.txt",check.names = F)
head(mRNA_miRNA)
table(mRNA_miRNA$gene_class)
mRNA_miRNA_res = mRNA_miRNA %>% filter(`Gene Symbol` %in% hub.genes)  %>%  filter(`miR Family` %in% gsub('hsa-','',rownames(GSE43732.DEmiRNA))) %>% filter(`Species ID` %in% 9606)  %>%
  dplyr::select(`miR Family`,`Gene Symbol`) %>% distinct()
head(mRNA_miRNA_res)
colnames(mRNA_miRNA_res)=c('miRNA','mRNA')
dim(mRNA_miRNA_res)


##4.2 mRNA-tf####
#从 HTFtarget 
target_tfs=read.delim('results/06.network/TF-Target-information.txt',check.names = F)
head(target_tfs)
mRNA_TFs = target_tfs %>% filter(`target` %in% hub.genes) %>% dplyr::select(TF,target) %>% distinct()
head(mRNA_TFs)
colnames(mRNA_TFs)=c('TF','mRNA')
dim(mRNA_TFs)




##4.3 ####
miRNA = intersect(unique(mRNA_miRNA_res$miRNA),gsub('hsa-','',rownames(GSE43732.DEmiRNA)))
mRNA = hub.genes
TF = unique(mRNA_TFs$TF)

nodes = tibble(gene=c(miRNA,mRNA,TF))
nodes$Type = ifelse(nodes$gene %in% miRNA,"miRNA",ifelse(nodes$gene %in% mRNA,"mRNA","TF"))
write.table(nodes,'results/06.network/nodes_anno.txt',quote = F,row.names = F,sep = '\t')

colnames(mRNA_miRNA_res) = c("TF/miRNA","mRNA")
colnames(mRNA_TFs) = c("TF/miRNA","mRNA")
edges = rbind(mRNA_miRNA_res[mRNA_miRNA_res$`TF/miRNA` %in% miRNA,],mRNA_TFs)
write.table(edges,'results/06.network/nodes.txt',quote = F,row.names = F,sep = '\t')


#05.######
dir.create('results/07.correlation')
##5.1 #####
GSE53625.cibersort=read.delim('results/GSE53625_CIBERSORT_Results.txt',row.names = 1,check.names = F)
head(GSE53625.cibersort)
pdf('results/07.correlation/GSE53625_cibersort.pdf',height = 6,width = 15)
my_mutiboxplot(GSE53625.cibersort[GSE53625.pheno$Samples,1:22], group = GSE53625.pheno$tissue, notch = F,
               group_cols =  c("#00AFBB",  "#FC4E07"),angle =30,ylab = 'score',size = 14)
dev.off()

library(ggbiplot)
before.pca <- prcomp(GSE53625.cibersort[GSE53625.pheno$Samples,1:22], scale=T)
pca_before <- ggbiplot(before.pca, scale=1, groups = GSE53625.pheno$tissue,
                       ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_manual(values =c("#00AFBB",  "#FC4E07")) + 
  theme_light() + xlab('PCA1') + ylab('PCA2') +
  theme(legend.direction = 'horizontal', legend.position = 'top',
        text = element_text(family = 'Times',size = 15),panel.grid = element_blank())
pca_before
ggsave('results/07.correlation/immune_PCA.pdf',pca_before,height = 7,width = 7)


cibersort_cor_res2 <- Hmisc::rcorr(as.matrix(GSE53625.cibersort[GSE53625.pheno$Samples[GSE53625.pheno$tissue=='cancer'],1:22]),type = 'spearman')
cibersort_cor_res2$P[is.na(cibersort_cor_res2$P)] <- 0

cibersort_cor_res <- cor(GSE53625.cibersort[GSE53625.pheno$Samples[GSE53625.pheno$tissue=='cancer'],1:22])

pdf('results/07.correlation/immune_cor.pdf',height = 8,width = 8,onefile = F)
corrplot(cibersort_cor_res, p.mat = cibersort_cor_res2$P, method = 'color', 
         col=colorRampPalette(c('blue', 'white','red'))(100),
         diag = FALSE, type = 'upper',   tl.col = 'black',
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, insig = 'label_sig', 
         pch.col = 'grey20',order = c("original", "AOE", "FPC", "hclust", "alphabet")[1])


dev.off()






##5.2 hub#####
library(ggcorrplot)
library(psych)
library(rcartocolor)

for (i in 1:length(hub.genes)) {
  hubgene_immr_cor <- corr.test(x = as.numeric(GSE53625.exp[hub.genes[i],GSE53625.pheno$Samples[GSE53625.pheno$tissue=='cancer']]),
                                y = GSE53625.cibersort[GSE53625.pheno$Samples[GSE53625.pheno$tissue=='cancer'],1:22],
                                method = "spearman",adjust = "BH",ci = F)
  GSE53625.cibersort_res=data.frame(immune_cell=colnames(GSE53625.cibersort)[1:22])
  GSE53625.cibersort_res$cor<-as.numeric(hubgene_immr_cor$r)
  GSE53625.cibersort_res$p.adj<-as.numeric(hubgene_immr_cor$p.adj)
  head(GSE53625.cibersort_res)
  table(GSE53625.cibersort_res$p.adj<0.05)
  GSE53625.cibersort_res=GSE53625.cibersort_res[order(GSE53625.cibersort_res$cor),]
  head(GSE53625.cibersort_res)
  p[[i]]=ggplot(data=GSE53625.cibersort_res,aes(x=cor,y=reorder(immune_cell,cor), color = -log10(p.adj))) +
    geom_point(aes(size=abs(cor)),show.legend = F) +
    scale_color_continuous(type = "gradient")+
    scale_color_gradient(low = "darkkhaki", high = "darkgreen")+
    geom_segment(aes(yend=immune_cell,xend=0),size=.5) +
    labs(x='spearman Correlation',y='immune cells')+theme_bw()+
    theme(text = element_text(family = 'Times'))+ggtitle(hub.genes[i])
  
}
pdf('results/07.correlation/GSE53625_cibersort_cor.pdf',height = 10,width = 15)
mg_merge_plot(p,ncol=3,nrow = 2,labels = 'B')
dev.off()


save.image(file = 'ESCA_PCD.RData')
