library(reshape2)
library(ManiNetCluster)
library(plyr)
library(RColorBrewer)
Dim_red = function(edata,gdata,method,cellnames,d,k_NN,k_medoid){
  #e-feature
  X = apply(edata,2,scale)
  #g-feature
  Y=t(log2(gdata+1))
  #Dim_red
  XY_corr=Correspondence(matrix=diag(nrow(X)))
  df=ManiNetCluster(X,Y,nameX='Ephys',nameY='Expr',corr=XY_corr,d=d,
                    method=method,k_NN=k_NN,k_medoids=k_medoid)
  df$cellnames = rep(unlist(cellnames),2)
  return(df[,-1])
}

# load data & gene markers
load('../data/102_visual.rda')
gene_type = read.csv("../data/DER-21_Single_cell_markergenes_UMI.csv",stringsAsFactors = F,header=T)

#feature selection elec
edata=edata[,!colnames(edata)%in%c('ADP (mV)','Latency (ms)')];nf = 12
cellnames=intersect(edata$`name sample`,colnames(gdata))
rownames(edata)=edata$`name sample`
edata=edata[,-c(1,ncol(edata),ncol(edata)-1)]
edata=edata[cellnames,]

#feature selection gene
gdata=gdata[!duplicated(gdata$gene),]
rownames(gdata)=gdata$gene
gene_type$Cluster = sapply(gene_type$Cluster,substr,start=1,stop=2)
genenames = unique(gene_type$Gene[gene_type$Cluster %in% c("Ex","In")])
genenames = na.omit(rownames(gdata)[match(genenames,toupper(rownames(gdata)))]) #1298
gdata = gdata[intersect(genenames,rownames(gdata)),cellnames]

n = nrow(edata)
cellnames = rownames(edata)
method = c('linear manifold','cca','manifold warping','nonlinear manifold aln','nonlinear manifold warp')

### MA ###
NMA_res = Dim_red(edata,gdata,method = method[4],cellnames =cellnames,d=3L,k_NN=2L,k_medoid=5L)
NMA_res_e = NMA_res[NMA_res$data=="Ephys",]
NMA_res_t= NMA_res[NMA_res$data=="Expr",]
### CCA ###
CCA_res = Dim_red(edata,gdata,method = method[2],cellnames =cellnames,d=3L,k_NN=2L,k_medoid=5L)
CCA_res_e = CCA_res[CCA_res$data=="Ephys",]
CCA_res_t= CCA_res[CCA_res$data=="Expr",]
### PCA ###
#clustering only use gdatas
PCA_res_t = prcomp(t(log10(gdata+1)),rank=3,retx=T)$x
PCA_res_t = data.frame(Val0 = as.numeric(PCA_res_t[,1]),Val1 = as.numeric(PCA_res_t[,2]),Val2 = as.numeric(PCA_res_t[,3]))
#clustering only use edatas
PCA_res_e = prcomp(edata,rank=3,retx=T)$x
PCA_res_e = data.frame(Val0 = as.numeric(PCA_res_e[,1]),Val1 = as.numeric(PCA_res_e[,2]),Val2 = as.numeric(PCA_res_e[,3]))
PCA_res = rbind(PCA_res_e[,1:3],PCA_res_t[,1:3])
head(NMA_res)

#Plot result
library(plot3D)
#Figure S1
#NMA
NMA_plot = -rbind(apply(NMA_res[1:n,2:4],2,scale),apply(NMA_res[(n+1):(2*n),2:4],2,scale))
NMA_plot[1:n,1] = NMA_plot[1:n,1] + 0.5
points3D(x=NMA_plot[,1], y=NMA_plot[,2], z=NMA_plot[,3],pch = 19,cex=0.5,bty="g",ticktype = "detailed", theta = 40, phi = 10,
         xlab = "",ylab = "",zlab = "",box=T, axes=F,
         colvar = as.numeric(mapvalues(NMA_res$data,names(table(NMA_res$data)),1:2)),col = brewer.pal(6,"Spectral")[c(1,6)],
         colkey = F)
#CCA
CCA_plot = rbind(apply(CCA_res[1:n,2:4],2,scale),apply(CCA_res[(n+1):(2*n),2:4],2,scale))
CCA_plot[1:n,1] = CCA_plot[1:n,1] + 0.5
points3D(x=CCA_plot[,1], y=CCA_plot[,2], z=CCA_plot[,3],pch = 19,cex=0.5,bty="g",ticktype = "detailed", theta = 40, phi = 10,
         xlab = "",ylab = "",zlab = "",box=T, axes=F,
         colvar = as.numeric(mapvalues(CCA_res$data,names(table(CCA_res$data)),1:2)),col = brewer.pal(6,"Spectral")[c(1,6)],
         colkey = F)
#PCA
PCA_plot = rbind(apply(PCA_res_e[,1:3],2,scale),apply(PCA_res_t[,1:3],2,scale))
PCA_plot[1:n,1] = PCA_plot[1:n,1] + 0.5
points3D(x=PCA_plot[,1], y=PCA_plot[,2], z=PCA_plot[,3],pch = 19,cex=0.5,bty="g",ticktype = "detailed", theta = 40, phi = 10,
         xlab = "",ylab = "",zlab = "",box=T, axes=F,
         colvar = as.numeric(mapvalues(NMA_res$data,names(table(NMA_res$data)),1:2)),col = brewer.pal(6,"Spectral")[c(1,6)],
         colkey = F)