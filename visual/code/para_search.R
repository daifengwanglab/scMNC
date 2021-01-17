library(reshape2)
library(ManiNetCluster)

Dim_red = function(ef,gf,method,d,k_NN,k_medoid){
  #e-feature
  X = apply(ef,2,scale)
  #g-feature
  Y=t(log10(gf+1))
  #Dim_red
  XY_corr=Correspondence(matrix=diag(nrow(X)))
  df=ManiNetCluster(X,Y,nameX='Ephys',nameY='Expr',corr=XY_corr,d=d,
                    method=method,k_NN=k_NN,k_medoids=k_medoid)
  return(df[,-1])
}

ef = read.csv("data/efeature_filtered.csv",header = T,stringsAsFactors = F)
gf = read.csv("data/expMat_filtered.csv",header = T,stringsAsFactors = F)
rownames(ef) = ef[,1];ef = ef[,-1]
rownames(gf) = gf[,1];gf = gf[,-1]
n=nrow(ef)
method = c('linear manifold','cca','manifold warping','nonlinear manifold aln','nonlinear manifold warp')

dist =c()
for (i in 1:5){
  res = Dim_red(ef,gf,method = method[i],d=3L,k_NN=10L,k_medoid=5L)
  res_e.sd = apply(res[res$data=="Ephys",2:4],2,scale)
  res_t.sd = apply(res[res$data=="Expr",2:4],2,scale)
  for (k in 1:n){
    dist = c(dist,sqrt(sum((res_e.sd[k,]-res_t.sd[k,])^2)))
  }
}

### PCA ###
#clustering only use gfs
PCA_res_t = prcomp(t(log10(gf+1)),scale=T,rank=3,retx=T)$x
PCA_res_t = data.frame(Val0 = as.numeric(PCA_res_t[,1]),Val1 = as.numeric(PCA_res_t[,2]),Val2 = as.numeric(PCA_res_t[,3]))
PCA_res_t.sd = apply(PCA_res_t[,1:3],2,scale)
#clustering only use efs
PCA_res_e = prcomp(ef,scale=T,rank=3,retx=T)$x
PCA_res_e = data.frame(Val0 = as.numeric(PCA_res_e[,1]),Val1 = as.numeric(PCA_res_e[,2]),Val2 = as.numeric(PCA_res_e[,3]))
PCA_res_e.sd = apply(PCA_res_e[,1:3],2,scale)
for (k in 1:n){
  dist = c(dist,sqrt(sum((PCA_res_e.sd[k,1:3]-PCA_res_t.sd[k,1:3])^2)))
}

#boxplot
distmat = data.frame(matrix(as.numeric(dist),nrow = nrow(ef)))
colnames(distmat) = c('LM','CCA','MW','NMA','NMW',"PCA")
distmat = distmat[,c(-5)]
distmat = melt(distmat,variable.name = "Method", value.name = "Pairwise distance between g & e features")
boxplot(`Pairwise distance between g & e features`~Method,data=distmat,
        ylab="Pairwise distance between G & E features", xlab="",outline =F,col ="white")
tapply(distmat$`Pairwise distance between g & e features`,INDEX = distmat$Method,FUN = mean)
