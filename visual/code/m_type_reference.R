library(plyr)
library(reshape2)
library(stringr)
library(dplyr)
library(ManiNetCluster)
library(ggplot2)
library(RColorBrewer)

Dim_red = function(ef,gf,method,cellnames,d,k_NN,k_medoid){
  #e-feature
  X = apply(ef,2,scale)
  #g-feature
  Y=t(log10(gf+1))
  #Dim_red
  XY_corr=Correspondence(matrix=diag(nrow(X)))
  df=ManiNetCluster(X,Y,nameX='Ephys',nameY='Expr',corr=XY_corr,d=d,
                    method=method,k_NN=k_NN,k_medoids=k_medoid)
  df$cellnames = rep(unlist(cellnames),2)
  return(df[,-1])
}

ef_all = read.csv("data/efeature.csv",stringsAsFactors = F)
gf_org = read.csv('data/expMat.csv',header=T,check.names = F)
match = read.csv('data/20200711_patchseq_metadata_mouse.csv',header=T,check.names = F,stringsAsFactors = F)
gf_type = read.csv("data/DER-21_Single_cell_markergenes_UMI.csv",stringsAsFactors = F,header=T)

ef_all$session_id = sapply(ef_all$ID,function(x){na.omit(unlist(strsplit(x, "[^0-9]+")))[3]})
ef_all$subject_id = sapply(ef_all$ID,function(x){na.omit(unlist(strsplit(x, "[^0-9]+")))[2]})
ef_all$session_idg = sapply(ef_all$session_id,FUN = function(x){match$transcriptomics_sample_id[match$ephys_session_id==x]})
match$genotype = sapply(match$t_type,function(x){na.omit(unlist(strsplit(x, " ")))[1]})
match$dendrite_type = mapvalues(match$dendrite_type,from = "sparsely spiny",to = "spiny")

ef_naomit <- ef_all[,c(1:11,
                       12:18,21:23,#ramp
                       24:30,33:35,#long
                       36:42,45:47,#short
                       51)]
ef_naomit = na.omit(ef_naomit)

match = match[match$transcriptomics_sample_id %in% ef_naomit$session_idg,]
match = match[match$dendrite_type %in% c("spiny","aspiny"),]
match = rbind(match[match$dendrite_type == "spiny",],
              match[match$genotype %in% c("Lamp5","Pvalb","Serpinf1","Sncg","Sst","Vip") & match$dendrite_type == "aspiny",])

cellnames = match$transcriptomics_sample_id
ef = ef_naomit[,-which(names(ef_naomit)=="session_idg")];
rownames(ef) = ef_naomit$session_idg;
ef = ef[cellnames,]                               

gf_org = gf_org[!duplicated(gf_org$gene),]
rownames(gf_org) = gf_org$gene;gf_org = gf_org[,-1]
gf_type$Cluster = sapply(gf_type$Cluster,substr,start=1,stop=2)
genenames = unique(gf_type$Gene[gf_type$Cluster %in% c("Ex","In")])
genenames = na.omit(rownames(gf_org)[match(genenames,toupper(rownames(gf_org)))])
#write.csv(data.frame(gene = genenames),"data3/neruonal_genes.csv",row.names = F)
gf = gf_org[intersect(genenames,rownames(gf_org)),cellnames]

dendrite_type = sapply(cellnames,FUN = function(x){match$dendrite_type[match$transcriptomics_sample_id==x]})
n = nrow(ef)
method = c('linear manifold','cca','manifold warping','nonlinear manifold aln','nonlinear manifold warp')
#NMA
NMA_res = Dim_red(ef,gf,method = method[4],cellnames =cellnames,d=3L,k_NN=2L,k_medoid=5L)
NMA_res_e = NMA_res[NMA_res$data=="Ephys",]
NMA_res_t= NMA_res[NMA_res$data=="Expr",]
#Fig S3
library(plot3D)
points3D(x=NMA_res_e$Val0, y=NMA_res_e$Val1, z=NMA_res_e$Val2,pch = 19,cex=0.5,bty="g",ticktype = "detailed", theta = 40, phi = 10,
         xlab = "",ylab = "",zlab = "",
         colvar = as.numeric(mapvalues(dendrite_type,names(table(dendrite_type)),1:2)),col =alpha(brewer.pal(6,"Spectral")[c(2,5)],0.8),
         colkey = F)
points3D(x=NMA_res_t$Val0, y=NMA_res_t$Val1, z=NMA_res_t$Val2,pch = 19,cex=0.5,bty="g",ticktype = "detailed", theta = 40, phi = 10,
         xlab = "",ylab = "",zlab = "",
         colvar = as.numeric(mapvalues(dendrite_type,names(table(dendrite_type)),1:2)),col =alpha(brewer.pal(6,"Spectral")[c(2,5)],0.8),
         colkey = F)
