library(plyr)
library(reshape2)
library(stringr)
library(dplyr)

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

match = match[match$t_type != "",]
match = match[match$transcriptomics_sample_id %in% ef_naomit$session_idg,]
match = match[match$genotype %in% c("Lamp5","Pvalb","Serpinf1","Sncg","Sst","Vip"),]
match = match[match$dendrite_type == "aspiny",]
#filter e
cellnames = match$transcriptomics_sample_id
ef = ef_naomit[,-which(names(ef_naomit)=="session_idg")];
rownames(ef) = ef_naomit$session_idg;
ef = ef[cellnames,]                               
#filter g
gf_org = gf_org[!duplicated(gf_org$gene),]
rownames(gf_org) = gf_org$gene;gf_org = gf_org[,-1]
gf_type$Cluster = sapply(gf_type$Cluster,substr,start=1,stop=2)
genenames = unique(gf_type$Gene[gf_type$Cluster %in% c("Ex","In")])
genenames = na.omit(rownames(gf_org)[match(genenames,toupper(rownames(gf_org)))])
#write.csv(data.frame(gene = genenames),"data3/neruonal_genes.csv",row.names = F)
gf = gf_org[intersect(genenames,rownames(gf_org)),cellnames]
#write.csv(gf,"data3/in_use/expMat_inuse.csv")
#write.csv(ef,"data3/in_use/efeature_inuse.csv")