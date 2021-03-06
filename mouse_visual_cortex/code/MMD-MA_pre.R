load(file = "../data/visual_data_filtered.rda")
edata = apply(edata,2,scale)
gdata = log10(t(gdata)+1)
write.table(as.matrix(edata) %*% as.matrix(t(edata)),file = "../mmd-ma/efeature_MMD.tsv",col.names = F,row.names = F)
write.table(as.matrix(gdata) %*% as.matrix(t(gdata)),file = "../mmd-ma/gene_MMD.tsv",col.names = F,row.names = F)
