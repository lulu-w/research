library(edgeR)

data<-read.table("./input/counts_matrix.txt",row.names=1,sep="\t")
data_combine<-cbind(as.matrix(rowSums(data[,1:2])),as.matrix(rowSums(data[,3:5])),as.matrix(rowSums(data[,6:7])), as.matrix(rowSums(data[,8:9])), as.matrix(rowSums(data[,10:12])), as.matrix(rowSums(data[,13:15])))

exon_length<-read.table("./input/exon_lengths.tsv",sep="\t")
exon_length<-t(exon_length)
data_combine<-data_combine[rownames(exon_length),]

#Fileter by cpm >= 1
data_nor<-cpm(data_combine)
data_cut<-data_combine[rownames(data_nor[data_nor[,1]>=1 | data_nor[,2]>=1 | data_nor[,3]>=1 | data_nor[,4]>=1 | data_nor[,5]>=1 | data_nor[,6]>=1,]),]

source("makeGFOLDdat.r")
ref<-read.table("./input/refgene_symbol.txt",header=T,sep="\t", stringsAsFactors = F)
rownames(ref)<-make.names(ref[,2],unique=T)


#Choose to normalize by one housekeeping gene ACTB
actb<-data_cut["NM_001101",]


#Norm_factor to adjust the read count of ACTB to be the same
norm_factor<-c(1,data_cut["NM_001101",2]/data_cut["NM_001101",1],data_cut["NM_001101",3]/data_cut["NM_001101",1],data_cut["NM_001101",4]/data_cut["NM_001101",1],data_cut["NM_001101",5]/data_cut["NM_001101",1],data_cut["NM_001101",6]/data_cut["NM_001101",1])


#Choose the transcript that has highest read count to represent gene.
ref1 = ref[, c("name", "name2")]
tmp = data.frame(id = rownames(data_cut), data_cut, sum = rowSums(data_cut))
mdf = merge(tmp, ref1, by.x = "id", by.y = "name")
genename = mdf[order(mdf$name2, mdf$sum, decreasing = T), ]
genename = as.character(genename$id[!duplicated(genename$name2)])

genename = as.data.frame(genename)
rownames(genename)<-genename[,1]

data_cut<-data_cut[rownames(genename),]
ref <- ref[rownames(genename), ]

#Add a fake gene to adjust the RPKM of ACTB to be the same
new_gene<-colSums(data_cut)[6]*data_cut["NM_001101",]/data_cut["NM_001101",6]-colSums(data_cut)
data_cut<-rbind(data_cut,new_gene)
exon_length<-rbind(exon_length,100)
rownames(exon_length)[dim(exon_length)[1]]<-"new_gene"
ref<-rbind(ref,rep(1, ncol(ref)))
rownames(ref)[nrow(ref)]<-"new_gene"
ref[nrow(ref), "name"] <- "new_gene"
ref[nrow(ref), "name2"] <- "new_gene"

condition <- c("sr_treated", "sr_untreated", "sen_treated", "sen_untreated", "rej_treated", "rej_untreated")

data_final<-round(t(t(data_cut)/norm_factor))
for ( i in 1:6)
{
	result<-makeGFOLDdat(count.dat=round(data_final),ref.dat=ref,exon.length=exon_length,i=i)
  print(which(rownames(result) == "NM_001017973"), )
	write.table(result,file=paste("./output/",condition[i],".read_cnt",sep=""),quote=F,sep="\t",col.names=F)
}

#Finish generating inputs for GFOLDS

#Read in outputs from GFOLDS
data1_de<-read.table("sen_untreatedVSsen_treated.diff",header=T,sep="\t")
data1_cnt<-read.table("sen_untreatedVSsen_treated.diff.ext",header=T,sep="\t")
rownames(data1_de)<-data1_de[,1]
rownames(data1_cnt)<-data1_cnt[,1]
data2_de<-read.table("sr_untreatedVSsr_treated.diff",header=T,sep="\t")
data2_cnt<-read.table("sr_untreatedVSsr_treated.diff.ext",header=T,sep="\t")
rownames(data2_de)<-data2_de[,1]
rownames(data2_cnt)<-data2_cnt[,1]
neg<-function(x) -x 

data1_de<-data1_de[rownames(genename),]
data2_de<-data2_de[rownames(genename),]
data1_cnt<-data1_cnt[rownames(genename),]
data2_cnt<-data2_cnt[rownames(genename),]

data1_de_mrna<-data1_de[grepl("NM_",rownames(data1_de),perl=T),]
data2_de_mrna<-data2_de[grepl("NM_",rownames(data2_de),perl=T),]
data1_cnt_mrna<-data1_cnt[grepl("NM_",rownames(data1_de),perl=T),]
data2_cnt_mrna<-data2_cnt[grepl("NM_",rownames(data2_cnt),perl=T),]

data1_up<-rownames(data1_de_mrna[data1_de_mrna[,3]>0.01,])
data1_down<-rownames(data1_de_mrna[data1_de_mrna[,3]<neg(0.01),])
data2_up<-rownames(data2_de_mrna[data2_de_mrna[,3]>0.01,])
data2_down<-rownames(data2_de_mrna[data2_de_mrna[,3]<neg(0.01),])
#data1_de<-data1_de[rownames(genename),]


data1<-cbind(data1_de_mrna,data1_cnt_mrna)
data2<-cbind(data2_de_mrna,data2_cnt_mrna)
write.csv(data1[,c(1,2,3,9,10)],quote=F,file="temp.csv")
