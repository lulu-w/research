##Read count matrix transform into read_cnt for GFOLD

makeGFOLDdat<-function(count.dat, ref.dat, exon.length,i)
{
	rpkm<-count.dat[,i]/exon.length[rownames(count.dat),1]/colSums(count.dat)[i]*1e9
	result<-cbind(as.character(ref.dat[, 13]),count.dat[, i], exon.length[rownames(count.dat),1], rpkm)	
	return(result)
}

#Usage: result<-makeGFOLDdat(count.dat=data_combine_rm,ref.dat=conversion,lib.size=lib.size,exon.length=exon_length,i=4)
#       write.table(result,file="sen_untreated.read_cnt",quote=F,sep="\t",col.names=F)

