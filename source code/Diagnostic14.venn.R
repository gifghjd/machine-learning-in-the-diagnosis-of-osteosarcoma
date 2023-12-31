library(venn)                  
outFile="interGenes.txt"        
setwd("C:\\Users\\gifghjd\\Desktop\\136Diagnostic\\14.venn")   
geneList=list()

rt=read.table("LASSO.gene.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])      
uniqGene=unique(geneNames)      
geneList[["LASSO"]]=uniqGene

rt=read.table("SVM-RFE.gene.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])      
uniqGene=unique(geneNames)       
geneList[["SVM-RFE"]]=uniqGene

mycol=c("blue2","red2")
pdf(file="venn.pdf", width=5, height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F,ilabels=F)
dev.off()

intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile, intersectGenes, sep="\t", quote=F, col.names=F, row.names=F)
