library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

inputFile="all.txt"         
gmtFile="c5.go.v7.4.symbols.gmt"     
setwd("C:\\Users\\gifghjd\\Desktop\\136Diagnostic\\11.GSEA")     

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
logFC=as.vector(rt[,2])
names(logFC)=as.vector(rt[,1])
logFC=sort(logFC, decreasing=T)

gmt=read.gmt(gmtFile)

kk=GSEA(logFC, TERM2GENE=gmt, pvalueCutoff = 1)
kkTab=as.data.frame(kk)
kkTab=kkTab[kkTab$p.adjust<0.05,]
write.table(kkTab,file="GSEA.result.txt",sep="\t",quote=F,row.names = F)

termNum=5      
kkUp=kkTab[kkTab$NES>0,]
if(nrow(kkUp)>=termNum){
	showTerm=row.names(kkUp)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Treat")
	pdf(file="GSEA.treat.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}

termNum=5    
kkDown=kkTab[kkTab$NES<0,]
if(nrow(kkDown)>=termNum){
	showTerm=row.names(kkDown)[1:termNum]
	gseaplot=gseaplot2(kk, showTerm, base_size=8, title="Enriched in Control")
	pdf(file="GSEA.con.pdf", width=7, height=5.5)
	print(gseaplot)
	dev.off()
}
