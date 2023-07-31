library(dplyr)
library(ggplot2)
library(ggrepel)

logFCfilter=2              
adj.P.Val.Filter=0.05       
inputFile="all.txt"        
setwd("C:\\Users\\gifghjd\\Desktop\\136Diagnostic\\07.vol")      

rt = read.table(inputFile, header=T, sep="\t", check.names=F)
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")

rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Sig))+
    scale_color_manual(values=c("green", "black","red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
    
p1=p+geom_label_repel(data=filter(rt, ((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter))),
                    box.padding=0.1, point.padding=0.1, min.segment.length=0.05,
                    size=1.8, aes(label=id)) + theme_bw()

pdf(file="vol.pdf", width=7, height=6.1)
print(p1)
dev.off()
