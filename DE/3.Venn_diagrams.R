library(ggplot2)
library(ggVennDiagram)


comps <- c("Samra","Metzelar")
GSEs <- c("GSE1159","GSE6891","GSE13204","GSE15434","GSE22845")

venn <- function(comp,GSE){
	compdf <- read.csv(paste0(comp,".csv"), header=TRUE,stringsAsFactors=FALSE)
	GSEdf <- read.table(Sys.glob(paste0("*",GSE,"*Low_DE.tsv")),
						sep="\t",header=TRUE, stringsAsFactors=FALSE)
	GSEdf <- GSEdf[(abs(GSEdf$logFC) > 1.5) & (GSEdf$P.Value < 0.05),]
	# venn diagram
	venn_list <- list(a=GSEdf$Gene.symbol,
					  b=compdf[[1]])
	names(venn_list) <- c(GSE, comp)
	g <- ggVennDiagram(venn_list)
	ggsave(paste0(GSE,"vs",comp,"_venn.pdf"),g)
}

for (i in comps){
	sapply(GSEs, function(x){venn(i,x)})
}
