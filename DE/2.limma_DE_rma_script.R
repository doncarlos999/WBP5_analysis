#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)
library(ggplot2)
library(affy)
library(pheatmap)
library(RColorBrewer)
mainDir  <- getwd() 
balancedHighLowFile <- "../survival_analysis/WBP5_217975_at_GSE6891_highvslow.csv"
balancedHighLowFile2 <- "../survival_analysis/WBP5_217975_at_GSE15434_highvslow.csv"
balancedHighLowFile3 <- "../survival_analysis/WBP5_217975_at_GSE1159_highvslow.csv"
gene <- stringr::str_split(balancedHighLowFile,"_", simplify=T)[1]
probe <- stringr::str_remove(balancedHighLowFile, paste0(gene,"_"))
probe <- stringr::str_extract(probe, ".*_at")
#GSE <- stringr::str_extract(balancedHighLowFile, "GSE\\d{4,5}")
date_time <- gsub("\\:",".",gsub(" ", "_",Sys.time()))
# subtype options "all_subs", "NN", "noNN"
subtype  <- "all_subs"
# first run fresh should be set as TRUE
# once raw data has been RMA processed once this can be set to FALSE
fresh <- TRUE
# main loop
# Kohlman = "GSE15434"
# Verhaak = "GSE6891"
# Valk = "GSE1159"
# Taskensen = "GSE22845"
for (i in c("GSE15434","GSE13204","GSE6891","GSE1159", "GSE22845")){
#for (i in c("GSE1159")){
    GSE <- i
    if (GSE == "GSE6891"){
            highLow  <- read.csv(balancedHighLowFile,
                 stringsAsFactors=F)
    }
    if (GSE == "GSE1159"){
            highLow  <- read.csv(balancedHighLowFile3,
                 stringsAsFactors=F)
    }
    if (GSE == "GSE15434"){
            highLow  <- read.csv(balancedHighLowFile2,
                 stringsAsFactors=F)
            highLow <- highLow[!duplicated(highLow$patient),]
    }

    gset <- getGEO(filename = Sys.glob(paste0("../survival_analysis/", GSE,"*.gz")),
               GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = '../')
    if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
    gset2 <- gset
    if (GSE %in% c("GSE6891","GSE15434")){
        sel  <- !grepl("score: FAB M3", gset2@phenoData@data$characteristics_ch1.6)
        gset2 <- gset2[ ,sel]
    }
    if (GSE == "GSE13204"){ 
        sel <- grepl("AML|healthy",gset2@phenoData@data$characteristics_ch1.1)
        gset2 <- gset2[ ,sel]
        # remove AML with t(15;17)
        sel <- !grepl("AML with t\\(15", gset2@phenoData@data$characteristics_ch1.1)
        gset2 <- gset2[ ,sel]
    }
    # chnge expression to rma version
	if (fresh){
		directory <- paste0(GSE,"_RAW/")
		setwd(directory)
		system(paste0("tar -xvf ",GSE,"_RAW.tar --directory stash"))
		if (GSE == "GSE6891"){
			system("for i in stash/*.cel.gz;do mv $i ${i%_*}.CEL.gz;done")
		}
		file_names <- paste0("stash/",colnames(gset2), ".CEL.gz")
		file.copy(file_names, ".")
		data <- ReadAffy()
		data.rma <- rma(data)
		rm(data)
		system("rm *.CEL.gz")
		system("rm stash/*.CEL.gz")
		setwd(mainDir)
		saveRDS(data.rma, file=paste0(GSE,"_data.rma.RDS"))
	}
	if (!fresh){
		data.rma <- readRDS(file=paste0(GSE,"_data.rma.RDS"))
	}
	if (GSE == "GSE1159"){
		colnames(gset2)<-gset2@phenoData@data[,2]
	}
	colnames(data.rma) <- stringr::str_remove(colnames(data.rma), ".CEL.gz")
	gset2 <- gset2[,colnames(gset2) %in% colnames(data.rma)]
	data.rma <- data.rma[row.names(data.rma) %in% row.names(gset2),]
	exprs(gset2) <- exprs(data.rma)
	rm(data.rma)
    ## output list of characteristics for potential box plots
    chars <- grep("characteristics", colnames(gset2@phenoData@data))
    out <- capture.output(for (i in chars){print(unique(gset2@phenoData@data[i]))})
    cat(paste0(GSE," unique characteristics"), out,
        file=paste0(GSE,date_time,"_unique_chars.txt"),sep="\n")
    ex_all <- exprs(gset2)
    # make subtype expression boxplots
    if (GSE == "GSE6891"){
        if (!subtype %in% c("NN","noNN")){
            probe_df <- data.frame(pheno=gset2@phenoData@data$characteristics_ch1.7,
                   probe=ex_all[row.names(ex_all) == probe,])
            probe_df <- probe_df[grepl("risk",probe_df$pheno),]
            medians <- vector("double", length(unique(probe_df$pheno)))
            for (i in seq_along(unique(probe_df$pheno))){
                medians[[i]] <- median(probe_df$probe[probe_df$pheno %in% unique(probe_df$pheno)[i]])
                names(medians)[i] <- as.character(unique(probe_df$pheno)[i])
            }
            medians <- sort(medians[!grepl("healthy", names(medians))], decreasing = F)
            probe_df$pheno <- factor(probe_df$pheno, levels =c(names(medians)))
            ggplot(probe_df, aes(pheno, probe)) + geom_boxplot() +
              labs(y=gene) +  coord_flip()
            ggsave(paste(gene,probe,date_time,subtype,GSE,'risk_expression.pdf', sep = "_"))  

            probe_df <- data.frame(pheno=gset2@phenoData@data$characteristics_ch1.9,
                   probe=ex_all[row.names(ex_all) == probe,])
            probe_df <- probe_df[grepl("karyotype",probe_df$pheno),]
            probe_df <- probe_df[!probe_df$pheno %in% c("karyotype: failure",
                                   "karyotype: ",
                                   "karyotype: NN"),]
            medians <- vector("double", length(unique(probe_df$pheno)))
            for (i in seq_along(unique(probe_df$pheno))){
                medians[[i]] <- median(probe_df$probe[probe_df$pheno %in% unique(probe_df$pheno)[i]])
                names(medians)[i] <- as.character(unique(probe_df$pheno)[i])
            }
            medians <- sort(medians[!grepl("healthy", names(medians))], decreasing = F)
            probe_df$pheno <- factor(probe_df$pheno, levels =c(names(medians)))
            g <- ggplot(probe_df, aes(pheno, probe)) + geom_boxplot() +
              labs(y=gene) +  coord_flip()
            ggsave(paste(gene,probe,date_time,subtype,GSE,'karyotype_expression.pdf', sep = "_"),g)  
               }
    }
    if (GSE == "GSE13204"){
        probe_df <- data.frame(pheno=gset2@phenoData@data$characteristics_ch1.1,
               probe=ex_all[row.names(ex_all) == probe,])
        rm(ex_all)
        medians <- vector("double", length(unique(probe_df$pheno)))
        for (i in seq_along(unique(probe_df$pheno))){
          medians[[i]] <- median(probe_df$probe[probe_df$pheno %in% unique(probe_df$pheno)[i]])
          names(medians)[i] <- as.character(unique(probe_df$pheno)[i])
        }
        medians <- sort(medians[!grepl("healthy", names(medians))], decreasing = F)
        probe_df$pheno <- factor(probe_df$pheno, levels =c(names(medians),
                          "leukemia class: Non-leukemia and healthy bone marrow"))
        ggplot(probe_df, aes(pheno, probe)) + geom_boxplot() +
          labs(y=gene) +  coord_flip()
        ggsave(paste(gene,probe,date_time,subtype,GSE,'subtype_expression.pdf', sep = "_"))  
    }
    rm(ex_all)
    ex <- exprs(gset2)
    # make high vs low boxplots for balanced datasets
    if (GSE %in% c("GSE6891", "GSE15434", "GSE1159")){
        gset2  <- gset2[,colnames(gset2) %in% highLow$patient]
        sml   <- vector("character",ncol(gset2))
        for (i in seq_along(colnames(gset2))){
            for (j in seq_along(highLow$patient)){
                if (colnames(gset2)[i] == highLow$patient[j]){
                    sml[[i]] <- highLow$highvlow[j]
                }
            }
        }
        ex  <- ex[,colnames(ex) %in% highLow$patient]
        probeDf <- data.frame(probe=ex[row.names(ex)==probe,],
                      patient=colnames(ex))
        probeDf <- dplyr::left_join(probeDf,highLow)
        probeDf$x <- "x"
        probeDf$highvlow <- factor(probeDf$highvlow, levels=c("high", "low"))
        colorPal <- c("#E41A1C","#377EB8")
        g<-ggplot(probeDf, aes(x, probe, fill=highvlow)) + geom_boxplot() +
          labs(y=gene) + theme_minimal() +
          scale_fill_manual(values = colorPal)
        ggsave(paste(gene,probe,date_time,subtype,GSE,'expression_high_low.pdf',
                 sep="_"),g)
    }
    # upper and lower quantile cut offs for unbalanced datasets
    if (GSE %in% c("GSE13204", "GSE22845")){
        quants <- quantile(ex[row.names(ex) %in% probe,])
        gset2 <- gset2[,(ex[row.names(ex) %in% probe,] < quants[2]) | (ex[row.names(ex) %in% probe,] > quants[4])]
        ex <- ex[,(ex[row.names(ex) %in% probe,] < quants[2]) | (ex[row.names(ex) %in% probe,] > quants[4])]
        # set group names
        sml <- ifelse(ex[row.names(ex) %in% probe,] < quants[2], "low", "high")
        ex  <- ex[,colnames(ex) %in% names(sml)]
        probeDf <- data.frame(probe=ex[row.names(ex)==probe,],
                  patient=colnames(ex))
        highLow <- data.frame(sml)
        colnames(highLow) <- "highvlow"
        highLow$patient <- row.names(highLow)
        probeDf <- dplyr::left_join(probeDf,highLow)
        probeDf$x <- "x"
        colorPal <- c("#E41A1C","#377EB8")
        probeDf$highvlow <- factor(probeDf$highvlow, levels=c("high", "low"))
        g<-ggplot(probeDf, aes(x, probe, fill=highvlow)) + geom_boxplot() +
            labs(y=gene) + theme_minimal() +
            scale_fill_manual(values = colorPal)
        ggsave(paste(gene,probe,date_time,subtype,GSE,'expression_high_low.pdf'),g)
    }
    # box plots of HOXA and HOXB expression
    prob2gene <- read.table("prob2gene.tsv",
                header=T)
    hoxab <- prob2gene[grepl("HOXA|HOXB|MEIS[123]|PBX[123]|FLT3|CEBPA|MYB|FOXC1|GATA2|CRNDE|CLU|CTSG|PRDM16|IGFBP2|CPA3|WT1", prob2gene$Gene.symbol) ,]
    hoxab <- hoxab[!grepl("-AS|MYB[PLB]|MEIS3P|FLT3LG",hoxab$Gene.symbol),]
    hoxab$Gene.symbol <- stringr::str_remove(hoxab$Gene.symbol,".*\\/\\/\\/") 
    hoxab <- dplyr::select(hoxab, ID, Gene.symbol)
    hoxEx <- as.data.frame(ex[row.names(ex) %in% hoxab$ID,])
    hoxEx$ID <- row.names(hoxEx) 
    hoxEx <- dplyr::left_join(hoxEx, hoxab)
    hoxEx <- dplyr::select(hoxEx, -ID)
    hoxEx <- aggregate(. ~ Gene.symbol, hoxEx, median)
    row.names(hoxEx) <- hoxEx$Gene.symbol
    hoxEx <- dplyr::select(hoxEx, -Gene.symbol)
    hoxEx <- as.data.frame(t(hoxEx))
    hoxEx[gene] <- sml
    hoxExMelt <- data.table::melt(hoxEx)
    hoxExMelt$group <- ifelse(grepl("HOX",hoxExMelt$variable),"1","o") 
    hoxExMelt$group <- ifelse(grepl("PBX|MEIS",hoxExMelt$variable),"2",hoxExMelt$group)
    hoxExMelt$group <- ifelse(grepl("GATA|MYB|FOXC1|FLT3|CEBPA",hoxExMelt$variable),"3",hoxExMelt$group)

    colorPal <- c("#E41A1C","#377EB8")
    hoxExMelt$GSE <- GSE
    readr::write_csv(hoxExMelt,paste(gene,probe,date_time,subtype,GSE,
			      'selected_gene_expression_high_low.csv',
			      sep ="_"))
    g <- ggplot(hoxExMelt, aes_string("variable", "value", fill=gene)) +
        geom_boxplot() + coord_flip() + theme_minimal()+
        labs(y="expression", x="gene") + 
        scale_fill_manual(values = colorPal) +
        facet_wrap(vars(group),scales="free_y")
    ggsave(paste(gene,probe,date_time,subtype,GSE,'selected_gene_expression_high_low.pdf',
             sep="_"),g,
           height=7,width=14)
    hox2 <- hoxExMelt[grepl("HOXA[9753]|HOXA10|HOXB[632]",hoxExMelt$variable),]
    g <- ggplot(hox2, aes_string("variable", "value", fill=gene)) +
        geom_boxplot() + coord_flip() + theme_minimal()+
        labs(y="expression", x="gene") + 
        scale_fill_manual(values = colorPal) 
    ggsave(paste(gene,probe,date_time,subtype,GSE,'hox_selected_gene_expression_high_low.pdf',
             sep="_"),g,
           height=7,width=7)
    # set up the data and proceed with DE analysis
    fl <- as.factor(sml)
    gset2$description <- fl
    design <- model.matrix(~ description + 0, gset2)
    colnames(design) <- levels(fl)
    fit <- lmFit(gset2, design)
    for (i in seq(2)){
        if (i == 1){
            contrast <- "LowMinusHigh"
            cont.matrix <- makeContrasts(low-high, levels=design)
        }
        if (i ==2){
            contrast <- "HighMinusLow"
            cont.matrix <- makeContrasts(high-low, levels=design)
        }
        fit2 <- contrasts.fit(fit, cont.matrix)
        fit2 <- eBayes(fit2, 0.01)
        tT <- topTable(fit2, adjust="fdr", sort.by="B", number=10000, p.value = 0.05)
        if (dim(tT)[2] > 0){
        tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
        write.table(tT, file=paste(gene,probe,date_time,subtype,GSE,contrast,'DE.tsv',sep = "_"), row.names=F, sep="\t")
        print("Done writing significant DE genes to file")

        full_table <- topTable(fit2, number=60000)
        full_table2 <- subset(full_table, select=c("ID","logFC","Gene.ID", "Gene.symbol","P.Value", "adj.P.Val"))
        full_table2$fcsign <- sign(full_table2$logFC)
        full_table2$logP <- -log10(full_table2$P.Value)
        full_table2$gsea.metric <- full_table2$logP / full_table2$fcsign
        write.table(full_table2, file=paste(gene,probe,date_time,subtype,GSE,contrast,'GSEA.tsv',sep = "_"), row.names=F, sep="\t")
        full_table2 <- dplyr::select(full_table2, Gene.symbol, gsea.metric)
        write.table(full_table2, file=paste(gene,probe,date_time,subtype,GSE,contrast,'GSEA.rnk',sep = "_"), row.names=F, sep="\t")
        print("Done writing all genes to GSEA file")

        full_table$minuslog10adjP <- -log10(full_table$adj.P.Val) 
        for (fc in c(1.5,2)){
            fccut <- fc
            full_table$color <- ifelse((full_table$minuslog10adjP > 1.3)& (full_table$logFC > fccut), "Up",
                           ifelse((full_table$minuslog10adjP > 1.3) & (full_table$logFC < -fccut),
                              "Down", "ns"))
            full_table$color <- factor(full_table$color, levels = c("Up","ns", "Down"))
            colorPal <- c("#E41A1C","#377EB8","#4DAF4A")
            names(colorPal) <- c("Up", "Down", "ns")
            ymax <- max(full_table$minuslog10adjP)+1
            xmax <- max(abs(full_table$logFC)) + 0.5
            ggplot(full_table, aes(logFC, minuslog10adjP, color=color,
                           label=ifelse(full_table$Gene.symbol %in% 
										c("HOXA9","HOXA7","HOXA10","HOXA3",
										  "HOXA5","HOXB2","HOXB3","HOXB6",
										  "MEIS1", "PBX1","PBX3","FOXC1",
										  "MYB","GATA2","FLT3","CEBPA"),
										full_table$Gene.symbol, NA))) + 
              geom_point(size=4) +
              geom_vline(xintercept = c(-fccut,fccut), linetype='dashed' ) + 
              geom_hline(yintercept = 1.3, linetype = 'dashed') +
              ggrepel::geom_label_repel() +
              labs(y="-log10(adjusted P value") +
			  theme_linedraw() +
              ylim(0,ymax) + 
              xlim(-xmax,xmax) + 
              scale_color_manual(values = colorPal)
            ggsave(paste(gene,probe,date_time,subtype,GSE,contrast,"logfc",fccut,'volcano.png',sep = "_"),
                   height = 10, width = 10)
               }
        } else {
        print("No significant genes")
        }
    }
    # correlation plots and heatmaps
    fccut <- c(1,0.5)
    for (i in seq_along(fccut)){
        fc <- fccut[i]
        deGenes <- full_table[(abs(full_table$logFC) > fc) & 
                      (full_table$adj.P.Val < 0.05),]
        hmapMat <- exprs(gset2)
        hmapMat <- hmapMat[row.names(hmapMat) %in% deGenes$ID,]
        corrmap <- Hmisc::rcorr(hmapMat)
	h <- 10
	w <- 10
        if (GSE == "GSE6891"){
            annoCol <- highLow
            row.names(annoCol) <- annoCol$patient
            info <-as.data.frame(gset2@phenoData@data)
            info <- info[c(2,55:68)]
            info <- dplyr::rename(info, patient=geo_accession)
            annoCol <- dplyr::left_join(annoCol, info)
            row.names(annoCol) <- annoCol$patient
            colnames(annoCol) <- stringr::str_remove(colnames(annoCol),":ch1")
            annoCol <- dplyr::select(annoCol, -`disease state`)
            annoCol <- dplyr::select(annoCol, -`cell type`)
            annoCol <- dplyr::select(annoCol, -gender)
            annoCol <- dplyr::select(annoCol, -cebpa)
            annoCol$`cebpa mutation` <- stringr::str_to_lower(annoCol$`cebpa mutation`)
            annoCol$karyotype[annoCol$karyotype==""] <- "Other"
            annoCol <- dplyr::select(annoCol,-patient)
        }
        if (GSE == "GSE15434"){
            annoCol <- highLow
            row.names(annoCol) <- annoCol$patient
            info <-as.data.frame(gset2@phenoData@data)
            info <- info[c(2,44,47,49)]
            info <- dplyr::rename(info, patient=geo_accession)
            annoCol <- dplyr::left_join(annoCol, info)
            row.names(annoCol) <- annoCol$patient
            colnames(annoCol) <- stringr::str_remove(colnames(annoCol),":ch1")
            annoCol <- dplyr::select(annoCol,-patient)
        }
        if (GSE == "GSE13204"){
            annoCol <- data.frame(sml)
            colnames(annoCol)<-"highvlow"
            annoCol$patient <- row.names(annoCol)
            info <-as.data.frame(gset2@phenoData@data)
            info <- info[c(2,34)]
            info <- dplyr::rename(info, patient=geo_accession)
            annoCol <- dplyr::left_join(annoCol, info)
            row.names(annoCol) <- annoCol$patient
            colnames(annoCol) <- stringr::str_remove(colnames(annoCol),":ch1")
            annoCol <- dplyr::select(annoCol,-patient)
            annoCol[[2]] <- stringr::str_remove(annoCol[[2]],"AML with")  
            annoCol[[2]] <- stringr::str_replace(annoCol[[2]],"AML.*","complex")
            annoCol[[2]] <- stringr::str_replace(annoCol[[2]],"normal.*","NN")
            annoCol[[2]] <- stringr::str_replace(annoCol[[2]],"Non.*","healthy")
        }
        if (GSE == "GSE1159"){
            annoCol <- highLow
            row.names(annoCol) <- annoCol$patient
	    info <- readr::read_csv("GSE1159_moreinfo.csv")
	    info <- dplyr::rename(info, "patient"="Array Samples")
	    info <- info[c(1,7:14)]
            annoCol <- dplyr::left_join(annoCol, info)
            row.names(annoCol) <- annoCol$patient
            annoCol <- dplyr::select(annoCol,-patient)

	}
        if (GSE == "GSE22845"){
            annoCol <- data.frame(sml)
            colnames(annoCol)<-"highvlow"
            annoCol$patient <- row.names(annoCol)
            info <-as.data.frame(gset2@phenoData@data)
            info <- info[c(1,2,30)]
            info <- dplyr::rename(info, patient=geo_accession)
            annoCol <- dplyr::left_join(annoCol, info)
            row.names(annoCol) <- annoCol$patient
            annoCol <- dplyr::select(annoCol,-patient)
	    annoCol$title <- stringr::str_remove(annoCol$title, " \\#.*")
	    colnames(annoCol) <- stringr::str_remove(colnames(annoCol),
						     "genotype/variation \\(")
	    colnames(annoCol) <- stringr::str_remove(colnames(annoCol),
						     "\\):ch1")
	    h <- 10
	    w <- 20
	}
        annoColColors <- list(highvlow=c(high="#E41A1C",low="#377EB8"))
        clust  <- hclust(dist(corrmap$r))
        pheatmap(corrmap$r,
             annotation_col = annoCol,
             color =colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(100),
             breaks=seq(0.5,1,length.out=100),
             filename=paste(gene,probe,date_time,subtype,GSE,"logfc",fc,"corplot.png",sep = "_"),
             cluster_rows=clust,
             cluster_cols=clust,
             dpi=300,
             annotation_colors = annoColColors,
             show_rownames =F, show_colnames=F,
             border_color=NA, width=10, height=10)
        hox <- grep("HOXA",full_table$Gene.symbol)
        hoxP <- full_table$ID[hox]
        names(hoxP) <- full_table$Gene.symbol[hox]
        hmapRownames <- ifelse(row.names(hmapMat) %in% hoxP, names(hoxP),"")
        pheatmap(hmapMat, annotation_col = annoCol,
             color =colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(100),
             filename=paste(gene,probe,date_time,subtype,GSE,"logfc",fc,"rawexpr.png",sep = "_"),
             cluster_cols=clust,
             show_colnames=F,
             annotation_colors = annoColColors,
             labels_row=hmapRownames,
             height=h, width=w)
        hmapMatZ <- scale(hmapMat)
        pheatmap(hmapMatZ, annotation_col = annoCol,
             color =colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100),
             filename=paste(gene,probe,date_time,subtype,GSE,"logfc",fc,"Zexpr.png",sep = "_"),
             annotation_colors = annoColColors,
             cluster_cols=clust,
             show_colnames=F,
             labels_row=hmapRownames,
	     height=h, width=w)
    }
}
# make combined dotplot
myfiles <- Sys.glob(paste0(gene,"_",probe,"*selected_gene*.csv"))
df_list <- vector("list", length(myfiles))
for (i in seq_along(myfiles)){
    df_list[[i]] <- readr::read_csv(myfiles[i])
}
dataFrame <- purrr::reduce(df_list,dplyr::bind_rows)
sel <-  c("HOXA9","HOXA7","HOXA10",
						  "HOXA3","HOXA5","HOXB2","HOXB3","HOXB6","MEIS1",
						  "PBX1","PBX3","FOXC1","MYB","GATA2","FLT3","CEBPA")
dataFrame2 <- dataFrame[dataFrame$variable %in% sel,]
dataFrame2$variable <- factor(dataFrame2$variable, levels=sel)
sel <- c("CRNDE","CLU","CTSG","IGFBP2","CPA3","PRDM16")
dataFrame3 <- dataFrame[dataFrame$variable %in% sel ,]
dataFrame3$variable <- factor(dataFrame3$variable, levels=sel)
sel <- c("HOXA1","HOXA2","HOXA3",
		  "HOXA4","HOXA5","HOXA6","HOXA7","HOXA9","HOXA10",
		  "HOXA11","HOXA13","HOXB1","HOXB2","HOXB3","HOXB4",
		  "HOXB5","HOXB6","HOXB7","HOXB8","HOXB9","HOXB12",
		  "MEIS1","MEIS2","MEIS3","PBX1","PBX2","PBX3")
dataFrame4 <- dataFrame[dataFrame$variable %in% sel,]
dataFrame4$variable <- factor(dataFrame4$variable, levels=sel)
colorPal <- c("#E41A1C","#377EB8")
df_list <- list(dataFrame,dataFrame2, dataFrame3, dataFrame4)
for (i in seq_along(df_list)){
    if (i == 1){h <- 20; w <- 20}
    if (i == 2){h <- 10; w <- 10}
    if (i == 3){h <- 10; w <- 10}
    if (i == 4){h <- 20; w <- 20}
    g <- ggplot(df_list[[i]], aes(GSE,value,colour=WBP5)) +
	geom_boxplot(#coef=0,outlier.shape=NA, 
				 size = 0.25) +
	#geom_jitter(size=0.25, alpha=0.75, width=0.25) +
	#ggbeeswarm::geom_quasirandom(size=0.25, alpha = 0.75) +
	facet_wrap(vars(variable)) +
	scale_colour_manual(values = colorPal) +
	theme_minimal() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1))
	#stat_summary(fun.y="median", geom="point", fill="black",
	#		 shape=95, size=15)
    ggsave(paste(gene,probe,date_time,"combined_boxplot",i,".pdf", sep="_"),g,
	   height=h, width=w)
}
info <- capture.output(sessionInfo())
cat(info, file=paste0(date_time,"sessionInfo", ".txt"), sep="\n")

system(paste0('cp limma_DE_rma_script.R ',
			                date_time,'_limma_DE_rma_script.R'))
#system(paste0('zip -m ',gene,probe,date_time,'_plots_data.zip *',date_time,'*'))
