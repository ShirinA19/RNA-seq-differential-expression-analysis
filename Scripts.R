#' Libs
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(Glimma))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(tximport))


#' Helper files
suppressMessages(source("~/Git/UPSCb/src/R/densityPlot.R"))
#suppressMessages(source("~/Git/UPSCb/src/R/plotMA.R"))
suppressMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))
suppressMessages(source("~/Git/UPSCb/src/R/plot.multidensity.R"))

#' Load saved data
#' Read the sample information
#' 
# make sample table according to tximpor order

samples <- read.table("~/Git/UPSCb/projects/spruce-cone-development/doc/Samples.wheat.salmon.correct.csv", , header = TRUE, sep = ";", quote = "\"'", dec = ".")

PathFilesSalmon <-"~/Git/UPSCb/projects/spruce-cone-development/doc/Salmon_wheat_quant_Dseq"

FilesSalmon <- list.files(PathFilesSalmon, recursive = TRUE, pattern = "quant.sf",
                          full.names = TRUE)



SpecifiedPartofPath <- 12
names(FilesSalmon) <- sapply(strsplit(FilesSalmon, "/"), .subset, SpecifiedPartofPath)


names(FilesSalmon) <- paste0(samples$ID)



txi_Wheat_Salmon <- tximport(FilesSalmon,type="salmon",txOut = T)


counts_salmon_wheat <- txi_Wheat_Salmon$counts
counts_salmon_wheat <- round(counts_salmon_wheat)

colnames(counts_salmon_wheat) <- gsub("\\.+","-",colnames(counts_salmon_wheat))


#' Setup graphics
pal=brewer.pal(8,"Dark2")
pal12 <- brewer.pal(12,"Paired")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' ```{r drop limma,echo=FALSE}
#' detach("package:limma")
#' ```
#' 
#' # Process
#' Re-order the samples
stopifnot(all(colnames(counts_salmon_wheat) == samples$ID))

#' ## DE at the gene level
#' A generic function to do the DE
#' default to a 1% FDR and abs(lfc) >=0.5


".de" <- function(counts_salmon_wheat,samples,design,nam,res.nam=NULL,alpha=0.01,lfc=0.5){
  dds <- DESeqDataSetFromMatrix(
    countData = counts_salmon_wheat,
    colData = samples,
    design = design)
  
  # Variance stabilisation
  vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
  vst <- assay(vsd)
  vst <- vst - min(vst)
  
  # Validate the VST 
  meanSdPlot(vst[rowSums(vst)>0,])
  
  # Differential Expression
  suppressMessages(dds <- DESeq(dds))
  
  # Dispersion Estimation
  #
  # The dispersion estimation is adequate
  plotDispEsts(dds)
  
  
  # Obtain the results 
  # resultsNames(dds)
  if(is.null(res.nam)){
    res <- results(dds)
  } else {
    res <- results(dds,name=res.nam)
  }
  
  
  # Plot the Median vs Average
  # There are many genes that appear differentially expressed at a 1% FDR cutoff
  # The log2 fold-change range is relatively broad, with three extreme
  # values
  #plotMA(res,alpha)
  
  # Plot the log10 odds (i.e. -log10 FDR) vs. log2 fold change
  #
  # The volcano plot shows the same results as the MA plot; a
  # large number of genes show significant fold-changes
  volcanoPlot(res,alpha=alpha,lfc = lfc)
  
  # take a look at some results
  # res[res$padj < 1e-150 & !is.na(res$padj),]
  
  
  # Plot the adjusted p-value histogram
  #
  # Which is almost evenly distributed, with an enrichment for
  # lower p-values (more significant)
  hist(res$padj,breaks=seq(0,1,.01))
  
  # Select genes below alpha
  #
  # Note the 0.5 cutoff on the log fold change is motivated by the
  # Schurch et al. RNA, 2016 publication
  sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
  message(sprintf("There are %s genes differentially expressed at a %s cutoff",
                  sum(sel),alpha))
  
  # Equally distributed between F (-1) and M (1)
  pander(table(sign(res[sel,"log2FoldChange"])))
  
  
  
  # Write them out
  outdir=file.path("analysis/DESeq2",nam)
  dir.create(outdir,
             showWarnings=FALSE,recursive=TRUE)
  write.csv(res,
            file=file.path(outdir,paste0(nam,".csv")))
  
  # write the de selected genes
  de.sel <- rownames(res)[sel]
  write(sub("\\.1$","",de.sel),file=file.path(outdir,paste0(nam,"-DE-geneIDs.txt")))
  
  # #### Cluster the VST expression of the genotype genes
  # The color code for the column next to the row dendogram (genes)
  # is red for significant positive fold change and blue for
  # significant negative fold change
  heatmap.2(as.matrix(vst[rownames(vst) %in% de.sel,]),
            scale="row",labRow=NA,trace="none",
            RowSideColors=ifelse(sign(res[sel,"log2FoldChange"])==-1,
                                 "blue",
                                 "red"),
            col=hpal,cexCol = 0.9, margins = c(6.1,0.1),
            labCol=samples[s.sel,"ID"])
  
  # #### Interactive MDS
  # dir.create(file.path(outdir,nam),showWarnings = FALSE)
  # res.df <- as.data.frame(res)
  # res.df$log10MeanNormCount <- log10(res.df$baseMean)
  # idx <- rowSums(counts_salmon_wheat(dds)) > 0
  # res.df <- res.df[idx,]
  # res.df$padj[is.na(res.df$padj)] <- 1
  # rownames(res.df) <- sub("\\.1$","",rownames(res.df))
  # glMDPlot(res.df,
  #          xval = "log10MeanNormCount",
  #          yval = "log2FoldChange",
  #          counts = counts_salmon_wheat(dds)[idx,],
  #          groups = as.character(dds$Type),
  #          samples = colnames(dds),
  #          status = sign(res[idx,"log2FoldChange"]) * sel[idx],
  #          display.columns = c("GeneID", "Confidence","padj","baseMean"),
  #          path = outdir,
  #          folder = nam,
  #          launch=FALSE)
  
  return(res)
}  



### wheat, BrKr , d-0, d-1

s1.sel <- which(samples$Sampling == "d-0" & samples$Type=="BrKr")
s2.sel <- which(samples$Sampling == "d-1" & samples$Type=="BrKr")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Sampling",drop=FALSE],~Sampling,nam="wheat.BrKr.d-0.vs.d-1")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.BrKr.d.0.vs.d.1 <- rownames(res)[sel]





### wheat, PBrC vs. BrKr , ni-0, d-0

s1.sel <- which(samples$Sampling == "ni-0" & samples$Type=="PBrC")
s2.sel <- which(samples$Sampling == "d-0" & samples$Type=="BrKr")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Type",drop=FALSE],~Type,nam="wheat.PBrC.vs.BrKr.ni-0.d-0_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.PBrC.vs.BrKr.ni.0.d.0.new <- rownames(res)[sel]



### wheat, PBrC vs. BrKr , ni-0, d-1

s1.sel <- which(samples$Sampling == "ni-0" & samples$Type=="PBrC")
s2.sel <- which(samples$Sampling == "d-1" & samples$Type=="BrKr")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Type",drop=FALSE],~Type,nam="wheat.PBrC.vs.BrKr.ni-0.d-1_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.PBrC.vs.BrKr.ni.0.d.1.new <- rownames(res)[sel]

### wheat, BrKr , d-0, d-2

s1.sel <- which(samples$Sampling == "d-0" & samples$Type=="BrKr")
s2.sel <- which(samples$Sampling == "d-2" & samples$Type=="BrKr")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Sampling",drop=FALSE],~Sampling,nam="wheat.BrKr.d-0.vs.d-2")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.BrKr.d.0.vs.d.2 <- rownames(res)[sel]



### wheat, PBrC vs. BrKr , ni-0, d-2

s1.sel <- which(samples$Sampling == "ni-0" & samples$Type=="PBrC")
s2.sel <- which(samples$Sampling == "d-2" & samples$Type=="BrKr")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Type",drop=FALSE],~Type,nam="wheat.PBrC.vs.BrKr.ni-0.d-2_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.PBrC.vs.BrKr.ni.0.d.2.new <- rownames(res)[sel]

### wheat, BrKr , d-0, d-3

s1.sel <- which(samples$Sampling == "d-0" & samples$Type=="BrKr")
s2.sel <- which(samples$Sampling == "d-3" & samples$Type=="BrKr")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Sampling",drop=FALSE],~Sampling,nam="wheat.BrKr.d-0.vs.d-3")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.BrKr.d.0.vs.d.3 <- rownames(res)[sel]




### wheat, PBrC vs. BrKr , ni-0, d-3

s1.sel <- which(samples$Sampling == "ni-0" & samples$Type=="PBrC")
s2.sel <- which(samples$Sampling == "d-3" & samples$Type=="BrKr")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Type",drop=FALSE],~Type,nam="wheat.PBrC.vs.BrKr.ni-0.d-3_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.PBrC.vs.BrKr.ni.0.d.3.new <- rownames(res)[sel]


### wheat, BrKr , d-0, d-8

s1.sel <- which(samples$Sampling == "d-0" & samples$Type=="BrKr")
s2.sel <- which(samples$Sampling == "d-8" & samples$Type=="BrKr")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Sampling",drop=FALSE],~Sampling,nam="wheat.BrKr.d-0.vs.d-8")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.BrKr.d.0.vs.d.8 <- rownames(res)[sel]




### wheat, PBrC vs. BrKr , ni-0, d-8

s1.sel <- which(samples$Sampling == "ni-0" & samples$Type=="PBrC")
s2.sel <- which(samples$Sampling == "d-8" & samples$Type=="BrKr")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Type",drop=FALSE],~Type,nam="wheat.PBrC.vs.BrKr.ni-0.d-8_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.PBrC.vs.BrKr.ni.0.d.8.new <- rownames(res)[sel]

##############
### wheat, PBrC vs. BrWar , ni-0, d-0

s1.sel <- which(samples$Sampling == "ni-0" & samples$Type=="PBrC")
s2.sel <- which(samples$Sampling == "d-0" & samples$Type=="BrWar")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Type",drop=FALSE],~Type,nam="wheat.PBrC.vs.BrWar.ni-0.d-0_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.PBrC.vs.BrWar.ni.0.d.0.new <- rownames(res)[sel]


### wheat, BrWar , d-0, d-1

s1.sel <- which(samples$Sampling == "d-0" & samples$Type=="BrWar")
s2.sel <- which(samples$Sampling == "d-1" & samples$Type=="BrWar")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Sampling",drop=FALSE],~Sampling,nam="wheat.BrWar.d-0.vs.d-1")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.BrWar.d.0.vs.d.1 <- rownames(res)[sel]




### wheat, PBrC vs. BrWar , ni-0, d-1

s1.sel <- which(samples$Sampling == "ni-0" & samples$Type=="PBrC")
s2.sel <- which(samples$Sampling == "d-1" & samples$Type=="BrWar")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Type",drop=FALSE],~Type,nam="wheat.PBrC.vs.BrWar.ni-0.d-1_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.PBrC.vs.BrWar.ni.0.d.1.new <- rownames(res)[sel]





### wheat, BrWar , d-0, d-1

s1.sel <- which(samples$Sampling == "d-0" & samples$Type=="BrWar")
s2.sel <- which(samples$Sampling == "d-1" & samples$Type=="BrWar")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Sampling",drop=FALSE],~Sampling,nam="wheat.BrWar.d-0.vs.d-1")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.BrWar.d.0.vs.d.1 <- rownames(res)[sel]


### wheat, BrWar , d-0, d-2

s1.sel <- which(samples$Sampling == "d-0" & samples$Type=="BrWar")
s2.sel <- which(samples$Sampling == "d-2" & samples$Type=="BrWar")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Sampling",drop=FALSE],~Sampling,nam="wheat.BrWar.d-0.vs.d-2")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.BrWar.d.0.vs.d.2 <- rownames(res)[sel]



### wheat, PBrC vs. BrWar , ni-0, d-2

s1.sel <- which(samples$Sampling == "ni-0" & samples$Type=="PBrC")
s2.sel <- which(samples$Sampling == "d-2" & samples$Type=="BrWar")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Type",drop=FALSE],~Type,nam="wheat.PBrC.vs.BrWar.ni-0.d-2_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.PBrC.vs.BrWar.ni.0.d.2.new <- rownames(res)[sel]


### wheat, BrWar , d-0, d-3

s1.sel <- which(samples$Sampling == "d-0" & samples$Type=="BrWar")
s2.sel <- which(samples$Sampling == "d-3" & samples$Type=="BrWar")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Sampling",drop=FALSE],~Sampling,nam="wheat.BrWar.d-0.vs.d-3")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.BrWar.d.0.vs.d.3 <- rownames(res)[sel]



### wheat, PBrC vs. BrWar , ni-0, d-3

s1.sel <- which(samples$Sampling == "ni-0" & samples$Type=="PBrC")
s2.sel <- which(samples$Sampling == "d-3" & samples$Type=="BrWar")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Type",drop=FALSE],~Type,nam="wheat.PBrC.vs.BrWar.ni-0.d-3_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.PBrC.vs.BrWar.ni.0.d.3.new <- rownames(res)[sel]


### wheat, BrWar , d-0, d-8

s1.sel <- which(samples$Sampling == "d-0" & samples$Type=="BrWar")
s2.sel <- which(samples$Sampling == "d-8" & samples$Type=="BrWar")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Sampling",drop=FALSE],~Sampling,nam="wheat.BrWar.d-0.vs.d-8")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.BrWar.d.0.vs.d.8 <- rownames(res)[sel]



### wheat, PBrC vs. BrWar , ni-0, d-8

s1.sel <- which(samples$Sampling == "ni-0" & samples$Type=="PBrC")
s2.sel <- which(samples$Sampling == "d-8" & samples$Type=="BrWar")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts_salmon_wheat[,s.sel],samples[s.sel,"Type",drop=FALSE],~Type,nam="wheat.PBrC.vs.BrWar.ni-0.d-8_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
wheat.PBrC.vs.BrWar.ni.0.d.8.new <- rownames(res)[sel]



#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.PBrC.vs.BrWar.ni.0.d.0.new=wheat.PBrC.vs.BrWar.ni.0.d.0.new,
  wheat.PBrC.vs.BrWar.ni.0.d.1.new=wheat.PBrC.vs.BrWar.ni.0.d.1.new,
  wheat.PBrC.vs.BrWar.ni.0.d.2.new=wheat.PBrC.vs.BrWar.ni.0.d.2.new,
  wheat.PBrC.vs.BrWar.ni.0.d.3.new=wheat.PBrC.vs.BrWar.ni.0.d.3.new,
  wheat.PBrC.vs.BrWar.ni.0.d.8.new=wheat.PBrC.vs.BrWar.ni.0.d.8.new),
  filename=NULL,
  col=pal[1:5],
  category.names=c("wheat.PBrC.vs.BrWar.ni.0.d.0.new","wheat.PBrC.vs.BrWar.ni.0.d.1.new","wheat.PBrC.vs.BrWar.ni.0.d.2.new","wheat.PBrC.vs.BrWar.ni.0.d.3.new","wheat.PBrC.vs.BrWar.ni.0.d.8.new")))


#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.PBrC.vs.BrKr.ni.0.d.0.new=wheat.PBrC.vs.BrKr.ni.0.d.0.new,
  wheat.PBrC.vs.Brkr.ni.0.d.1.new=wheat.PBrC.vs.BrKr.ni.0.d.1.new,
  wheat.PBrC.vs.BrKr.ni.0.d.2.new=wheat.PBrC.vs.BrKr.ni.0.d.2.new,
  wheat.PBrC.vs.BrKr.ni.0.d.3.new=wheat.PBrC.vs.BrKr.ni.0.d.3.new,
  wheat.PBrC.vs.BrKr.ni.0.d.8.new=wheat.PBrC.vs.BrKr.ni.0.d.8.new),
  filename=NULL,
  col=pal[1:5],
  category.names=c("wheat.PBrC.vs.BrKr.ni.0.d.0.new","wheat.PBrC.vs.BrKr.ni.0.d.1.new","wheat.PBrC.vs.BrKr.ni.0.d.2.new","wheat.PBrC.vs.BrKr.ni.0.d.3.new","wheat.PBrC.vs.BrKr.ni.0.d.8.new")))



#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.PBrC.vs.BrKr.ni.0.d.0.new=wheat.PBrC.vs.BrKr.ni.0.d.0.new,
  wheat.PBrC.vs.BrWar.ni.0.d.0.new=wheat.PBrC.vs.BrWar.ni.0.d.0.new),
  filename=NULL,
  col=pal[1:2],
  category.names=c("wheat.PBrC.vs.BrKr.ni.0.d.0.new","wheat.PBrC.vs.BrWar.ni.0.d.0.new")))


#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.PBrC.vs.BrKr.ni.0.d.1.new=wheat.PBrC.vs.BrKr.ni.0.d.1.new,
  wheat.PBrC.vs.BrWar.ni.0.d.1.new=wheat.PBrC.vs.BrWar.ni.0.d.1.new),
  filename=NULL,
  col=pal[1:2],
  category.names=c("wheat.PBrC.vs.BrKr.ni.0.d.1.new","wheat.PBrC.vs.BrWar.ni.0.d.1.new")))



#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.PBrC.vs.BrKr.ni.0.d.2.new=wheat.PBrC.vs.BrKr.ni.0.d.2.new,
  wheat.PBrC.vs.BrWar.ni.0.d.2.new=wheat.PBrC.vs.BrWar.ni.0.d.2.new),
  filename=NULL,
  col=pal[1:2],
  category.names=c("wheat.PBrC.vs.BrKr.ni.0.d.2.new","wheat.PBrC.vs.BrWar.ni.0.d.2.new")))



#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.PBrC.vs.BrKr.ni.0.d.3.new=wheat.PBrC.vs.BrKr.ni.0.d.3.new,
  wheat.PBrC.vs.BrWar.ni.0.d.3.new=wheat.PBrC.vs.BrWar.ni.0.d.3.new),
  filename=NULL,
  col=pal[1:2],
  category.names=c("wheat.PBrC.vs.BrKr.ni.0.d.3.new","wheat.PBrC.vs.BrWar.ni.0.d.3.new")))


#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.PBrC.vs.BrKr.ni.0.d.8.new=wheat.PBrC.vs.BrKr.ni.0.d.8.new,
  wheat.PBrC.vs.BrWar.ni.0.d.8.new=wheat.PBrC.vs.BrWar.ni.0.d.8.new),
  filename=NULL,
  col=pal[1:2],
  category.names=c("wheat.PBrC.vs.BrKr.ni.0.d.8.new","wheat.PBrC.vs.BrWar.ni.0.d.8.new")))


#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.PBrC.vs.BrKr.ni.0.d.0.new=wheat.PBrC.vs.BrKr.ni.0.d.0.new,
  wheat.PBrC.vs.BrKr.ni.0.d.1.new=wheat.PBrC.vs.BrKr.ni.0.d.1.new,
  wheat.PBrC.vs.BrKr.ni.0.d.2.new=wheat.PBrC.vs.BrKr.ni.0.d.2.new,
  wheat.PBrC.vs.BrKr.ni.0.d.3.new=wheat.PBrC.vs.BrKr.ni.0.d.3.new,
  wheat.PBrC.vs.BrKr.ni.0.d.8.new=wheat.PBrC.vs.BrKr.ni.0.d.8.new,
  wheat.PBrC.vs.BrWar.ni.0.d.0.new=wheat.PBrC.vs.BrWar.ni.0.d.0.new,
  wheat.PBrC.vs.BrWar.ni.0.d.1.new=wheat.PBrC.vs.BrWar.ni.0.d.1.new,
  wheat.PBrC.vs.BrWar.ni.0.d.2.new=wheat.PBrC.vs.BrWar.ni.0.d.2.new,
  wheat.PBrC.vs.BrWar.ni.0.d.3.new=wheat.PBrC.vs.BrWar.ni.0.d.3.new,
  wheat.PBrC.vs.BrWar.ni.0.d.8.new=wheat.PBrC.vs.BrWar.ni.0.d.8.new),
  filename=NULL,
  col=pal[1:10],
  category.names=c("wheat.PBrC.vs.BrKr.ni.0.d.0.new","wheat.PBrC.vs.BrKr.ni.0.d.1.new","wheat.PBrC.vs.BrKr.ni.0.d.2.new","wheat.PBrC.vs.BrKr.ni.0.d.3.new","wheat.PBrC.vs.BrKr.ni.0.d.8.new","wheat.PBrC.vs.BrWar.ni.0.d.0.new","wheat.PBrC.vs.BrWar.ni.0.d.1.new","wheat.PBrC.vs.BrWar.ni.0.d.2.new","wheat.PBrC.vs.BrWar.ni.0.d.3.new","wheat.PBrC.vs.BrWar.ni.0.d.8.new")))




#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.PBrC.vs.BrKr.ni.0.d.0.new=wheat.PBrC.vs.BrKr.ni.0.d.0.new,
  wheat.PBrC.vs.BrKr.ni.0.d.1.new=wheat.PBrC.vs.BrKr.ni.0.d.1.new,
  wheat.PBrC.vs.BrKr.ni.0.d.2.new=wheat.PBrC.vs.BrKr.ni.0.d.2.new,
  wheat.PBrC.vs.BrKr.ni.0.d.3.new=wheat.PBrC.vs.BrKr.ni.0.d.3.new,
  wheat.PBrC.vs.BrKr.ni.0.d.8.new=wheat.PBrC.vs.BrKr.ni.0.d.8.new),
  filename=NULL,
  col=pal[1:5],
  category.names=c("wheat.PBrC.vs.BrKr.ni.0.d.0.new","wheat.PBrC.vs.BrKr.ni.0.d.1.new","wheat.PBrC.vs.BrKr.ni.0.d.2.new","wheat.PBrC.vs.BrKr.ni.0.d.3.new","wheat.PBrC.vs.BrKr.ni.0.d.8.new")))


#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.BrKr.d.0.vs.d.1=wheat.BrKr.d.0.vs.d.1,
  wheat.BrKr.d.0.vs.d.2=wheat.BrKr.d.0.vs.d.2,
  wheat.BrKr.d.0.vs.d.3=wheat.BrKr.d.0.vs.d.3,
  wheat.BrKr.d.0.vs.d.8=wheat.BrKr.d.0.vs.d.8),
  filename=NULL,
  col=pal[1:4],
  category.names=c("wheat.BrKr.d.0.vs.d.1","wheat.BrKr.d.0.vs.d.2","wheat.BrKr.d.0.vs.d.3","wheat.BrKr.d.0.vs.d.8")))

#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.BrWar.d.0.vs.d.1=wheat.BrWar.d.0.vs.d.1,
  wheat.BrWar.d.0.vs.d.2=wheat.BrWar.d.0.vs.d.2,
  wheat.BrWar.d.0.vs.d.3=wheat.BrWar.d.0.vs.d.3,
  wheat.BrWar.d.0.vs.d.8=wheat.BrWar.d.0.vs.d.8),
  filename=NULL,
  col=pal[1:4],
  category.names=c("wheat.BrWar.d.0.vs.d.1","wheat.BrWar.d.0.vs.d.2","wheat.BrWar.d.0.vs.d.3","wheat.BrWar.d.0.vs.d.8")))


############################

#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.BrKr.d.0.vs.d.1=wheat.BrKr.d.0.vs.d.1,
  wheat.BrWar.d.0.vs.d.1=wheat.BrWar.d.0.vs.d.1),
  filename=NULL,
  col=pal[1:2],
  category.names=c("wheat.BrKr.d.0.vs.d.1","wheat.BrWar.d.0.vs.d.1")))


#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.BrKr.d.0.vs.d.2=wheat.BrKr.d.0.vs.d.2,
  wheat.BrWar.d.0.vs.d.2=wheat.BrWar.d.0.vs.d.2),
  filename=NULL,
  col=pal[1:2],
  category.names=c("wheat.BrKr.d.0.vs.d.2","wheat.BrWar.d.0.vs.d.2")))


#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.BrKr.d.0.vs.d.3=wheat.BrKr.d.0.vs.d.3,
  wheat.BrWar.d.0.vs.d.3=wheat.BrWar.d.0.vs.d.3),
  filename=NULL,
  col=pal[1:2],
  category.names=c("wheat.BrKr.d.0.vs.d.3","wheat.BrWar.d.0.vs.d.3")))


#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  wheat.BrKr.d.0.vs.d.8=wheat.BrKr.d.0.vs.d.8,
  wheat.BrWar.d.0.vs.d.8=wheat.BrWar.d.0.vs.d.8),
  filename=NULL,
  col=pal[1:2],
  category.names=c("wheat.BrKr.d.0.vs.d.8","wheat.BrWar.d.0.vs.d.8")))






## wheat, BrKr vs. BrWar. d-0

s.sel <- which(samples$Sampling == "d-0" & samples$Type %in% c("BrKr","BrWar") & samples$Plant_path=="Br")
s.sel.df <- samples[s.sel,"Type",drop=FALSE]
res <- .de(counts_salmon_wheat[,s.sel],s.sel.df,~Type,nam="Puc_BrKr_BrWar_d_0_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
Puc.BrKr.BrWar.d.0.new <- rownames(res)[sel]



## wheat, BrKr vs. BrWar. d-1

s.sel <- which(samples$Sampling == "d-1" & samples$Type %in% c("BrKr","BrWar") & samples$Plant_path=="Br")
s.sel.df <- samples[s.sel,"Type",drop=FALSE]
res <- .de(counts_salmon_wheat[,s.sel],s.sel.df,~Type,nam="Puc_BrKr_BrWar_d_1_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
Puc.BrKr.BrWar.d.1.new <- rownames(res)[sel]

## wheat, BrKr vs. BrWar. d-2

s.sel <- which(samples$Sampling == "d-2" & samples$Type %in% c("BrKr","BrWar") & samples$Plant_path=="Br")
s.sel.df <- samples[s.sel,"Type",drop=FALSE]
res <- .de(counts_salmon_wheat[,s.sel],s.sel.df,~Type,nam="Puc_BrKr_BrWar_d_2_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
Puc.BrKr.BrWar.d.2.new <- rownames(res)[sel]

## wheat, BrKr vs. BrWar. d-3

s.sel <- which(samples$Sampling == "d-3" & samples$Type %in% c("BrKr","BrWar") & samples$Plant_path=="Br")
s.sel.df <- samples[s.sel,"Type",drop=FALSE]
res <- .de(counts_salmon_wheat[,s.sel],s.sel.df,~Type,nam="Puc_BrKr_BrWar_d_3_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
Puc.BrKr.BrWar.d.3.new <- rownames(res)[sel]


## wheat, BrKr vs. BrWar. d-8

s.sel <- which(samples$Sampling == "d-8" & samples$Type %in% c("BrKr","BrWar") & samples$Plant_path=="Br")
s.sel.df <- samples[s.sel,"Type",drop=FALSE]
res <- .de(counts_salmon_wheat[,s.sel],s.sel.df,~Type,nam="Puc_BrKr_BrWar_d_8_new")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
Puc.BrKr.BrWar.d.8.new <- rownames(res)[sel]




#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  Puc.BrKr.BrWar.d.0.new=Puc.BrKr.BrWar.d.0.new,
  Puc.BrKr.BrWar.d.1.new=Puc.BrKr.BrWar.d.1.new,
  Puc.BrKr.BrWar.d.2.new=Puc.BrKr.BrWar.d.2.new,
  Puc.BrKr.BrWar.d.3.new=Puc.BrKr.BrWar.d.3.new,
  Puc.BrKr.BrWar.d.8.new=Puc.BrKr.BrWar.d.8.new),
  filename=NULL,
  col=pal[1:5],
  category.names=c("Puc.BrKr.BrWar.d.0.new","Puc.BrKr.BrWar.d.1.new","Puc.BrKr.BrWar.d.2.new","Puc.BrKr.BrWar.d.3.new","Puc.BrKr.BrWar.d.8.new")))



#' ### Venndiagrams
plot.new()
grid.draw(venn.diagram(list(
  Puc.BrKr.BrWar.d.0.new=Puc.BrKr.BrWar.d.0.new,
  Puc.BrKr.BrWar.d.1.new=Puc.BrKr.BrWar.d.1.new),
  filename=NULL,
  col=pal[1:1],
  category.names=c("Puc.BrKr.BrWar.d.0.new","Puc.BrKr.BrWar.d.1.new")))




#' ### Puc BrKr vs. BrWar, d-1
s.sel <- which(samples$Sampling == "d-0" & samples$Type %in% c("BrKr","BrWar") & samples$Plant_path=="Br")
s.sel.df <- samples[s.sel,"Type",drop=FALSE]
s.sel.df$Type %<>% relevel("BrKr")
res <- .de(counts_salmon[,s.sel],s.sel.df,~Type,nam="Puc_BrKr_BrWar_d_1")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
Puc.BrKr.BrWar.d.1 <- rownames(res)[sel]




s.sel <- which(samples$Type == "BrKr" & samples$Sampling %in% c("d-0","d-1") & samples$Plant_path=="Br")
s.sel.df <- samples[s.sel,"Sampling",drop=FALSE]
s.sel.df$Sampling %<>% relevel("d-1")
res <- .de(counts_salmon[,s.sel],s.sel.df,~Sampling,nam="Puc_BrKr_BrWar_d_0_1_test")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
Puc.BrKr.BrWar.d.0.vs.1 <- rownames(res)[sel]


s.sel <- which(samples$Type == "BrWar" & samples$Sampling %in% c("d-0","d-2") & samples$Plant_path=="Br")
s.sel.df <- samples[s.sel,"Sampling",drop=FALSE]
s.sel.df$Sampling %<>% relevel("d-0")
res <- .de(counts_salmon[,s.sel],s.sel.df,~Sampling,nam="Puc_BrWar_d_0_2_test")
alpha=0.01
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
Puc.BrWar.d.0.vs.2 <- rownames(res)[sel]



s1.sel <- which(samples$Sampling == "d-0" & samples$Type=="BrKr")
s2.sel <- which(samples$Sampling == "d-1" & samples$Type=="BrWar")
s.sel <- c(s1.sel,s2.sel)
res <- .de(counts[,s.sel],samples[s.sel,"Type",drop=FALSE],~Type,nam="Puc_BrWar_BrKr_d_0_1_test")
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
Puc.BrWar.BrKr.d.0.vs.1 <- rownames(res)[sel]
