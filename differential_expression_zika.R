### Differential expression of Zika
### Andrew R Gross, original file: 2016-12-28; Revised file: 2017-04-12
### Calculate differential expression of Zika proteomics data for figures
### Input: Zika proteomics data
### Output: Diff. Ex. data; Diff. Ex. heatmap; line graph comparisons of specific genes

########################################################################
### Header

library(DESeq2)
library(ggplot2)
library(reshape)
library(biomaRt)
#ensembl = useMart(host="www.ensembl.org")
ensembl = useMart(host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
listMarts(host="www.ensembl.org")
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL')
listDatasets(ensembl)
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)
grep('uniprot', attributes[,1])
filters[grep('uniprot', filters[,1]),]

########################################################################
### Functions

addGene <- function(dataframe) {
  genes <- getBM(attributes=c('uniprot_genename','external_gene_name'), filters='uniprot_genename', values=row.names(dataframe), mart=ensembl)
  genes <- genes[match(row.names(dataframe),genes[,1]),]
  Gene <- c()
  for (rowNumber in 1:length(genes[,1])) {
    newGene <- genes[rowNumber,][,2]
    Gene <- c(Gene, newGene)
  }
  dataframe[length(dataframe)+1] <- Gene
  names(dataframe)[ncol(dataframe)] <- "Gene"
  return(dataframe)
}

uniprot_gn.to.gene.name <- function(dataframe) {
  ### Generate initial conversion table
  conversion.df <- getBM(attributes=c('uniprot_gn','external_gene_name'), filters='uniprot_gn', values= row.names(dataframe), mart=ensembl)
  unique.ids <- unique(conversion.df$uniprot_gn)
  ### Identify the redundancy of each ID
  redundancy.df <- data.frame(unique.ids, 'no.of.results' = rep('',length(unique.ids)), stringsAsFactors = FALSE)
  for (row.num in 1:nrow(redundancy.df)) {
    id = redundancy.df[row.num,][,1]
    hits <- grep(id,conversion.df[,1])
    redundancy.df[row.num,][2] <- length(hits)
  }
  ### Generate a df containing single-result IDs and multiple-result IDs
  single.result.df <- data.frame('uniprot_gn' = c(), 'external_gene_name' = c(), stringsAsFactors = FALSE)
  multiple.result.df <- data.frame('uniprot_gn' = c(), 'external_gene_name' = c(), stringsAsFactors = FALSE)
  
  for (row.num in 1:nrow(redundancy.df)) {
    if (redundancy.df$no.of.results[row.num] == 1) {
      row.to.add <- conversion.df[grep(redundancy.df$unique.ids[row.num], conversion.df$uniprot_gn),]
      single.result.df <- rbind(single.result.df,row.to.add)
    }
    if (redundancy.df$no.of.results[row.num] > 1) {
      all.gene.names <- conversion.df$external_gene_name[grep(redundancy.df$unique.ids[row.num], conversion.df$uniprot_gn)]
      all.gene.names <- paste(all.gene.names, collapse = '-')
      row.to.add <- data.frame('uniprot_gn' = redundancy.df$unique.ids[row.num], 'external_gene_name' = all.gene.names)
      multiple.result.df <- rbind(multiple.result.df,row.to.add)    
    }
  }
  ### Generate a df containing missing IDs
  missing.ids <- setdiff(row.names(dataframe), conversion.df$uniprot_gn)
  missing.result.df <- data.frame('uniprot_gn' = missing.ids, 'external_gene_name' = missing.ids)
  ### Join dfs and reorder
  conversion.df <- rbind(single.result.df, multiple.result.df, missing.result.df)
  conversion.df <- conversion.df[match(row.names(dataframe),conversion.df$uniprot_gn),]
  ### Convert old IDs to new IDs
  converted.df <- dataframe
  row.names(converted.df) <- conversion.df$external_gene_name
  return(converted.df)
}

########################################################################
### Import Data

zika.df <- read.csv('Z:/Data/zika/Andrew/ZIKA-iPS_TMT10plex_Report_ARG.csv')

########################################################################
### Format

row.names(zika.df) <- zika.df$Accession
zika.df$calc..pI <- 1
zika.df <- zika.df[c(11,12,13,14,15,16,17,18,19,20)]
names(zika.df) <- c('Mock.2DPI.1' , 'Mock.2DPI.2' , 'Mock.5DPI.1' , 'Mock.5DPI.2' , 'Mock.5DPI.3' , 'Zika.2DPI.1' , 'Zika.2DPI.2' , 'Zika.5DPI.1' , 'Zika.5DPI.2' , 'Zika.5DPI.3')
zika.df <- round(zika.df*100)

column.metadata <- data.frame('Day' = c(2,2,5,5,5,2,2,5,5,5), 'Treatment' = c('mock', 'mock', 'mock', 'mock', 'mock', 'zika', 'zika', 'zika', 'zika', 'zika'), 'full.category' = c('2DPI.mock', '2DPI.mock', '5DPI.mock', '5DPI.mock', '5DPI.mock', '2DPI.zika', '2DPI.zika', '5DPI.zika', '5DPI.zika' , '5DPI.zika' ))
row.names(column.metadata) <- names(zika.df)

### Assign new names

zika.df <- uniprot_gn.to.gene.name(zika.df) # ~11 seconds

### Separate by days post infection

zika2.df <- zika.df[row.names(column.metadata[column.metadata$Day == 2 ,])]
zika5.df <- zika.df[row.names(column.metadata[column.metadata$Day == 5 ,])]

### Convert to Summerized experiment format

zika2.se <- DESeqDataSetFromMatrix(countData = as.matrix(zika2.df), colData = column.metadata[column.metadata$Day == 2 ,], design = ~ Treatment)
zika5.se <- DESeqDataSetFromMatrix(countData = as.matrix(zika5.df), colData = column.metadata[column.metadata$Day == 5 ,], design = ~ Treatment)

########################################################################
### Calc Diff. Ex. with DESeq2

zika2.de <- DESeq(zika2.se)
res2 <- results(zika2.de)
sum(res2$padj < 0.1, na.rm=TRUE)
#res2 <- results(zika2.de, alpha = 0.1)
res2 <- res2[order(res2$padj),]
res2.df <- as.data.frame(subset(res2, padj < 0.1))

zika5.de <- DESeq(zika5.se)
res5 <- results(zika5.de)
sum(res5$padj < 0.01, na.rm=TRUE)
#res5 <- results(zika5.de, alpha = 0.01)
res5 <- res5[order(res5$padj),]
res5.df <- as.data.frame(subset(res5, padj < 0.01))

########################################################################
### Filter columns based on DE results

res5.names <- row.names(row.names(res5.df))
length(row.names(res2.df))
length(row.names(res5.df))
overlap.genes <- intersect(row.names(res2.df),row.names(res5.df))
length(overlap.genes)

### Subset original data by main genes

overlap.df <- zika.df[overlap.genes,]

### Reshape data

reshaped.data <- data.frame(t(overlap.df[6,]))
reshaped.data <- cbind(reshaped.data, column.metadata)
print(names(reshaped.data[1]))

### Plot

ggplot(data = reshaped.data, aes(x = Day, y = HIBADH, color = Treatment)) +
  geom_point(size = 5) +
  labs(title = 'Expression in Zika and Mock by Day') +
  scale_color_manual(values = c('blue','red')) +
  xlim(c(1.5,5.5)) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(color = "black", fill = NA)) 




########################################################################
### Calc Diff. Ex. Manually

zika.avgs.df <- data.frame('Mock.2DPI' = apply(zika.df[1:2],1,mean),
                           'Mock.5DPI' = apply(zika.df[3:5],1,mean),
                           'Zika.2DPI' = apply(zika.df[6:7],1,mean),
                           'Zika.5DPI' = apply(zika.df[8:10],1,mean),
                           'M2.SD' = apply(zika.df[1:2],1,sd),
                           'M5.SD' = apply(zika.df[3:5],1,sd),
                           'Z2.SD' = apply(zika.df[6:7],1,sd),
                           'Z5.SD' = apply(zika.df[8:10],1,sd)
                           )
row.names(zika.avgs.df) <- row.names(zika.df)

### zika over mock, 2 DPI

zika.2DPI <- zika.avgs.df[c(3,1,7,5)]
zika.2DPI$z2.over.m2 <- zika.2DPI[1]/zika.2DPI[2]
zika.2DPI <- zika.2DPI[order(zika.2DPI$z2.over.m2, decreasing = TRUE),]
rows.to.keep <- zika.2DPI$z2.over.m2 > 2 | zika.2DPI$z2.over.m2 < 0.75
zika.2DPI.filtered <- zika.2DPI[rows.to.keep,]

### zika over mock, 5 DPI

zika.5DPI <- zika.avgs.df[c(4,2,8,6)]
zika.5DPI$z5.over.m5 <- zika.5DPI[1]/zika.5DPI[2]
zika.5DPI <- zika.5DPI[order(zika.5DPI$z5.over.m5, decreasing = TRUE),]
rows.to.keep <- zika.5DPI$z5.over.m5 >2 | zika.5DPI$z5.over.m5 < 0.5
zika.5DPI.filtered <- zika.5DPI[rows.to.keep,]

########################################################################
### Plot heat map

### 2DPI

### Filter and order all samples

zika.2 <- log2(zika2.df[row.names(zika.2DPI.filtered),]/100)
zika.2$gene <- row.names(zika.2)
#zika.2 <- data.frame(gene = row.names(zika.2DPI.filtered), zika.2DPI = log2(zika.2DPI.filtered$Zika.2DPI), mock2DPI = log2(zika.2DPI.filtered$Mock.2DPI))
zika.2 <- zika.2[rev(row.names(zika.2)),]
#zika.2$gene <- levels(factor(zika.2$gene))
zika.2$gene <- factor(zika.2$gene, levels = zika.2$gene, ordered = TRUE)
zika.2.m <- melt(zika.2)

### Plot
ggplot(zika.2.m, aes(x = variable, y = gene)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colors = c('blue', 'white', 'red'), limits = c(-3.3,3.3)) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(position = "top") +
  labs(title = 'Differential Expression, 2 DPI') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(color = "black", fill = NA)) 



### 5DPI
zika.5 <- log2(zika5.df[row.names(zika.5DPI.filtered),]/100)
zika.5$gene <- row.names(zika.5)
zika.5$gene <- factor(zika.5$gene, levels = zika.5$gene, ordered = TRUE)
zika.5 <- zika.5[rev(row.names(zika.5)),]
zika.5.m <- melt(zika.5)
limits <- c(-max(abs(zika.5.m$value)), max(abs(zika.5.m$value)))

### Plot
ggplot(zika.5.m, aes(x = variable, y = gene)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colors = c('blue', 'white', 'red'), limits = limits) +
  scale_y_discrete(position = "right") +
  scale_x_discrete(position = "top") +
  labs(title = 'Differential Expression, 5 DPI') +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5),
        legend.text = element_text(size = 12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(color = "black", fill = NA)) 


### Identify most disregulated genes

goi.2dpi <- zika.2$gene
goi.5dpi <- zika.5$gene
goi.full <- intersect(goi.2dpi,goi.5dpi)

zika.avgs.filtered <- log2(zika.avgs.df[goi.full,])
zika.avgs.filtered$gene <- row.names(zika.avgs.filtered)
zika.avgs.filtered <- zika.avgs.filtered[c(9,3,1,4,2)]

zika.avgs.filtered.m <- melt(zika.avgs.filtered)

ggplot(zika.avgs.filtered.m, aes(x = variable, y = gene)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = 'blue', high = 'yellow') +
  scale_y_discrete(position = "right") +
  scale_x_discrete(position = "top")






nba <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")
nba$Name <- with(nba, reorder(Name, PTS))
nba.m <- melt(nba)
nba.m <- ddply(nba.m, .(variable), transform, rescale = rescale(value))

(p <- ggplot(nba.m, aes(variable, Name)) + geom_tile(aes(fill = rescale), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue"))

########################################################################
### Header
########################################################################

library(gplots)
library(RColorBrewer)
library(biomaRt)
library(DESeq2)
library(Hmisc)
library(grid)
library(gridExtra)

library("mzR")
library("mzID")
library("MSnID")
library("MSnbase")
library("rpx")
library("MLInterfaces")
library("pRoloc")
library("pRolocdata")
library("MSGFplus")
library("rols")
library("hpar")

ensembl = useMart(host="www.ensembl.org",dataset="hsapiens_gene_ensembl")
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL')
ensembl = useMart(host='www.ensembl.org',biomart='ENSEMBL_MART_ENSEMBL',dataset="hsapiens_gene_ensembl")
listMarts(host="www.ensembl.org")
listDatasets(ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

########################################################################
### Functions
########################################################################

addMedSD <- function(dataframe) {
  median <- apply(dataframe,1,median)
  sd <- apply(dataframe,1,sd)
  return(data.frame(dataframe,median,sd))  
}

sortByMed <- function(dataframe) {
  order <- order(dataframe$median,decreasing=TRUE)
  return(dataframe[order,])
}

convertIDs <- function(dataframe) {
  ensemblIDs <- c()
  for (rowName in row.names(dataframe)) {
    ensemblID <- strsplit(rowName,"\\.")[[1]][1]
    ensemblIDs <- c(ensemblIDs,ensemblID)
  }
  row.names(dataframe)<-ensemblIDs
  return(dataframe)
}
addGene <- function(dataframe) {
  genes <- getBM(attributes=c('ensembl_gene_id','external_gene_name'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  genes <- genes[match(row.names(dataframe),genes[,1]),]
  Gene <- c()
  for (rowNumber in 1:length(genes[,1])) {
    newGene <- genes[rowNumber,][,2]
    Gene <- c(Gene, newGene)
  }
  dataframe[length(dataframe)+1] <- Gene
  names(dataframe)[ncol(dataframe)] <- "Gene"
  return(dataframe)
}
generate.bargraphs <- function(dataframe) {
  plot.list <- list()
  for(gene.column in 1:(ncol(dataframe)-2)){
    current.plot <- ggplot(data = dataframe, 
                           aes_string(x = 'Sample', y = names(dataframe)[gene.column], fill = "Type")) +
      geom_bar(stat = 'identity') +
      labs(title = names(dataframe[gene.column])) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x = element_blank(),
            axis.title.y = element_blank(), legend.position = 'none',
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            panel.border = element_rect(color = "black", fill = NA)) +
      scale_fill_brewer(palette = 'Set1')
    plot.list[[length(plot.list)+1]] <- current.plot
  }
  print(paste("List contains",length(plot.list),"plots"))
  return(plot.list)
}
generate.boxplots <- function(dataframe) {
  plot.list <- list()
  for(gene.column in 1:(ncol(dataframe)-2)){
    current.plot <- ggplot(data = dataframe, 
                           aes_string(x = 'Type', y = names(dataframe)[gene.column], fill = "Type")) +
      geom_boxplot(varwidth = FALSE) +
      scale_y_continuous(limits = c(0,NA)) + geom_jitter(width = 0, size = 2) +
      labs(title = names(dataframe[gene.column])) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x = element_blank(),
            axis.title.y = element_blank(), legend.position = 'none',
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            panel.border = element_rect(color = "black", fill = NA)) +
      scale_fill_brewer(palette = 'Set1')
    plot.list[[length(plot.list)+1]] <- current.plot
  }
  print(paste("List contains",length(plot.list),"plots"))
  return(plot.list)
}
########################################################################
### Import Data
########################################################################

# Normalized refseq data in units of TPM
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)

# DESeq data
deseq.data <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/DE_gene_list after_filter/OBS_vs_CTR/OBS-iHT_vs_CTR-iHT_filtered_DE_csv_results.csv", row.names=1)

# Metadata
metaData <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)
#metaData <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/HT_plus_reference_metadata.csv", row.names=1)

# Import genes of interest
uthras.genes <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/Uthras_genes_of_interest.txt")

# Import housekeeping genes
housekeeping.genes <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/Housekeeping_genes.txt")

# Import wang primers
ht.genes_lit <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/Wang_primers.txt")

# Import hypothalamic genes at pSI 1e-4
ht.genes_0.01.df <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/ht.genes_0.01.csv")
ht.genes_0.005.df <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/ht.genes_0.005.csv")
ht.genes_0.0005.df <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/ht.genes_0.0005.csv")
ht.genes_0.0001.df <- read.csv("z:/Data/RNAseq HT neurons and tissue/Andrews_files/pSI_files/ht.genes_0.0001.csv")

# Import references
references <- read.table("c://Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct",sep="\t",header=TRUE,row.names=1)

########################################################################
### Format
########################################################################

#references <- convertIDs(references)  # Remove decimal values from Ensembl IDs
#ht.reference <- references[c(1,17)]  # Subselect from the references for just hypothalamus
#ht.reference <- ht.reference[order(ht.reference[2],decreasing = TRUE),]  # Sort hypothalamus from low to high
#genes.of.interest <- ht.reference[1:20,]  # Define the genes of interest as top hypothalamic genes
#genes.of.interest$Ensembl.ID <- row.names(genes.of.interest)

### Convert references to TPM
ht.reference$Brain...Hypothalamus <- ht.reference$Brain...Hypothalamus/sum(ht.reference$Brain...Hypothalamus)*1000000

metaData$Type <- as.character(metaData$Group)
metaData$Type[grep('Adult',metaData$Source)] <- 'aHT'
metaData$Type[intersect(grep('iHT',metaData$Group),grep('CTR',metaData$Disease))] <- 'CTR'
metaData$Type[intersect(grep('iHT',metaData$Group),grep('OBS',metaData$Disease))] <- 'OBS'

### Reorder columns
TPMdata <- TPMdata[c(7,8,9,10,12,11,1,16,18,14,15,17,4,20,19,3,2,5,6)]

### Remove empty columns, trim ID decimals
deseq.data <- deseq.data[1:16]
deseq.data <- convertIDs(deseq.data)

########################################################################
### Select comparison set
########################################################################

genes.of.interest <- uthras.genes           ;title <- "Uthra's Genes of Interest"

genes.of.interest <- housekeeping.genes     ;title <- "Housekeeping Genes"

genes.of.interest <- ht.genes_0.0001.df     ;title <- "Hypothalamus genes, pSI = 1e-4"

genes.of.interest <- ht.genes_lit           ;title <- "Selected Hypothalamus genes"

########################################################################
### Subsample rows
########################################################################

### Subsample for TPM
gene.positions <- match(genes.of.interest$Ensembl.ID,row.names(TPMdata))  # Declare the row numbers in the rnaseq data which correspond to genes in or list of genes of interest
genes.df <- TPMdata[gene.positions,]  # Make a dataframe containing just the rows of RNAseq data corresponding to genes of interest
### Rename the rows with genes
row.names(genes.df) <- genes.of.interest$Gene  # Replace the row names with the gene names
genes.df <- genes.df[complete.cases(genes.df),]

### Subsample for DEseq
gene.positions <- match(genes.of.interest$Ensembl.ID,row.names(deseq.data))  # Declare the row numbers in the rnaseq data which correspond to genes in or list of genes of interest
genes.df <- deseq.data[gene.positions,]  # Make a dataframe containing just the rows of RNAseq data corresponding to genes of interest
### Rename the rows with genes
row.names(genes.df) <- genes.of.interest$Gene  # Replace the row names with the gene names
genes.df <- genes.df[complete.cases(genes.df),]


########################################################################
### Normalize (optional)
########################################################################

dataframe <- addMedSD(genes.df)

for(column in 1:(ncol(dataframe)-2)){
  current.vector <- dataframe[,column] / dataframe$median
  correction.factor <- mean(current.vector)
  #print(correction.factor)
  dataframe[column] <- round(dataframe[column]/correction.factor,2)
}

genes.df <- dataframe[1:(ncol(dataframe)-2)]

#for (column in 1:13) {
#  barplot.df[column] <- round(barplot.df[column]*1000/barplot.df$GAPDH,3)
#}

########################################################################
### Prepare plot data
########################################################################

#genes.df$Reference <- ht.reference$Brain...Hypothalamus[1:nrow(genes.df)]
plot.df <- data.frame(t(genes.df))
plot.df$Type <- factor(metaData$Type[match(row.names(plot.df),row.names(metaData))], levels = c('aHT','CTR','OBS','iMN'))
plot.df$Sample <- factor(row.names(plot.df), levels = names(genes.df))

########################################################################
### Generate plots
########################################################################
###### Bars

bar.list <- generate.bargraphs(plot.df)

###### Boxes
box.list <- generate.boxplots(plot.df)

########################################################################
### Select data to plot
########################################################################

plot.list <- bar.list

plot.list <- box.list[13:24]

########################################################################
### Display tiled plots
########################################################################

(tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],
                              plot.list[[7]],plot.list[[8]],plot.list[[9]],plot.list[[10]],plot.list[[11]],plot.list[[12]],
                              ncol = 4, top = title))

(tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],ncol = 4, top = title))

(tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],
                              plot.list[[7]],plot.list[[8]],plot.list[[9]],
                              ncol = 3, top = title))


########################################################################
### Export plot
########################################################################

setwd("z:/Uthra/HT paper/Bioinformatics figures/Plots of gene expression/")

tiff(filename=paste0(substr(title,1,6),'_2-',strftime(Sys.time(),"%a%b%d%H%M"),".tiff"), 
     type="cairo",
     units="in", 
     width=14, 
     height=14, 
     pointsize=12, 
     res=300)
tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],
                             plot.list[[7]],plot.list[[8]],plot.list[[9]],plot.list[[10]],plot.list[[11]],plot.list[[12]],
                             ncol = 4, top = title)
dev.off()

png(filename=paste0(substr(title,1,6),'_02-',strftime(Sys.time(),"%a%b%d%H%M"),".png"), 
     type="cairo",
     units="in", 
     width=14, 
     height=14, 
     pointsize=12, 
     res=300)
tiled.figure <- grid.arrange(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],plot.list[[6]],
                             plot.list[[7]],plot.list[[8]],plot.list[[9]],plot.list[[10]],plot.list[[11]],plot.list[[12]],
                             ncol = 4, top = title)
dev.off()

########################################################################
########################################################################
### Generate Single Plots
########################################################################
########################################################################




########################################################################
### Select geneset
########################################################################

genes.of.interest <- uthras.genes           ;title <- "Selected Hypothalamus genes"

genes.of.interest <- ht.genes_lit           ;title <- "Selected Hypothalamus genes"

genes.of.interest <- ht.genes_0.0001.df     ;title <- "Hypothalamus genes, pSI = 1e-4"

########################################################################
### Subsample
########################################################################
### Subsample rows

gene.positions <- match(genes.of.interest$Ensembl.ID,row.names(TPMdata))  # Declare the row numbers in the rnaseq data which correspond to genes in or list of genes of interest
genes.df <- TPMdata[gene.positions,]  # Make a dataframe containing just the rows of RNAseq data corresponding to genes of interest
### Rename the rows with genes
row.names(genes.df) <- genes.of.interest$Gene  # Replace the row names with the gene names
genes.df <- genes.df[complete.cases(genes.df),]

genes.with.med <- addMedSD(genes.df)
genes.df <- genes.df[genes.with.med$median>5,]

  
########################################################################
### Format data
########################################################################

#genes.df$Reference <- ht.reference$Brain...Hypothalamus[1:nrow(genes.df)]
plot.df <- data.frame(t(genes.df))
plot.df$Type <- as.character(metaData$Type[match(row.names(plot.df),row.names(metaData))], levels = c('aHT','CTR','OBS','iMN'))
plot.df$Sample <- factor(row.names(plot.df), levels = names(genes.df))

### Subsample columns
plot.df <- plot.df[plot.df$Type != 'iMN',]

### Reclassify type
plot.df$Type[plot.df$Type == 'CTR'] <- 'iHT'
plot.df$Type[plot.df$Type == 'OBS'] <- 'iHT'

########################################################################
### Loop through columns and output plots
########################################################################

setwd("z:/Uthra/HT paper/Bioinformatics figures/Plots of gene expression/Individual plots/")

for(col.numb in 1:(ncol(plot.df)-2)){
  gene <- names(plot.df)[col.numb]
  png(filename=paste0(gene,'-HT-',strftime(Sys.time(),"%a%b%d%H%M"),".png"), 
      type="cairo",
      units="in", 
      width=10, 
      height=10, 
      pointsize=20, 
      res=300)
  
    plot <-ggplot(data = plot.df, 
         aes_string(x = 'Type', y = gene, fill = "Type")) +
    geom_boxplot(varwidth = FALSE) +
    scale_y_continuous(limits = c(0,NA)) + geom_jitter(width = 0, size = 2) +
    labs(title = gene) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20),axis.title.x = element_blank(),
          axis.text.y = element_text(size = 20), axis.title.y = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          panel.border = element_rect(color = "black", fill = NA)) +
    scale_fill_brewer(palette = 'Set1')
    print(plot)
  dev.off()
  print(paste(gene,'exported'))
}


########################################################################
### Export
########################################################################










####
plot.list <- list()

for (gene in 1:13){
  print(gene)
  #barplot.temp <- barplot.df[c(gene,ncol(boxplot.df)-1,ncol(boxplot.df))]
  #print(barplot.temp)
  barplot.temp <- barplot.df
  f = ggplot(data = barplot.temp, aes_string(x = "Sample", y = names(barplot.df)[gene], fill="disease")) + 
    geom_bar(stat = "identity") +
    ggtitle(names(barplot.df[gene])) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none")
  #Sys.sleep(1)
  print(f$labels$title)
  print(f$data[,1])
  plot.list[[length(plot.list)+1]] = f
  #print(barplot.df[gene][1,])
}

plot.list[[1]]
plot.list[[2]]
plot.list[[3]]
plot.list[[4]]
plot.list[[5]]
plot.list[[6]]
plot.list[[7]]
plot.list[[8]]
plot.list[[9]]
plot.list[[10]]
plot.list[[11]]
plot.list[[12]]
plot.list[[13]]





