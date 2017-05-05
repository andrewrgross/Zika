######################################################################################
###### Proteomic analysis for networks, 2016-02-17
######################################################################################
### Functions

convert.char.column.to.num <- function(dataframe,column) {
  new.vector <- c()
  for (value in dataframe[,column]) {
    new.value <- as.numeric(strsplit(as.character(value),"/")[[1]][1])
    new.vector <- c(new.vector,new.value)
  }
  dataframe[column] <- new.vector
  return(dataframe)
}
######################################################################################
### Upload data

#setwd(dir = "Z:/Uthra/HT paper/Bioinformatics figures/IPA analysis/network_files_for_cytoscape/")
setwd(dir = "Z:/Data/zika/Andrew/IPA results/transcription_factor_figure/")

#networks <- read.csv("z:/Data/C9-ALS MNs and bioenergetics/Mass Spec/Dhruv Sareen/Report 2015.02/Andrew/IPA_pathways_V2-1.csv")
#networks <- read.csv("Features_disease.csv")
networks <- read.csv("Features_pathways.csv")

#gene.fc <- read.csv('Z:/Data/C9-ALS MNs and bioenergetics/Mass Spec/Dhruv Sareen/Report 2015.02/Andrew/Gene_fold_changes.csv')    # Import genes and their fold changes
gene.fc <- read.csv('5DPI ZKVvsMOCK-all_pep_IPA.csv')

#upstream <- read.csv('Z:/Data/C9-ALS MNs and bioenergetics/Mass Spec/Dhruv Sareen/Report 2015.02/Andrew/IPA_upstream_regulators-SELECT_FOR_FIG_V2-1.csv')
upstream <- read.csv('Upstream_regulators_B.csv')

######################################################################################
### Format new columns

### Adjust files
#networks <- networks[-5]                                                     # Remove unwanted column from network table
length <- nrow(networks)                                                     # Declare number of rows in network table [DEPRECIATED?]
upstream <- upstream[order(upstream$p.value.of.overlap),]                    # Reorder upstream regulators
#upstream <- upstream[1:16,]                                                  # Filter upstream regulators
upstream <- upstream[-2]

### Drop unwanted rows
#networks <- networks[-c(6,7,9),]

######################################################################################
### Generate gene and pathway nodes and blank network table

### Define empty data frames
output.gene.node <- data.frame(sources = character(), label = character(), group = character(), size = double(), fc = double(), type = character(), shape = character(), border.color = character(), stringsAsFactors=FALSE)
output.gene.net <- data.frame(sources=character(),target=character(),interaction=character(),boolean=character(), string=character(),shared=numeric(), edge.color = character(), stringsAsFactors=FALSE)
full.gene.list <- c('GENES!')

### Add each pathway and gene
for(pathway.num in 1:length) {                              # pathway.num = 1 (Troubleshooting)
  genes <- as.character(networks$Molecules[pathway.num])                      # Call the genes from the current pathway
  genes <- strsplit(genes,',')[[1]]                                           # Break the list into individual genes
  pathway = as.character(networks$Ingenuity.Canonical.Pathways[pathway.num])  # Declare the name of the current pathway
  total.genes <- length(genes)                             # Define the number of possible genes in the pathway
  new.node.row <- data.frame(sources = as.character(pathway), label = as.character(pathway), group=as.character("PATHWAY"), size = (total.genes*10)+20, fc = -4, type = 'Pathway', shape = 'Round Rectangle', border.color = '#000000', stringsAsFactors = FALSE)
  output.gene.node <- rbind(output.gene.node, new.node.row)                   # Add new pathway node to node table
  new.net.row <- data.frame(sources = as.character(pathway), target = as.character(''), interaction = 'cooccurrence',boolean = "TRUE",string = "ABC", shared = 0, edge.color = '', stringsAsFactors = FALSE)
  output.gene.net <- rbind(output.gene.net, new.net.row)                      # Add new pathway net row to table
  for (gene in genes) {                                                          # Loop through all genes in a pathway
    gene.row <- grep(gene,gene.fc$Symbol)                                        # Identify the row of the current gene in the fold change table
    fold.change <- gene.fc$Expr.Log.Ratio[gene.row]                              # Declare the fold change of the current gene
    gene.type <- as.character(gene.fc$Type.s.[gene.row])                         # Declare the gene type
    shape <- ''
    boarder.color <- ''
    source.name <- gene                                                          # Identify gene as the source
    target <- ''                                                                 # Assume the gene has no connections.  Could I omit this?
    count <- 0                                                                 # Declare count to be zero
    number.of.duplicates <- sum(grepl(gene, full.gene.list))
    if(number.of.duplicates > 0){                        # Check if the source name is in the full.gene.list
      target <- source.name                                                         # If a gene has been seen before, declare it as a target
      roman.numerals <- paste0(rep('i',number.of.duplicates), collapse = '')
      source.name <- paste0(gene, '_', roman.numerals, collapse = '')                                       # Modify the name to distinguish it from the previous one
      print(source.name)
      count <- 1                                                                  # Declare count = 1
    }
    new.node.row <- data.frame(sources = as.character(source.name), label= as.character(gene), group=as.character(pathway), size = 10, fc = fold.change, type = gene.type, shape = '', border.color = '', stringsAsFactors = FALSE)
    output.gene.node <- rbind(output.gene.node, new.node.row)
    new.net.row <- data.frame(sources = as.character(source.name), target = as.character(target), interaction = 'cooccurrence', boolean = "TRUE",string = "ABC",shared = count, edge.color = '#ffffff', stringsAsFactors = FALSE)
    output.gene.net <- rbind(output.gene.net, new.net.row)
    #print(data.frame(sources = as.character(source.name, target = as.character(target))))
    full.gene.list <- c(full.gene.list,gene)
  }
}


######################################################################################
### Generate upstream regulator node table

upstream.gene.node <- data.frame(sources = character(), label = character(), group = character(), size = double(), fc = double(), type = character(), shape = character(), border.color = character(), stringsAsFactors=FALSE)
upstream.gene.net <- data.frame(sources=character(),target=character(),interaction=character(),boolean=character(), string=character(),shared=numeric(), edge.color = character(), stringsAsFactors=FALSE)

for (row.num in 1:nrow(upstream)) {               # row.num = 1
  row <- upstream[row.num,]
  key <- row[,1]
  type <- row[,2]
  fc <- row[,4]
  pval <- row[,5]
  new.row <- data.frame(sources = key, label = key, group = 'upstream', size = 1, fc = fc, type = type, shape = '', border.color = '', stringsAsFactors = FALSE)
  upstream.gene.node <- rbind(upstream.gene.node, new.row)
  ### Loop through downstream genes to add network connections
  genes <- as.character(row$Target.molecules.in.dataset)                      # Call the genes from the current pathway
  genes <- strsplit(genes,',')[[1]]                                           # Break the list into individual genes
  genes <- intersect(genes, output.gene.node$label)
  for (gene in genes) {           # gene = genes[1]
    targets <- grep(gene,output.gene.node$label)
    #print(length(targets))
    pathways <- as.character(output.gene.node$group[targets])
    targets <- as.character(output.gene.node$sources[targets])
    for(target.pos in 1:length(targets)) {
      target <- targets[target.pos]
      pathway <- pathways[target.pos]
      shared <- length(unique(pathways))
      #shared <- length(targets)
      new.net.row <- data.frame(sources = key, target = target, interaction = 'regulates',boolean = "TRUE",string = pathway, shared = shared, edge.color = '', stringsAsFactors = FALSE)
      upstream.gene.net <- rbind(upstream.gene.net, new.net.row)
    }
  }
}

######################################################################################
### Append upstream regulator nodes to full node table

output.gene.node <- rbind(output.gene.node, upstream.gene.node)
output.gene.net <- rbind(output.gene.net, upstream.gene.net)

### TROUBLE SHOOTING: This breaks it
output.gene.net$edge.color <- as.character(output.gene.net$edge.color)
#output.gene.node[7:8] <- as.character(output.gene.node[7:8])

######################################################################################
### Add nodes for legend

sources <- c('Pathway', 'Transcription Factor', 'Enzyme', 'Kinase', 'Transporter', 'Transmembrae Receptor', 'Ligand-dependent Nuclear Receptor', 'Other Protein', 'Chemical')
label <- c('Pathway', 'Transcription Factor', 'Enzyme', 'Kinase', 'Transporter', 'Transmembrae Receptor', 'Ligand-dependent Nuclear Receptor', 'Other Protein', 'Chemical')
group <- rep('Legend', 9)
size <- c(110, 10, 10, 10, 10, 10, 10, 10, 10)
fc <- rep(0.65, 9)
type <- sources
legend.gene.node <- data.frame(sources, label, group, size, fc, type, stringsAsFactors = FALSE)

######################################################################################
### Assign edge colors

save.net <- output.gene.net
output.gene.net <- save.net
save.node <- output.gene.node
output.gene.node <- save.node

### Edge color conversion table
colors <- c('#d53e4f','#f46d43', '#fdae61', '#fee08b','#ffffbf','#e6f598','#abdda4','#66c2a5','#3288bd')
edge.color.table <- data.frame(networks$Ingenuity.Canonical.Pathways, colors)

for(row.pos in 1:nrow(output.gene.net)) {
  pathway <- as.character(output.gene.net[,5][row.pos])
  edge.color.table.row <- match(pathway, edge.color.table[,1])
  color <- as.character(edge.color.table$colors[edge.color.table.row])
  output.gene.net[,7][row.pos] <- color
}

######################################################################################
### Assign border colors

### Border color conversion table
gene.types <- as.character(unique(output.gene.node$type))
border.colors <- c('#000000','#3F1669','#1F6060','#7A135A','#336699','#505050','#996600','#969696','#7A1327')
shapes <- c('Round Rectangle', 'Triangle', 'Diamond', 'Rectangle', 'Diamond', 'Ellipse', 'v', 'Hexagon', 'V')
node.shape.color.table <- data.frame(gene.types, border.colors, shapes)

### Specify shapes and border colors
for(row.pos in 1:nrow(output.gene.node)) {
  gene.type <- as.character(output.gene.node[,6][row.pos])
  shape.color.row <- match(gene.type, node.shape.color.table[,1])
  border.color <- as.character(node.shape.color.table$border.colors[shape.color.row])
  shape <- as.character(node.shape.color.table$shapes[shape.color.row])
  output.gene.node[,7][row.pos] <- shape
  output.gene.node[,8][row.pos] <- border.color
}


######################################################################################
### Output shared genes

#setwd(dir = "Z:/Data/C9-ALS MNs and bioenergetics/Mass Spec/Dhruv Sareen/Report 2015.02/Andrew/network files for cytoscape/")
#names(outputNetwork) <- c("source","target","interaction","boolean attribute","string attribute","floating point attribute")
#write.table(output.node,paste0("NODE--C9-pathways_for_cytoscape_",substr(weekdays(Sys.Date()),1,4),"-",format(Sys.Date(),"%b-%d"),".txt"),row.names=FALSE,sep="\t",quote=FALSE)
#write.table(output.net,paste0("NET--C9-pathways-for_cytoscape_",substr(weekdays(Sys.Date()),1,4),"-",format(Sys.Date(),"%b-%d"),".txt"),row.names=FALSE,sep="\t",quote=FALSE)

### Output alternate plot files
write.table(output.gene.node, paste0("GENE-NODES--Zika_cytoscape_",substr(weekdays(Sys.Date()),1,4),"-",format(Sys.Date(),"%b-%d"),".txt"),row.names=FALSE,sep="\t",quote=FALSE)
write.table(output.gene.net, paste0("GENE-NET--Zika_cytoscape_",substr(weekdays(Sys.Date()),1,4),"-",format(Sys.Date(),"%b-%d"),".txt"),row.names=FALSE,sep="\t",quote=FALSE)



### Scratchwork

### The old way
for(row.num in 1:nrow(networks)) {
  row <- networks[row.num,]
  pathway <- as.character(unlist(row[1]))
  size <- sum(row[4:7])
  fc <- row[,3]
  pval <- row[,2]
  new.row <- data.frame(key = pathway, label = pathway, type = 'pathway', group = 'pathway', size = size, fc = 0, pval = pval)
  output.node <- rbind(output.node,new.row)
  genes <- as.character(unlist(row[8]))
  genes <- strsplit(genes,',')[[1]]
  for(gene in genes) {
    key <- paste0(substr(pathway,1,3),'-',gene)
    label <- gene
    gene.row <- grep(gene,gene.fc$Symbol)
    fc <- gene.fc$Expr.Log.Ratio[gene.row]
    type <- as.character(gene.fc$Type.s.[gene.row])
    new.row <- data.frame(key = key, label = label, type = type, group = pathway, size = 1, fc = fc, pval = 0)
    output.node <- rbind(output.node,new.row) 
  }
}

######################################################################################
### Format alternate gene plot

### Generate an empty data frame for holding a list of genes and groups
output.gene.node <- data.frame(sources = character(), label = character(), group = character(), size = double(), fc = double(), type = character(), stringsAsFactors=FALSE)
output.gene.net <- data.frame(sources=character(),target=character(),interaction=character(),boolean=character(), string=character(),shared=numeric(), stringsAsFactors=FALSE)

full.gene.list <- c('GENES!')

### Loop through the pathways and add genes
for(pathway.num in 1:length) {
  genes <- as.character(networks$Molecules[pathway.num])
  genes <- strsplit(genes,',')[[1]]
  pathway = as.character(networks$Ingenuity.Canonical.Pathways[pathway.num])
  total.genes <- sum(networks[pathway.num,][4:7])
  new.node.row <- data.frame(sources = as.character(pathway), label = as.character(pathway), group=as.character("PATHWAY"), size = total.genes, fc = 0, type = 'pathway')
  output.gene.node <- rbind(output.gene.node, new.node.row)                   # Add new pathway node to node table
  #output.gene.node <- rbind(output.gene.node,data.frame(sources = as.character(pathway), label = as.character(pathway), group=as.character("PATHWAY"), size = total.genes, fc = 0, type = 'pathway'))
  new.net.row <- data.frame(sources = as.character(pathway), target = as.character(''), interaction = 'cooccurrence',boolean = "TRUE",string = "ABC", shared = 0)
  output.gene.net <- rbind(output.gene.net, new.net.row)                      # Add new pathway net row to table
  #output.gene.net <- rbind(output.gene.net, data.frame(sources = as.character(pathway), target = as.character(''), interaction = 'cooccurrence',boolean = "TRUE",string = "ABC", shared = 0))
  for (gene in genes) {                                      # Loop through all genes in a pathway
    gene.row <- grep(gene,gene.fc$Symbol)
    fold.change <- gene.fc$Expr.Log.Ratio[gene.row]
    gene.type <- as.character(gene.fc$Type.s.[gene.row])
    source.name <- gene
    target <- ''
    count <- 0
    if(TRUE %in% grepl(source.name,full.gene.list) == TRUE){
      target <- source.name
      source.name <- paste0(source.name,'_i')
      count <- 1
      }
    output.gene.node <- rbind(output.gene.node,data.frame(sources = as.character(source.name), label= as.character(gene), group=as.character(pathway), size = 1, fc = fold.change, type = gene.type))
    output.gene.net <- rbind(output.gene.net, data.frame(sources = as.character(source.name), target = as.character(target), interaction = 'cooccurrence', boolean = "TRUE",string ="ABC",shared = count))
    #print(data.frame(sources = as.character(source.name, target = as.character(target))))
    full.gene.list <- c(full.gene.list,source.name)
  }
}
######################################################################################
### Calculate shared genes

outputNetwork <- data.frame(sources=character(),target=character(),
                            interaction=character(),boolean=character(),
                            string=character(),shared=double(),
                            stringsAsFactors=FALSE)

for (i in 1:(length-1)) {                                           # Run through each node in turn
  sourceName <- toString(networks$Ingenuity.Canonical.Pathways[i])        # Define current source
  molecules <- networks$Molecules[i]                            # Generate a string of molecules in the current node
  moleculesA <- unlist(strsplit(toString(molecules),","))       # Convert string to list
  for (j in 1:(length-i)) {                                     # Run through each of the nodes following the current one
    targetName <- toString(networks$Ingenuity.Canonical.Pathways[i+j])      # Define current source
    molecules <- networks$Molecules[i+j]                        # Generate a string of molecules
    moleculesB <- unlist(strsplit(toString(molecules),","))     # Convert string to list
    count <- length(intersect(moleculesA,moleculesB))           # Calculate the number of genes shared between the nodes
    row <- c(sourceName,targetName,"cooccurrence","TRUE","ABC",count)
    outputNetwork[length(outputNetwork[,1])+1,] <- row
  }
}

outputNode <- networks[1:7]
names(outputNode) <- c("sources","pVal","ratio","Downregulated","No.change","Upregulated","No.overlap")
