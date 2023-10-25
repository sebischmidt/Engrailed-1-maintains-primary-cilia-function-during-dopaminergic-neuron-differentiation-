######Load packages######
library(biomaRt)
library(org.Hs.eg.db)
library(RCy3)
library(RColorBrewer)
library(xlsx)
#install.packages('writexl')
cytoscapePing()  #this will tell you if you're able to successfully connect to Cytoscape or not
#installApp('WikiPathways') 
#installApp('CyTargetLinker') 
#installApp('stringApp') 



######Load data######
setwd('C:/Users/sschm/OneDrive - Helmholtz Zentrum München/Manuscript Engrailed/sina_transcriptome/R_En1ko')


unstim <- read.delim("results/tables/DEGs_unstim.txt", sep = "\t")
stim <- read.delim("results/tables/DEGs_stim.txt", sep = "\t")



######Prepare data######

genes <- as.character(unstim$symbol)
mart = useDataset("hsapiens_gene_ensembl", useEnsembl(biomart="ensembl", version=98))
G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),values=genes,mart=mart)
unstim$ENTREZID = G_list[match(genes, G_list$hgnc_symbol),]$entrezgene_id
unstim$ensemblid = G_list[match(genes, G_list$hgnc_symbol),]$ensembl_gene_id
colnames(unstim)[7] <- "SYMBOL"
colnames(unstim)[2] <- "log2fc"

genes <- as.character(stim$symbol)
G_list <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),values=genes,mart=mart)
stim$ENTREZID = G_list[match(genes, G_list$hgnc_symbol),]$entrezgene_id
stim$ensemblid = G_list[match(genes, G_list$hgnc_symbol),]$ensembl_gene_id
colnames(stim)[7] <- "SYMBOL"
colnames(stim)[2] <- "log2fc"


######Plot figures in Cytoscape - unstim######
##Wnt signaling pathway
RCy3::commandsRun('wikipathways import-as-pathway id=WP363') 
loadTableData(unstim, data.key.column = "ENTREZID", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/unstim_Wnt signaling pathway.pdf", type='PDF')

##Wnt signaling
RCy3::commandsRun('wikipathways import-as-pathway id=WP428') 
loadTableData(unstim, data.key.column = "ENTREZID", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/unstim_Wnt signaling.pdf", type='PDF')

##BDNF signaling
RCy3::commandsRun('wikipathways import-as-pathway id=WP2380') 
loadTableData(unstim, data.key.column = "ENTREZID", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/unstim_BDNF signaling.pdf", type='PDF')

##Neurodegeneration
RCy3::commandsRun('wikipathways import-as-pathway id=WP4577') 
loadTableData(unstim, data.key.column = "ensemblid", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/unstim_neurodegeneration.pdf", type='PDF')

##G-protein signaling
RCy3::commandsRun('wikipathways import-as-pathway id=WP35') 
loadTableData(unstim, data.key.column = "ENTREZID", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/unstim_G-protein signaling.pdf", type='PDF')

##Ciliary landscape
RCy3::commandsRun('wikipathways import-as-pathway id=WP4352') 
loadTableData(unstim, data.key.column = "ensemblid", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/unstim_Ciliary landscape.pdf", type='PDF')

##Ciliopathies
RCy3::commandsRun('wikipathways import-as-pathway id=WP4803') 
loadTableData(unstim, data.key.column = "ensemblid", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/unstim_Ciliopathies.pdf", type='PDF')


##Cilia development
RCy3::commandsRun('wikipathways import-as-pathway id=WP4536') 
loadTableData(unstim, data.key.column = "ENTREZID", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/unstim_Cilia development.pdf", type='PDF')


##one-carbon metabolism
RCy3::commandsRun('wikipathways import-as-pathway id=WP241') 
loadTableData(unstim, data.key.column = "ENTREZID", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/unstim_one-carbon metabolism.pdf", type='PDF')


######Plot figures in Cytoscape - stim######
##Wnt signaling pathway
RCy3::commandsRun('wikipathways import-as-pathway id=WP363') 
loadTableData(stim, data.key.column = "ENTREZID", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/stim_Wnt signaling pathway.pdf", type='PDF')

##Wnt signaling
RCy3::commandsRun('wikipathways import-as-pathway id=WP428') 
loadTableData(stim, data.key.column = "ENTREZID", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/stim_Wnt signaling.pdf", type='PDF')

##BDNF signaling
RCy3::commandsRun('wikipathways import-as-pathway id=WP2380') 
loadTableData(stim, data.key.column = "ENTREZID", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/stim_BDNF signaling.pdf", type='PDF')

##Neurodegeneration
RCy3::commandsRun('wikipathways import-as-pathway id=WP4577') 
loadTableData(stim, data.key.column = "ensemblid", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/stim_neurodegeneration.pdf", type='PDF')

##G-protein signaling
RCy3::commandsRun('wikipathways import-as-pathway id=WP35') 
loadTableData(stim, data.key.column = "ENTREZID", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/stim_G-protein signaling.pdf", type='PDF')

##Ciliary landscape
RCy3::commandsRun('wikipathways import-as-pathway id=WP4352') 
loadTableData(stim, data.key.column = "ensemblid", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/stim_Ciliary landscape.pdf", type='PDF')

##Ciliopathies
RCy3::commandsRun('wikipathways import-as-pathway id=WP4803') 
loadTableData(stim, data.key.column = "ensemblid", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/stim_Ciliopathies.pdf", type='PDF')


##Cilia development
RCy3::commandsRun('wikipathways import-as-pathway id=WP4536') 
loadTableData(stim, data.key.column = "ENTREZID", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/stim_Cilia development.pdf", type='PDF')


##one-carbon metabolism
RCy3::commandsRun('wikipathways import-as-pathway id=WP241') 
loadTableData(stim, data.key.column = "ENTREZID", table.key.column = "XrefId")
data.values = c(-6,-0.3,-0.2,-0.1,-0.05,0,0.05,0.1,0.2,0.3,6)
node.colors <- c(rev(brewer.pal(length(data.values), "RdBu")))
setNodeColorMapping("log2fc", data.values, node.colors, default.color = "#FFFFFF", style.name = "WikiPathways")
exportImage(filename="results/cytoscape/stim_one-carbon metabolism.pdf", type='PDF')



