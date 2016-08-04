### RNA Seq preliminary analysis, Andrew R Gross, 2016-05-16
### This script is intended to upload normalized expression data and plot a variety of features in order to make general assessments of the data

########################################################################
### Header
########################################################################

library(gplots)
library(RColorBrewer)
library(biomaRt)
library(DESeq2)
library(Hmisc)
library(ggbiplot)
library(rgl)

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

addDescription <- function(dataframe) {
  descriptions <- getBM(attributes=c('ensembl_gene_id','description'), filters='ensembl_gene_id', values=row.names(dataframe), mart=ensembl)
  #print(descriptions)
  descriptions <- descriptions[match(row.names(dataframe),descriptions[,1]),]
  Descrip <- c()
  for (rowNumber in 1:length(descriptions[,1])) {
    fullDescr <- descriptions[rowNumber,][,2]
    shortDescr <- strsplit(fullDescr," \\[Source")[[1]][1]
    Descrip <- c(Descrip, shortDescr)
  }
  #print(Descrip)
  dataframe[length(dataframe)+1] <- Descrip
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
  return(dataframe)
}

########################################################################
### Import Data
########################################################################

# Metadata
metadata.sam <- read.csv("z://Data/RNAseq HT neurons and tissue/2nd rerun/samples.csv", row.names=1)
metadata.ref <- read.csv("c://Users/grossar/Bioinform/DATA/rna_seq/reference_metadata.csv",row.names=1)
metadata.all <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/HT_plus_reference_metadata.csv",row.names=1)

# Normalized
TPMdata <- read.csv("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.tpm.csv", row.names=1)

# Raw, Non-normalized
counts <- read.table("z://Data/RNAseq HT neurons and tissue/Andrews_files/20160317_Sareen_rerun-29299281.all.count", row.names=1, header=TRUE)

# References
references <- read.table("c://Users/grossar/Bioinform/DATA/rna_seq/Reference_transcriptomes/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct",sep="\t",header=TRUE,row.names=1)

# Example data
#dds <- DESeqDataSet(se, design = ~ cell + dex)

########################################################################
### Formating
########################################################################

TPMdata <- convertIDs(TPMdata)

### Rename rows and columns
names(counts) <- row.names(metadata.sam)                    # Correct column names
counts <- convertIDs(counts)                            # Correct row names

### Remove unwanted row, correct row names
references <- references[2:length(references)]
references <- convertIDs(references)

### Rescale RPKM to TPM
for (columnNumber in 1:length(references)) {
  references[columnNumber] <- references[columnNumber]*1000000/sum(references[columnNumber])
}

### Round values to two decimal places
references <- round(references,2)


########################################################################
### Joining 
########################################################################

### Generate list of all ensemble IDs in each of the two sets
referenceIDs <- row.names(references)
sampleIDs <- row.names(TPMdata)

### Genereate list of the IDs not present in each
notSampleIDs <- setdiff(referenceIDs,sampleIDs)
notReferenceIDs <- setdiff(sampleIDs,referenceIDs)

### Create empty rows and append them to sample DF
sampleAddDF <- data.frame(mat.or.vec(length(notSampleIDs),length(TPMdata)))
names(sampleAddDF) <- names (TPMdata)
row.names(sampleAddDF) <- notSampleIDs
TPMdata.plus <- rbind(TPMdata,sampleAddDF)

### Create empty rows and append them to reference dataframe
referenceAddDF <- data.frame(mat.or.vec(length(notReferenceIDs),length(references)))
names(referenceAddDF) <- names(references)
row.names(referenceAddDF) <- notReferenceIDs
reference.plus <- rbind(references,referenceAddDF)

### Reorder sample DF
samplesDF <- TPMdata.plus[order(row.names(TPMdata.plus)),]
referencesDF <- reference.plus[order(row.names(reference.plus)),]

### Append references to sample DF
sam.plus.ref.df <- cbind(samplesDF,referencesDF)

########################################################################
### Subsetting
########################################################################

### Normalized
DF.norm.iHT <- TPMdata[grep("iHT",metadata.sam$Group)]
DF.norm.HT <- TPMdata[grep("HT",metadata.sam$Type)]
DF.norm.iPSC <- TPMdata[grep("iPSC",metadata.sam$Source)]
DF.norm.woref <- TPMdata[intersect(grep("M",metadata.sam$Sex),grep("F",metadata.sam$Sex))]
DF.norm <- TPMdata

### Raw
DF.raw.iHT <- counts[grep("iHT",metadata.sam$Group)]
DF.raw.HT <- counts[grep("HT",metadata.sam$Type)]
DF.raw.iHT.fem <- counts[intersect(grep("iHT",metadata.sam$Group),grep("F",metadata.sam$Sex))]
DF.raw.HT.fem <- counts[intersect(grep("HT",metadata.sam$Type),grep("F",metadata.sam$Sex))]
DF.raw.iPSC <- counts[grep("iPSC",metadata.sam$Source)]
DF.raw.woref <- counts[intersect(grep("M",metadata.sam$Sex),grep("F",metadata.sam$Sex))]
DF.raw <- counts

### Refernces
DF.references <- references
DF.ref.brain <- references[grep("Brain",metadata.ref$Class)]

### Both
DF.both.brain <- sam.plus.ref.df[grep("Brain",metadata.all$Organ)]

########################################################################
### Filtering
########################################################################

metadata.sam.iHT <- metadata.sam[metadata.sam$Group == 'iHT',]
metadata.sam.iHT <- metadata.sam.iHT[c(1,2,3,5,6,7,8,9,10,11),]
countMatrix <- as.matrix(iHTcounts[match(row.names(metadata.sam.iHT),names(iHTcounts))])

metadata.sam.iHT.f <- metadata.sam.iHT[grep("F",metadata.sam.iHT$Sex),]
metadata.sam.iHT.f <- metadata.sam.iHT.f[c(1,2,4,5,6,7),]
countMatrix <- as.matrix(iHTcounts[match(row.names(metadata.sam.iHT.f),names(iHTcounts))])

metadata.sam.HT <- metadata.sam[metadata.sam$Type == "HT",]
countMatrix <- as.matrix(HTcounts[match(row.names(metadata.sam.HT),names(HTcounts))])

countMatrix <- as.matrix(references)
countMatrix <- round(countMatrix)
head(countMatrix)

metadata.sam.selected <- metadata.sam.HT
#metadata.sam.selected <- metadata.sam.iHT.f


########################################################################
### Generate DESeq dataset object
########################################################################

### Generate DESeq object from HT counts
deHT <- DESeqDataSetFromMatrix(countData = countMatrix, colData = metadata.sam.selected, design = ~ Source + Sex + Disease)
nrow(deHT)                                  # Count rows in deHT
deHT <- deHT[rowSums(counts(deHT)) > 1, ]   # Retain only rows with at least one count
nrow(deHT)                                  # Count rows again
deHTtrans <- rlog(deHT,blind=FALSE)         # Perform rlog transformation



### Examine transformation effect
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(deHT)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1), pch=16, cex=0.3)
plot(assay(deHTtrans)[,1:2], pch=16, cex=0.3)

######################################################################################
### Perform PCA
######################################################################################

### Plot PCA using DESeq method
plotPCA(deHTtrans, intgroup= c("Source", "Disease", "Sex"))

### Collect plot details
(data <- plotPCA(deHTtrans, intgroup = c("Source", "Disease", "Sex"), returnData=TRUE))
percentVar <- round(100 * attr(data, "percentVar"))

### Plot manually
ggplot(data, aes(PC1, PC2, color=Source, shape=Disease)) + geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  geom_text(aes(label=row.names(data)),hjust=0,vjust=0)

### Plot PCA pc method

nrow(countMatrix)                                  # Count rows
countMatrix <- countMatrix[rowSums(countMatrix) > 1, ]   # Retain only rows with at least one count
nrow(countMatrix)                                  # Count rows again
countMatrix.rlog <- rlog(countMatrix,blind=FALSE)  # Perform rlog transformation
count.trans <- t(countMatrix.rlog)                 # Transpose matrix to place samples as rows and genes as columns

### Examine transformation effect
par( mfrow = c( 1, 2, 3 ) )                           # Format output plots to show side-by-side?
dds <- estimateSizeFactors(deHT)                   # ?
plot(countMatrix[,1:2], pch=16, cex=0.3)           # Plot the rlog transformed values
plot(log2(countMatrix[,1:2]), pch=16, cex=0.3)     # Plot the log transformed values
plot(countMatrix.rlog[,1:2], pch=16, cex=0.3)           # Plot the rlog transformed values

### Testing various conditions
test <- countMatrix[rowSums(countMatrix) > 17, ]
test <- t(test)
log.count <- log(count.trans)
log.test <- log(test)
rlog.test <- t(test)
#ir.species <- iris[, 5]

### Perform PCA with prcomp
count.pca <- prcomp(t(assay(deHTtrans)), center = TRUE, scale. = TRUE) # scale. = TRUE is highly advisable, but default is FALSE.
count.pca <- prcomp(count.trans, center = TRUE, scale. = TRUE) # scale. = TRUE is highly advisable, but default is FALSE.


### print method
#print(count.pca)

### plot method
plot(count.pca, type = "l")

### summary method
summary(count.pca)

### Predict PCs
#??? predict(ir.pca,  newdata=tail(log.ir, 2))

### Plot using biplot
biplot(count.pca,var.axes=FALSE,ylabs=NULL)
#biplot(data[1:2],var.axes=FALSE,ylabs=NULL)
biplot(ir.pca)

### Plot using ggbiplot
ggbiplot(count.pca, obs.scale = 1, var.scale = 1)

count.pca.df <- data.frame(count.pca$x)[1:2]
ggplot(count.pca.df, aes(PC1, PC2)) + geom_point(size=3) +
  geom_text(aes(label=row.names(count.pca.df)),hjust=0,vjust=0)
  

### Plot 3D ################

test <- countMatrix[1:20,1:17]
test.trans <- t(test)

pc <- princomp(iris[,1:4], cor=TRUE, scores=TRUE)
#pc2 <- prcomp(iris[,1:4], center = TRUE, scale. = TRUE) # scale. = TRUE is highly advisable, but default is FALSE.
pc2 <- prcomp(t(assay(deHTtrans)), center = TRUE, scale. = TRUE) # scale. = TRUE is highly advisable, but default is FALSE.
pc2 <- prcomp(count.trans, center = TRUE, scale. = TRUE) # scale. = TRUE is highly advisable, but default is FALSE.

pc2 <- count.pca

#pc <- princomp(test.trans, cor=TRUE, scores=TRUE)

### Calculate vectors for PC
coords <- NULL
for (i in 1:nrow(pc$loadings)) {
  coords <- rbind(coords, rbind(c(0,0,0),pc$loadings[i,1:3]))
}

coords <- NULL
for (i in 1:nrow(pc2$rotation)) {
  coords <- rbind(coords, rbind(c(0,0,0),pc$rotation[i,1:3]))
}

# Plot points, labels, component labels, and component vectors
plot3d(pc$scores[,1:3])
text3d(pc$scores[,1:3],texts=rownames(test.trans))
text3d(pc$loadings[,1:3], texts=rownames(pc$loadings), col="red")
lines3d(coords, col="red", lwd=4)

plot3d(pc2$x[,1:3])
text3d(pc2$x[,1:3],texts=rownames(count.trans))

###################################################################




### Plotting
plotMA(res, main="DESeq2", ylim=c(-2,2))
plotCounts(deHT, gene=which.min(res$padj), intgroup="Disease")




### Format for output
resultsDF <- as.data.frame(resOrdered)

### Fold change cutoffs
logFoldChange <- resultsDF$log2FoldChange
absLFC <- abs(logFoldChange)
lowVals <- absLFC<0.33
lowPvals <- resultsDF$padj>0.1
logFoldChange[lowVals] <- 0
logFoldChange[is.na(logFoldChange)] <- 0
logFoldChange[lowPvals] <- 0
resultsDF$log2FoldChange <- logFoldChange

# Reorder orginal counts
#reorderedHT <- iHTcounts[match(row.names(resultsDF),row.names(iHTcounts)),]
#reorderedOBS <- obsHTdata2[match(row.names(resultsDF),row.names(obsHTdata2)),]
#reorderedCTR <- ctrHTdata2[match(row.names(resultsDF),row.names(ctrHTdata2)),]

#head(reorderedOBS)
#head(reorderedCTR)
head(resultsDF)
#head(reorderedHT)
resultsDF <- addGene(resultsDF)
resultsDF <- addDescription(resultsDF)
names(resultsDF) <- c("baseMean","log2FoldChange", "lfcSE","stat","pvalue","padj","Gene","Descrip")

#outputDataFrame <- data.frame(row.names(resultsDF),reorderedOBS$median,reorderedCTR$median,reorderedOBS$sd,reorderedCTR$sd,resultsDF$log2FoldChange,resultsDF$padj)
outputDataFrame <- data.frame(row.names(resultsDF),resultsDF$Gene,resultsDF$log2FoldChange,resultsDF$padj,resultsDF$Descrip)
names(outputDataFrame) <- c("Ensembl-ID","Gene","Log2FoldChange","p-adj","Description")


write.table(as.data.frame(outputDataFrame[1:8000,]), file="z://Data/RNAseq HT neurons and tissue/2nd rerun/DE_genes_ANDREW/Obs_iHT_vs_Ctr_iHT--8000.txt",sep="\t",row.names = FALSE)


