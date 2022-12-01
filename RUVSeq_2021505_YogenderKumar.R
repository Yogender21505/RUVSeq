# load the libraries
library(RUVSeq)
library(zebrafishRNASeq)

# load the data
data(zfGenes)

# glimpse of the data
head(zfGenes)
tail(zfGenes)

# filter out the genes with little to no expression
filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
head(filter)

filtered <- zfGenes[filter,]

genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

head(genes)

# convert to S4 object:
x <- as.factor(rep(c("Ctl", "Trt"), each=3))
set <- newSeqExpressionSet(as.matrix(filtered),
                           phenoData = data.frame(x, row.names=colnames(filtered)))
set


# Check for variability in the replicates:
library(RColorBrewer)

colors <- brewer.pal(3, "Set2")

plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x],main="variability_RLE")
plotPCA(set, col=colors[x], cex=1.2,main="variability_PCA")


# Upper quantile normalization to reduce batch effects:
set <- betweenLaneNormalization(set, which="upper")

plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x],main="Normalized_RLE")
plotPCA(set, col=colors[x], cex=1.2,main="Normalized_PCA")

# Spike ins to detect unwanted variations:
set1 <- RUVg(set, spikes, k=1)
pData(set1)

plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x],main="spike_RLE_plot")
plotPCA(set1, col=colors[x], cex=1.2,main="spike_PCA_plot")

# Differential Gene Expression Analysis ####
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts(set1),
                              colData = pData(set1),
                              design = ~ W_1 + x)

dds <- DESeq(dds)
res <- results(dds)
res

# convert to dataframe
output <- as.data.frame(res)


# volcano plots:
final <- data.frame(row.names(output))
final[,c(2,3)] <- output[,c(2,5)]


names(final)[1] <- "Gene"

#Top50 upregulation genes pvalue<0.5
upreg <- final[final$pvalue < 0.05 & final$log2FoldChange>0,]
names(upreg)[1] <- "Upregulatedgenes"
data_new1<-upreg
rownames(data_new1) <- 1:nrow(data_new1)
Top50upregulatedgenes <- data_new1[1:50,]



#Top50 downregulation genes pvalue<0.5
downreg <- final[final$pvalue < 0.05 & final$log2FoldChange<0,]
names(downreg)[1] <- "Downregulatedgenes"
data_new2<-downreg
rownames(data_new2) <- 1:nrow(data_new2)
Top50Downregulatedgenes <- data_new2[1:50,]

#Exporting data into CSV
write.csv(Top50Downregulatedgenes,"C:\\Users\\kumar\\OneDrive\\Documents\\Semester\\SEM3\\Books\\GMB\\2021505_YogenderKumar\\2021505_YogenderKumar\\Downregulated_Top50_2021505_Yogender.csv", row.names = FALSE)
write.csv(Top50upregulatedgenes,"C:\\Users\\kumar\\OneDrive\\Documents\\Semester\\SEM3\\Books\\GMB\\2021505_YogenderKumar\\2021505_YogenderKumar\\Upregulated_Top50_2021505_Yogender.csv", row.names = FALSE)
# Make a basic volcano plot
with(final, plot(log2FoldChange, -log10(pvalue), pch = 20, main = "2021505_Yogender_DEGs"))

with(subset(final, pvalue < .05 & abs(log2FoldChange) > 1), points(log2FoldChange, -log10(pvalue), pch = 20, col = "blue"))