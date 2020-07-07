# Proteomics

```
---
title: "FOUNDIN_proteomics"
author: "Mark Cookson"
date: "4/11/2019"
output: html_document
notes: Analysis of TMT proteomics of 10 PPMI cell lines run across three replicates. Data supplied was as ratio for each protein relative to a pool of all samples combined. 
---
```

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
#set libraries and directory
library("impute")
library("multtest")
library("ggplot2")
library("grid")
library("pheatmap")
library("biomaRt")
previous_theme <- theme_set(theme_bw())
Dir <- "~/PATH/TO/WORK/DIR/"
```

```{r}
#read in data, label rows with accession number then tidy up to get a matrix of expression ratios and a colData object for anotation
TMT0 <- read.table(file=paste(Dir,"DAneuron_proteomics-PPMI-ID.txt",sep=""),header=T, stringsAsFactors=F)
row.names(TMT0) <- TMT0$Accession
TMT0 <- TMT0[,-1]
sample <- paste0("PPMI",as.vector(t(as.data.frame(strsplit(as.vector(colnames(TMT0)),"X")))[,2]))
colnames(TMT0) <- sample
colData <- as.data.frame(sample)
colData$PPMIID <- as.vector(t(as.data.frame(strsplit(as.vector(sample),"[.]")))[,1])
colData$replicate <- as.vector(t(as.data.frame(strsplit(as.vector(sample),"[.]")))[,2])
row.names(colData) <- colData$sample
```


```{r}
#remove proteins with >50% missingness
TMT <- TMT0[rowSums(is.na(TMT0)) <ncol(TMT0)/2,]
dim(TMT)
any(is.na(TMT))


#look at distribution of values, note some at 100 or 0.01
distribution <- stack(as.data.frame(TMT))
a <- ggplot(data=distribution, aes(x=ind, y=log10(values)))
a <- a + geom_boxplot()
a <- a + ggtitle("TMT ratio distribution raw")
a

#simple filter to remove proteins with extreme values
id <- which(apply(log10(TMT), 1, function (x) all(abs(x) <= 1)))
TMT <- TMT[id, ]

#simple scale data in columns
center_colmeans <- function(x) {
    xcenter = colMeans(x)
    x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

unscaled <- stack(as.data.frame(TMT))
a <- ggplot(data=unscaled, aes(x=ind, y=values))
a <- a + geom_boxplot()
a <- a + ggtitle("TMT ratio distribution filtered unscaled")
a

TMT <- center_colmeans(TMT)
scaled <- stack(as.data.frame(TMT))
b <- ggplot(data=scaled, aes(x=ind, y=values))
b <- b + geom_boxplot()
b <- b + ggtitle("TMT ratio distribution filtered scaled")
b

#impute missing values
temp1 <- as.matrix(as.data.frame(lapply(TMT, as.numeric)))
row.names(temp1) <- row.names(TMT)
temp1 <- as.matrix(impute.knn(temp1)$data)
TMT <- temp1
dim(TMT)
any(is.na(TMT))

```
```{r}
#PCA and plot
PCs <- prcomp(TMT)
eigs <- PCs$sdev^2
pctvarexp <- round(eigs/sum(eigs),3)*100

data2plot <- PCs$rotation
data2plot <- cbind(data2plot, colData)

a <- ggplot(data=data2plot, aes(x=PC1, y=PC2, color=PPMIID, shape=replicate))
a <- a + geom_point(size=3)
a <- a + xlab(paste("PC1 (",pctvarexp[1],"% variance)"))
a <- a + ylab(paste("PC2 (",pctvarexp[2],"% variance)"))
a <- a + ggtitle("PCs without regression for batch")
a

#regress out run effects
temp <- matrix(0, nrow=nrow(TMT), ncol=ncol(TMT))
row.names(temp) <- row.names(TMT)
colnames(temp) <- colnames(TMT)

for (i in 1:length(row.names(TMT)))
{
  thisProtein <- as.vector(row.names(TMT)[i])
  data2test <- colData
  data2test$value <- as.vector(TMT[thisProtein,])
  data2test$valAdj <- glm(value ~ replicate, data = data2test)$resid
  data2test$adjZ <- (data2test$valAdj - mean(data2test$valAdj))/sd(data2test$valAdj)
  temp[i,] <- data2test$adjZ
}

TMT<- temp

```

```{r}
#final cleaned dataset heatmap
pheatmap(mat=TMT, annotation_col = colData[,2:3], show_rownames = F)

#PCA and plot
PCs <- prcomp(TMT)
eigs <- PCs$sdev^2
pctvarexp <- round(eigs/sum(eigs),3)*100

data2plot <- PCs$rotation
data2plot <- cbind(data2plot, colData)

a <- ggplot(data=data2plot, aes(x=PC1, y=PC2, color=PPMIID, shape=replicate))
a <- a + geom_point(size=3)
a <- a + xlab(paste("PC1 (",pctvarexp[1],"% variance)"))
a <- a + ylab(paste("PC2 (",pctvarexp[2],"% variance)"))
a <- a + ggtitle("PCs after regression for batch")
a
#looks like line-to-line variability drives much of PC1, PC2 may capture technical replicates
```

```{r}
#identify variable proteins by ANOVA
results <- data.frame("protein"=NA,"F"=NA,"p"=NA)

for (i in 1:length(row.names(TMT)))
{
  thisProtein <- as.vector(row.names(TMT)[i])
  data2test <- colData
  data2test$value <- as.vector(TMT[thisProtein,])
  model <- aov(value~PPMIID, data=data2test)
  results[i,] <- c(thisProtein, summary(model)[[1]][["F value"]][1], summary(model)[[1]][["Pr(>F)"]][1])
}
results$F <- as.numeric(results$F)
results$p <- as.numeric(results$p)
results <- results[order(results$p),]
results$adjp <- mt.rawp2adjp(results$p, c("Bonferroni"))$adjp
head(results)

#bring in annotations
annotations <- read.table(file=paste0(Dir, "annotations.txt"), header=T, sep="\t")
results <- cbind(results, annotations[match(results$protein, annotations$Accession),])
results[,5:10][results[,5:10]==""] <- NA
write.table(results, file=paste0(Dir, "FOUNDIN_TMT_preliminary_ANOVA.txt"), quote=F, row.names = F, sep="\t")
mat <-TMT[match(results[!is.na(results$Gene.Symbol),]$Accession[1:20], row.names(TMT)),]
row.names(mat) <- as.vector(results$Gene.Symbol[match(row.names(mat), results$Accession)])
  
pheatmap(mat=mat, annotation_col = colData[,2:3], show_rownames = T)
```


