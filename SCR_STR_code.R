library("edgeR")
library("ggplot2")
library("limma")
library(xlsx)
setwd("D:\\New folder")
rawCountTable <-read.delim("countfile.txt",  header=T, sep="\t", row.names=1)
sampleInfo <- read.table("design.csv", header=TRUE, sep=",", row.names=1)
head(rawCountTable)
nrow(rawCountTable)
dgeFull <- DGEList(rawCountTable, group=sampleInfo$condition)
dgeFull

pseudoCounts <- log2(dgeFull$counts+1)
head(pseudoCounts)
hist(pseudoCounts[,"HISAT2_STR"],col = "orange",xlab = "STR")
hist(pseudoCounts[,"HISAT2_SCR"],col = "pink",xlab = "SCR")
boxplot(pseudoCounts, main = "Boxplot", border = "brown",names = c("HISAT_STR","HISAT2_SCR"),col = c("lightgreen","lightblue"),notch = TRUE,varwidth = TRUE,)
sampleDists <- as.matrix(dist(t(pseudoCounts)))
sampleDists
plot(sampleDists)

dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)
head(dgeFull$counts)
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples
head(dgeFull$counts)
eff.lib.size <- dgeFull$samples$lib.size*dgeFull$samples$norm.factors
normCounts <- cpm(dgeFull)
pseudoNormCounts <- log2(normCounts + 1)
boxplot(pseudoNormCounts, col="gray", las=3)
bcv <- 0.1
nRows <- 32230
nCols <- 2
counts <- matrix( rnbinom(nRows * nCols,size=1/bcv^2,mu=10),nrow = nRows, ncol = nCols)
counts <- rawCountTable
dgeFull <- DGEList(counts=counts, group=1:2)
dgeFull$common.dispersion=0.1
et <- exactTest(dgeFull)
dgeFull
et
group <- factor(c("STR","SCR"))
design <- model.matrix(~group)
design
fit <- glmFit(dgeFull, design ,robust=TRUE)
lrt <- glmLRT(fit)
lrt
topTags(lrt)
colnames(design)
o <- order(lrt$table$PValue)
cpm(dgeFull)[o[1:10],]
summary(decideTests(lrt))
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
table3 <- (topTags(lrt, n = 32228)[["table"]])
table3
write.xlsx(table3, "table3.xlsx")
