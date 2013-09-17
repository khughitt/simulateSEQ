
library(preprocessCore)
library(limma)
library(sva)
library(tukeyGH)
library(matrixStats)

#-------------- Define Trends --------------------------------#
# trend fcts
trend1 <- function(x, k, a) (1 / (1 + exp(-a*(x - k))) - .5)*5
trend2 <- function(x, k, a)  (exp(-((x - k)/a)^2) - .5)*5

# time points
x <- c(0, 2, 4, 6, 8, 10, 12, 14, 16, 18)
# mean trend of group1
mean.grp1 <-  trend1(x, k=5, a=1)
plot(x, mean.grp1, col="gray40", pch=16, xlab="Time (Hrs)",
     ylab="Standardized Expression Level (SEL)", xaxt="n",
     yaxt="n", ylim=c(-3,3), main="3 main trends")
axis(1, x, labels=x)
axis(2, -3:3, labels=-3:3)
# mean trend of group2
mean.grp2 <- trend1(x, k=12, a=-2)
points(x, mean.grp2, col="gray70", pch=16)
# mean trend of group3
mean.grp3 <- trend2(x, k=15, a=5)
points(x, mean.grp3, col="red")
# plot curves
curve(trend1(x, k=5, a=1), col="gray40", add=TRUE, lty=2)
curve(trend1(x, k=12, a=-2), col="gray70", add=TRUE, lty=2)
curve(trend2(x, k=15, a=5), col="red", add=TRUE, lty=2)


#-------------------- Data prep  -----------------------#

library(GEOquery)
gse <- getGEO("GSE11300")

gse[[1]]

expr <- log(exprs(gse[[1]]))

# quantile normalize
rn <- rownames(expr)
cn <- colnames(expr)
expr <- normalize.quantiles(as.matrix(expr))
rownames(expr) <- rn
colnames(expr) <- cn
head(expr)

# compute rowmeans and rowsds
sds <- rowSds(expr)
ms <- rowMeans(expr)

# plot sd by ave
plot(ms, sds, xlab="ave", ylab="sd", 
     col=rgb(77, 77, 77, 50, maxColorValue=255), pch=19,
     cex=0.8, main=paste(nrow(expr), "genes"))

grp1 <- 1:2000
grp2 <- 2001:5000
grp3 <- 5001:7000
ns <- length(x) # number of bio conds
G <- nrow(expr) # number of genes

stdExpr <- (expr - ms)/sds

rSEL <- replicate(ns, 1.5*sample(stdExpr, G, replace=TRUE))
dim(rSEL)
colnames(rSEL) <- x


# set first 3000 genes to follow trend1
rSEL[grp1,] <- t(matrix(rep(mean.grp1, length(grp1)), ns))
# set next 3000 genes to follow trend2
rSEL[grp2,] <- t(matrix(rep(mean.grp2, length(grp2)), ns))
# set next 3000 genes to follow trend3
rSEL[grp3,] <- t(matrix(rep(mean.grp3, length(grp3)), ns))

matplot(t(rSEL[-c(grp1, grp2, grp3),]), type="l", lty=3,
        col=rgb(0, 0, 0, 5, maxColorValue=255), 
        xaxt="n", yaxt="n", ylab="SEL", xlab="Time (Hrs)",
        ylim=c(-3,3))
axis(1, 1:ncol(rSEL), labels=x)
axis(2, -3:3, labels=-3:3)
matplot(t(rSEL[grp1,]), type="l", lty=3,
        col="red", add=T)
matplot(t(rSEL[grp2,]), type="l", lty=3,
        col="red", add=T)
matplot(t(rSEL[grp3,]), type="l", lty=3,
        col="red", add=T)

natureGeneProfile <- ms
hist(natureGeneProfile, freq=F, breaks=40)


source("/Users/kokrah/Dropbox/library_size/code/setupFunctions.R")

natureGeneProfileMat <- rSEL*sds + natureGeneProfile
hist(natureGeneProfileMat, breaks=40)

genData <- function(x){
    res <-  replicate(2, generateArray(natureGeneProfile=x, 
                       dispersion=sample(sds, G, replace=T)))
    as.data.frame(res)
}

simExpr <- apply(natureGeneProfileMat, 2, genData) 
simExpr <- as.data.frame(simExpr)
rownames(simExpr) <- paste("gene_", 1:G, sep="")
head(simExpr)
dim(simExpr)


boxplot(simExpr, las=3)
hist(as.matrix(simExpr), freq=F, breaks=40)

stdmat <- (simExpr - rowMeans(simExpr)) / rowSds(simExpr)
hld <- list()
stdmat <- for(i in 1:(ncol(simExpr)/2)){
  k <- i*2 - 1
  hld[[i]] <- rowMeans(stdmat[,k:(k+1)])
}  
stdmat <- as.data.frame(hld)
colnames(stdmat) <- x

matplot(t(stdmat[-c(grp1, grp2, grp3),]), type="l", lty=3,
        col=rgb(0, 0, 0, 5, maxColorValue=255), 
        xaxt="n", yaxt="n", ylab="SEL", xlab="Time (Hrs)",
        ylim=c(-3,3))
axis(1, 1:ncol(stdmat), labels=x)
axis(2, -3:3, labels=-3:3)

matplot(t(stdmat[grp1,]), type="l", lty=3,
        col=rgb(0, 0, 0, 50, maxColorValue=255), 
        xaxt="n", yaxt="n", ylab="SEL", xlab="Time (Hrs)",
        ylim=c(-3,3))
axis(1, 1:ncol(stdmat), labels=x)
axis(2, -3:3, labels=-3:3)

matplot(t(stdmat[grp2,]), type="l", lty=3,
        col=rgb(0, 0, 0, 50, maxColorValue=255), 
        xaxt="n", yaxt="n", ylab="SEL", xlab="Time (Hrs)",
        ylim=c(-3,3))
axis(1, 1:ncol(stdmat), labels=x)
axis(2, -3:3, labels=-3:3)

matplot(t(stdmat[grp3,]), type="l", lty=3,
        col=rgb(0, 0, 0, 50, maxColorValue=255), 
        xaxt="n", yaxt="n", ylab="SEL", xlab="Time (Hrs)",
        ylim=c(-3,3))
axis(1, 1:ncol(stdmat), labels=x)
axis(2, -3:3, labels=-3:3)

#--------------------------------------------------------------------------#



res <- qFit(simExpr)


head(res)

hist(res$g)

ms <- rowMeans(simExpr)
sds <- rowSds(simExpr)

plot(res$g, sds, cex=0.3)

sum(sel1 <- res$g > .8)

y <- stdmat[sel1,]

hist(as.matrix(y), breaks=40)

matplot(t(y), type="l", lty=3,
        col=rgb(0, 0, 0, 50, maxColorValue=255), 
        xaxt="n", yaxt="n", ylab="SEL", xlab="Time (Hrs)",
        ylim=c(-3,3))
axis(1, 1:ncol(stdmat), labels=x)
axis(2, -3:3, labels=-3:3)



sum(sel2 <- res$g < -.8)

y <- stdmat[sel2,]

hist(as.matrix(stdmat[sel2,]), breaks=40)

matplot(t(stdmat[sel2,]), type="l", lty=3,
        col=rgb(0, 0, 0, 50, maxColorValue=255), 
        xaxt="n", yaxt="n", ylab="SEL", xlab="Time (Hrs)",
        ylim=c(-3,3))
axis(1, 1:ncol(stdmat), labels=x)
axis(2, -3:3, labels=-3:3)

matplot(t(stdmat[grp1,]), type="l", lty=3, add=T,
        col=rgb(100, 0, 0, 50, maxColorValue=255))

matplot(t(stdmat[grp2,]), type="l", lty=3, add=T,
        col=rgb(100, 100, 0, 50, maxColorValue=255))

matplot(t(stdmat[grp3,]), type="l", lty=3, add=T,
        col=rgb(100, 100, 150, 50, maxColorValue=255))

sum(sel2 <- abs(res$g) < .01)

y <- stdmat[sel2,]

hist(as.matrix(stdmat[sel2,]), breaks=40)

matplot(t(stdmat[sel2,]), type="l", lty=3,
        col=rgb(0, 0, 0, 50, maxColorValue=255), 
        xaxt="n", yaxt="n", ylab="SEL", xlab="Time (Hrs)",
        ylim=c(-3,3))
axis(1, 1:ncol(stdmat), labels=x)
axis(2, -3:3, labels=-3:3)

matplot(t(stdmat[grp1,]), type="l", lty=3, add=T,
        col=rgb(100, 0, 0, 50, maxColorValue=255))

matplot(t(stdmat[grp2,]), type="l", lty=3, add=T,
        col=rgb(100, 100, 0, 50, maxColorValue=255))

matplot(t(stdmat[grp3,]), type="l", lty=3, add=T,
        col=rgb(100, 100, 150, 50, maxColorValue=255))