# dist, hclust
x <- matrix(rnorm(100),nrow=5)
dim(x)

args(dist)

dist(x)
d <- dist(x)

as.matrix(d)

hc <- hclust(d, method = 'single')
plot(hc)



# heatmap

mtcars
heatmap(scale(as.matrix(mtcars)))

library(gplots)

heatmap.2(scale(as.matrix(mtcars)))

library(pheatmap)
pheatmap(scale(as.matrix(mtcars)))

# ALL data
library(ALL)
data(ALL)

ALL

pData(ALL)[1:5,]

expr <- exprs(ALL)
head(expr)

# mad 유전자 발현량이 차이가 있는 값을 찾기 위해 사용한다.

mad_ <- apply(expr,1,mad)
quantile(mad_)

f.expr <- expr[mad_>0.5,]
dim(f.expr)

## dist, hclust

head(f.expr)
head(t(f.expr))[,1:5]

d <- dist(t(f.expr))
d
head(as.matrix(d))[,1:5]

hc <- hclust(d)
plot(hc)

hc$order
hc$labels[hc$order]


table(cutree(hc, k=3))

ALL$hc.group = cutree(hc, k=3)
head(pData(ALL))

## heatmap

hist(f.expr)


gene.mean <- apply(f.expr,1,mean)
f.cen.expr = f.expr - gene.mean

hist(f.cen.expr)

heatmap.2(f.cen.expr)

heatmap.2(f.cen.expr,trace='none')

heatmap.2(f.cen.expr,trace='none',scale='none')

# color 
breaks <- seq(-2,2,0.01)
library(RColorBrewer)
col <- colorRampPalette(brewer.pal(10,'RdBu'))(length(breaks)-1)
col

heatmap.2(f.cen.expr,trace='none',scale='none',breaks = breaks,col=rev(col))

h2 <- heatmap.2(f.cen.expr,trace='none',scale='none',breaks = breaks,col=rev(col))

ALL$h2.group <- cutree(as.hclust(h2$colDendrogram),k=3)

table(ALL$hc.group)
table(ALL$h2.group)

#colorbar
ifelse(grepl('B',ALL$BT),'B','T')
col.bt <- ifelse(grepl('B',ALL$BT),'blue','green')

heatmap.2(f.cen.expr,trace='none',scale='none',breaks = breaks,col=rev(col)
          ,ColSideColors = col.bt)

## consensus clustering
library(ConsensusClusterPlus)
cc.res <- ConsensusClusterPlus(f.cen.expr,maxK = 5,plot='png',title='ALL')

ALL$cc.group <- cc.res[[2]]$consensusClass

## nmf clustering
library(NMF)

hist(f.expr) # only >0 value

nmf.res <- nmf(f.expr,rank = 2:5)
nmf.res$fit
consensusmap(nmf.res$fit)
args(write)
save(x = nmf.res, file='ALL/ALL_nmf.rda')

plot(nmf.res)

ALL$nmf.group <- as.numeric(predict(nmf.res$fit[[1]]))


# extract feature
c1 <- rownames(f.expr)[extractFeatures(nmf.res$fit[[1]])[[1]]]
c2 <- rownames(f.expr)[extractFeatures(nmf.res$fit[[1]])[[2]]]

f.cen.expr[c(c1, c2),order(ALL$nmf.group)]

heatmap.2(f.cen.expr[c(c1, c2),order(ALL$nmf.group)], scale = 'none', breaks=breaks, col=col, trace = 'none', Colv=F, Rowv=F, dendrogram = "none")

### compare clusters

class.col=brewer.pal(n=3, name="Set1")

hc.col=class.col[ALL$h2.group]
cc.col=class.col[ALL$cc.group]
nmf.col=class.col[ALL$nmf.group]

par(mfrow=c(3,1), mar=c(0,0,0,0))
barplot(rep(1, ncol(ALL)), col=hc.col)
barplot(rep(1, ncol(ALL)), col=cc.col)
barplot(rep(1, ncol(ALL)), col=nmf.col)
