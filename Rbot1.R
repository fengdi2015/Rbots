########## fun heatmap ###################
x  <- as.matrix(mtcars)
plot(x)
data=t(scale(t(x)))
hr <- hclust(as.dist(1-cor(t(data),method="pearson")),method="complete");
hc <- hclust(as.dist(1-cor(data,method="spearman")), method="complete");
##cut tree
rcl=cutree(hr,h=max(hr$height/2));rowcl=sample(rainbow(16))[as.vector(rcl)];
ccl=cutree(hc,h=max(hc$height/2));colcl=sample(rainbow(16))[as.vector(ccl)];
colorpanel <- brewer.pal(11,"RdBu");
hmcols<-colorRampPalette(c("red","white","blue"))(256);
#heat.colors(10,alpha=1)
heatmap(data,Rowv = as.dendrogram(hr),Colv = as.dendrogram(hc),
        scale="row",cexCol = 1,cexRow = 1,RowSideColors = rowcl,ColSideColors = colcl,
        col=hmcols);

## for ggplot version
library(ggplot2);
library(reshape2);
data2 <- melt(data);
p <- ggplot(data2, aes(Var2, Var1))+
     geom_tile(aes(fill = as.numeric(value), colour = "white"))+
     scale_fill_gradient(low = "white",high = "steelblue")
p

