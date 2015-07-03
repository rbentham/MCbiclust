SilhoetteHclust <- function(gem,seed,max.groups = 20){
require(cluster)
cor.mat <- cor(t(gem[,seed]))

cor.hclust <- hclust(dist(cor.mat))
cor.hclust.dist <- dist(CCLE_mito_cor)

sil.value <- seq(length = max.groups-1)
for(i in 2:10){
  si2 <- silhouette(x = cutree(cor.hclust,k = i),dist = cor.hclust.dist)
  sil.value[i-1] <- mean(si2[,3])}

plot(seq(length = max.groups-1)+1, sil.value, xlab="k",ylab="Mean silhoette width")

k1 <- which.max(sil.value) + 1

si2 <- silhouette(x = cutree(cor.hclust,k = k1), dist = cor.hclust.dist)
row.names(si2) <- row.names(gem)

return(si2)
}