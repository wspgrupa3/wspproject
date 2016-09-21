klast <- function(dane)
{
library(MASS); 
library(cluster)

# hclust- hierarchiczne(average)
X=exprs(dane)
d=as.dist(1-cor(X[,],method='spearman')) 
h=hclust(d, method = "average", members=NULL)
h2=hclust(d, method = "median", members=NULL)
h3=hclust(d, method = "centroid", members=NULL)
h4=hclust(d, method='ward', members=NULL)

#Dendrogram
png('dendrogram-average.png')
plot(as.dendrogram(h),
     main="Dendrogram'hclust-average'", sub="expressionset ")
dev.off()

png('dendrogram-median.png')
plot(as.dendrogram(h2),
     main="Dendrogram'hclust-median'", sub="expressionset")
dev.off()

png('dendrogram_centroid.png')
plot(as.dendrogram(h3),
     main="Dendrogram'hclust-centroid'", sub="expressionset")
dev.off()

## Dendrogram Diana
di=diana(t(d))
png('Dendrogram_diana.png')
plot(di,main="Dendrogram 'diana'", sub="ExpressionSet ", ask=F)
dev.off()

a=numeric(10)
for (k in 2:10) {
  a[[k]]=pam(t(d), k)$silinfo$avg.width 
k.best=which.max(a)+1 
cat("Nalepsza liczba klastrów dla danych:", k.best, "\n")
}

# Klastrowanie Pam
p=pam(t(d), k.best)
png('Pam_naj-l-klastów.png')
plot(p,ask=T,main='Wyniki klastrowania PAM z najlepsz¹ liczb¹ klastrów')
dev.off()
png('PAM-5klastrów.png')
paalt=pam(t(d), 5)
pamp=plot(paalt,ask=T, main='Wyniki klastrowania PAM')
dev.off
result<-list(k.best)
return (result)
}