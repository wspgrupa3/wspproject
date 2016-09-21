pca_MS <- function(expset,nrKl){
  
  Dane=exprs(expset)
  
  # X=Dane
  X=Dane[1:100,1:100] # Obliczenia d³ugo trwaj¹ dla ca³oœci, wiêc biorê tylko czêœæ danych aby pokazaæ jak to dzia³a. 
  
  # Wprowadzenie danych
  # Sprawdzimy wp³yw skalowania danych na wyniki analizy
  X2=t(t(X)/apply(X,2,sd))
  X3=log(X)
    
  # Analiza Sk³adowych G³ównych (Principal Component Analysis): ortogonalizujemy
  # kolumny X (cechy) i ustawiamy je w porzadku malejšcych wariancji.
  pca=prcomp(X); pca2=prcomp(X2); pca3=prcomp(X3)
  
  # Skalowanie wielowymiarowe (Multidimensional Scaling): rzutujemy obserwacje
  # na p³aszczyznê (podprzestrzeñ) najlepiej zachowujšcš odleg³o??ci miêdzy
  # obserwacjami. Pomijamy przypadek minimalizacji odleg³o??ci euklidesowej
  # (classical multidimensional scaling = metrical scaling  = principal coordinate
  # analysis), poniewa¿ jest równowa¿ny PCA: cmd=cmdscale(dist(X)).
  # Skalowanie Sammona minimalizuje wzglêdne odleg³o??ci miêdzy danymi.
  sa=sammon(dist(X),trace=FALSE)$points
  sa2=sammon(dist(X2),trace=FALSE)$points
  sa3=sammon(dist(X3),trace=FALSE)$points
  
  # Wykonanie 6 obrazków na 1 okienku wymaga zmiany parametrów rysowania.
  oldpar=par(no.readonly=TRUE)
  
  par(las=1, mfcol=c(3,2), oma=c(0,0,1,0))
  klasy=kmeans(X, nrKl)$cluster
  plot(pca$x[,1], pca$x[,2], main="original data", col=klasy)
  plot(pca2$x[,1],pca2$x[,2],main="scaled data",   col=klasy)
  plot(pca3$x[,1],pca3$x[,2],main="log data",      col=klasy)
  plot(sa[,1], sa[,2], main="original data", col=klasy)
  plot(sa2[,1],sa2[,2],main="scaled data",   col=klasy)
  plot(sa3[,1],sa3[,2],main="log data",      col=klasy)
  title("True classes",outer=TRUE)
  
  # Podzia³ dany przez klasteryzacjê metodš k-??rednich
  klasy= kmeans(X, nrKl)$cluster
  klasy2=kmeans(X2,nrKl)$cluster
  klasy3=kmeans(X3,nrKl)$cluster;
  x11(); par(las=1, mfcol=c(3,2), oma=c(0,0,1,0))
  plot(pca$x[,1], pca$x[,2], main="original data", col=klasy)
  plot(pca2$x[,1],pca2$x[,2],main="scaled data",   col=klasy2)
  plot(pca3$x[,1],pca3$x[,2],main="log data",      col=klasy3)
  plot(sa[,1], sa[,2], main="original data", col=klasy)
  plot(sa2[,1],sa2[,2],main="scaled data",   col=klasy2)
  plot(sa3[,1],sa3[,2],main="log data",      col=klasy3)
  title("K-means",outer=TRUE)
  
  # Podzia³ dany przez aglomeracyjnš klasteryzacjê hierarchicznš
  hcl= hclust(dist(X));   klasy = as.vector(cutree(hcl,nrKl))
  hcl2=hclust(dist(X2)); klasy2 = as.vector(cutree(hcl2,nrKl))
  hcl3=hclust(dist(X3)); klasy3 = as.vector(cutree(hcl3,nrKl))
  x11(); par(las=1, mfcol=c(3,2), oma=c(0,0,1,0))
  plot(pca$x[,1], pca$x[,2], main="original data", col=klasy)
  plot(pca2$x[,1],pca2$x[,2],main="scaled data",   col=klasy2)
  plot(pca3$x[,1],pca3$x[,2],main="log data",      col=klasy3)
  plot(sa[,1], sa[,2], main="original data", col=klasy)
  plot(sa2[,1],sa2[,2],main="scaled data",   col=klasy2)
  plot(sa3[,1],sa3[,2],main="log data",      col=klasy3)
  title("Hierarch clust",outer=TRUE)
  
  par(oldpar)
  
  x11(); plot(hcl) # dendrogram

}