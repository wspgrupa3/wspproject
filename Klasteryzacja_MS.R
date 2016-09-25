Klasteryzacja<- function(expset, nrKl, folder_path){
  
  
  # Magdalena Szeremet
  # Funkcja skaluje dane mikromacierzowe, przeprowadza PCA i klasteryzuje metod? k-means i metod? hierarchiczn?.
  
  ### SET THIS FOR DATA MAX CAP ###
  max_cols = 4;
  max_genes = 10;
  #################################
  library('Biobase')
  library('BiocGenerics')
  library('parallel')
  library('MASS')
  library('nlme')
  
  Dane=exprs(expset)
  
  data_size = dim(Dane);
  last_col = min(max_cols, data_size[2]);
  X=Dane[1:max_genes,1:last_col]; # Obliczenia d?ugo trwaj? dla ca?o?ci, wi?c bior? tylko cz??? danych aby pokaza? jak to dzia?a. 
  
  # Wprowadzenie danych
  # Sprawdzimy wp?yw skalowania danych na wyniki analizy
  X2=t(t(X)/apply(X,2,sd))
  X3=log(X)
    
  # Analiza Sk?adowych G??wnych (Principal Component Analysis): ortogonalizujemy
  # kolumny X (cechy) i ustawiamy je w porzadku malej?cych wariancji.
  pca=prcomp(X); pca2=prcomp(X2); pca3=prcomp(X3)
  
  # Skalowanie wielowymiarowe (Multidimensional Scaling): rzutujemy obserwacje
  # na p?aszczyzn? (podprzestrze?) najlepiej zachowuj?c? odleg?o??ci mi?dzy
  # obserwacjami. Pomijamy przypadek minimalizacji odleg?o??ci euklidesowej
  # (classical multidimensional scaling = metrical scaling  = principal coordinate
  # analysis), poniewa? jest r?wnowa?ny PCA: cmd=cmdscale(dist(X)).
  # Skalowanie Sammona minimalizuje wzgl?dne odleg?o??ci mi?dzy danymi.
  sa=sammon(dist(X),trace=FALSE)$points
  sa2=sammon(dist(X2),trace=FALSE)$points
  sa3=sammon(dist(X3),trace=FALSE)$points
  
  # Wykonanie 6 obrazk?w na 1 okienku wymaga zmiany parametr?w rysowania.
  oldpar=par(no.readonly=TRUE)
  
  di = get('plot_device_1', envir = .GlobalEnv);
  dev.set(di);
  par(las=1, mfcol=c(3,2), oma=c(0,0,1,0))
  klasy=kmeans(X, nrKl)$cluster
  plot(pca$x[,1], pca$x[,2], main="original data", col=klasy)
  plot(pca2$x[,1],pca2$x[,2],main="scaled data",   col=klasy)
  plot(pca3$x[,1],pca3$x[,2],main="log data",      col=klasy)
  plot(sa[,1], sa[,2], main="original data", col=klasy)
  plot(sa2[,1],sa2[,2],main="scaled data",   col=klasy)
  plot(sa3[,1],sa3[,2],main="log data",      col=klasy)
  title("True classes",outer=TRUE)
  dev.copy(png, file.path(folder_path, 'True_classes.png'));
  dev.off();
  
  # Podzia? dany przez klasteryzacj? metod? k-?rednich
  klasy= kmeans(X, nrKl)$cluster
  klasy2=kmeans(X2,nrKl)$cluster
  klasy3=kmeans(X3,nrKl)$cluster;
  di = get('plot_device_2', envir = .GlobalEnv);
  dev.set(di);
  par(las=1, mfcol=c(3,2), oma=c(0,0,1,0))
  plot(pca$x[,1], pca$x[,2], main="original data", col=klasy)
  plot(pca2$x[,1],pca2$x[,2],main="scaled data",   col=klasy2)
  plot(pca3$x[,1],pca3$x[,2],main="log data",      col=klasy3)
  plot(sa[,1], sa[,2], main="original data", col=klasy)
  plot(sa2[,1],sa2[,2],main="scaled data",   col=klasy2)
  plot(sa3[,1],sa3[,2],main="log data",      col=klasy3)
  title("K-means",outer=TRUE)
  dev.copy(png, file.path(folder_path, 'k_means.png'));
  dev.off();
  
  # Podzia? dany przez aglomeracyjn? klasteryzacj? hierarchiczn?
  hcl= hclust(dist(X));   klasy = as.vector(cutree(hcl,nrKl))
  hcl2=hclust(dist(X2)); klasy2 = as.vector(cutree(hcl2,nrKl))
  hcl3=hclust(dist(X3)); klasy3 = as.vector(cutree(hcl3,nrKl))
  
  png(file.path(folder_path, 'Hierarch_clust.png'));
  par(las=1, mfcol=c(3,2), oma=c(0,0,1,0))
  plot(pca$x[,1], pca$x[,2], main="original data", col=klasy)
  plot(pca2$x[,1],pca2$x[,2],main="scaled data",   col=klasy2)
  plot(pca3$x[,1],pca3$x[,2],main="log data",      col=klasy3)
  plot(sa[,1], sa[,2], main="original data", col=klasy)
  plot(sa2[,1],sa2[,2],main="scaled data",   col=klasy2)
  plot(sa3[,1],sa3[,2],main="log data",      col=klasy3)
  title("Hierarch clust",outer=TRUE)
  dev.off();
  
  png(file.path(folder_path, 'hcl_dendrogram.png'));
  par(oldpar)
  plot(hcl) # dendrogram
  dev.off();

}