### opis danych po normalizacji metodami graficznymi

#Marta Zarebska

graficzne<- function(data, folder_path)
    
  {

library(affy)
library(estrogen) 
#odczytanie warto?ci ekspresji
data=exprs(data_eSet)

#metody graficzne dla znormalizowanych danych z eSEt

#wykres pudelkowy
png(file.path(folder_path, 'wykres_pudelkowy.png'));
boxplot(data, main="Wykres pude?kowy dla danych znormalizowanych");
dev.off()

#histogram
png(file.path(folder_path, 'histogram.png'));
plotDensity(data, main="Histogram dla danych znormalizowanc", xlab="Sygna? sondy", ylab="G?sto??");
dev.off()

#wykres degradacj
degradation=AffyRNAdeg(data)
di = get('plot_device_1', envir = .GlobalEnv);
dev.set(di);
plotAffyRNAdeg(degradation);
dev.copy(png, file.path(folder_path, 'degradacja.png'));
dev.off();

#heatmap
rs=apply(exprs(data), 1, sd)
sel_=order(rs, decreasing = TRUE)[1:50]
png(filename="D:/dysk D/studia/biotechnologia/4rok/VIII semestr/wsp_laboratorium/lab_2") #sciezka docelowa
di = get('plot_device_2', envir = .GlobalEnv);
dev.set(di);
heatmap(exprs(data)[sel_, ], col = gentlecol(256))
dev.copy(png, file.path(folder_path, 'heatmap.png'));
dev.off();

}
