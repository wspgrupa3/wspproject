klas_hier<- function(data) #dane po normalizacji u¿yte do klasyfikacji hierarchicznej
{
#method: average, centroid, complete
#metoda 1
data_1=data
distance = dist(t(data_1))
image(as.matrix(distance))
hc_l = hclust(distance, method = "average")
#(chyba) mo¿na normalnie: png("average.png")
png(filename="D:/dysk D/studia/biotechnologia/4rok/VIII semestr/wsp_laboratorium/lab_2") #sciezka docelowa
plot(hc_l)
dev.off()

#metoda 2
distance = dist(t(data_1))
image(as.matrix(distance))
hc_l = hclust(distance, method = "centroid")
#png("centroid.png")
png(filename="D:/dysk D/studia/biotechnologia/4rok/VIII semestr/wsp_laboratorium/lab_2")
plot(hc_l)
dev.off()

#metoda 3
distance = dist(t(data_1))
image(as.matrix(distance))
hc_l = hclust(distance, method = "complete")
#png("complete.png")
png(filename="D:/dysk D/studia/biotechnologia/4rok/VIII semestr/wsp_laboratorium/lab_2")
plot(hc_1)
dev.off()

}
