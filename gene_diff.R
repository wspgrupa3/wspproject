### PORÓWNANIE DANYCH METODAMI STATYSTYCZNYMI
 
gene_diff <- function(dane1,dane2,thr,thr2) { # dane1- dane ekspresji po normalizacja dla grupy 1; dane2- dane ekspresji po normalizacja dla grupy 2, thr - próg dla p-wartoœci, thr2 - próg dla FC (niezlogarytmowanego)
  ## SPRAWDZENIE NORMALNOŒCI ROZK£ADU DANYCH (Test Lillieforsa)
  # Grupa 1:
  pLil1=matrix(data=NA,nrow=nrow(dane1),ncol=1)
  for (i in 1:nrow(dane1)) {
    temp1=lillie.test(dane1[i,]) # test
    pLil1[i,]=temp1$p.value # p-wartoœci 
  }
  # Obliczenie zawartoœci danych o rozk³adzie normalnym:
  l_norm1=length(which(pLil1[]>=0.05))
  zaw_norm1=l_norm1/nrow(pLil1)
  
  # Grupa 2:
  pLil2=matrix(data=NA,nrow=nrow(dane2),ncol=1)
  for (i in 1:nrow(dane2)) {
    temp2=lillie.test(dane2[i,]) # test
    pLil2[i,]=temp2$p.value # p-wartoœci 
  }
  # Obliczenie zawartoœci danych o rozk³adzie normalnym:
  l_norm2=length(which(pLil2[]>=0.05))
  zaw_norm2=l_norm2/nrow(pLil2)
  
#   ## HISTOGRAMY P-WARTOŒCI Z TESTU LILLEFORSA:
#   hist(pLil1[],main='Grupa 1',xlab='p-value') # grupa 1
#   hist(pLil2[],main='Grupa 2',xlab='p-value') # grupa 2
  
  ## SPRAWDZENIE JEDNORODNOŒCI WARIANCJI (F test)
  pF=matrix(data=NA,nrow=nrow(dane1),ncol=1) # Macierz wartoœci logicznych
  pF2=matrix(data=NA,nrow=nrow(dane1),ncol=1) # Macierz p-wartoœci
  for (i in 1:nrow(dane1)) {
    temp1=var.test(dane1[i,],dane2[i,]) # test F 
    pF2[i]=temp1$p.value # p-wartoœæ 
    if ((temp1$statistic>=temp1$conf.int[1])&&(temp1$statistic<=temp1$conf.int[2]))pF[i,1]=1
    else pF[i,1]=0
  }
  rownames(pF2)=rownames(dane1)
  rownames(pF)=rownames(dane1)

  # Obliczenie zawartoœci danych o jednorodnej wariancji:
  l_var=length(which(pF[]==1))
  zaw_var=l_var/nrow(pF)
  
  
  ## TEST t
  # Hipoteza zerowa: œrednie s¹ sobie równe.
  # Hipoteza alternatywna: œrednie nie s¹ sobie równe.
  # im mniejsza p-wartoœæ, tym wiêksze ró¿nice w ekspresji genów
  if ((zaw_norm1>=0.5)&& (zaw_norm2>=0.5)){ # jeœli co najmniej po³owa ka¿dej grupy ma rozk³ad normalny, to test t
    if (zaw_var>.5) {
      var.equal=T
      print("Test t")
    } else {
      var.equal=F
      print("Test Welcha") }
      
    pval_p=matrix(data=NA,nrow=nrow(dane1),ncol=1)
    for (i in 1:nrow(dane1)) {
      temp1=t.test(dane1[i,],dane2[i,],var.equal=var.equal)
      pval_p[i]=temp1$p.value # p-wartoœci
    }
    #pt.kor=mt.rawp2adjp(pt,proc=c("BH")) # korekcja Benjaminiego-Hochberga
    pval=p.adjust(pval_p,method="BH") # korekcja Benjaminiego-Hochberga
  
  
  ## TEST MANNA-WHITNEYA 
  # Hipoteza zerowa: mediany s¹ sobie równe.
  # Hipoteza alternatywna: mediany nie s¹ sobie równe.
  # im mniejsza p-wartoœæ, tym wiêksze ró¿nice w ekspresji genów
  } else { # test U Manna-Whitneya
    print('Test U Manna-Whitneya')
    pval_p=matrix(data=NA,nrow=nrow(dane1),ncol=1)
    for (i in 1:nrow(dane1)) {
      temp1=wilcox.test(dane1[i,],dane2[i,],adjust.method='BH')
      pval_p[i]=temp1$p.value # p-wartoœci
    }
    pval=p.adjust(pval_p,method="BH") # korekcja Benjaminiego-Hochberga
  }
  
  pval=cbind(pval)
  rownames(pval)=rownames(dane1)
  
   ### FOLD CHANGE
  ## UŒREDNIENIE KOLUMN
  mn1=rowMeans(dane1)
  mn2=rowMeans(dane2)

  ## FC
  FC=mn2-mn1
  changes=cbind(FC)
  row.names(changes)<-rownames(dane1)



  ## ZESTAWIENIE WYNIKÓW
  changes=cbind(changes,pval)
  il.pval=c()
  j=1;
  wek=seq(0,1,by=0.001)
  for (i in wek) {
    il.pval[j]=length(which(pval<=i)) # ile p-wartoœci jest mniejszych od i?
    j=j+1;
  }
  # Wykres dystrybuanty
  plot(wek,il.pval,type="l", col='violet', axes=T, ann=T,
       main='Wykres liczby genów o p-wartoœci mniejszej ni¿ poszczególne wartoœci',
       xlab="p-wartoœæ",
       ylab="Liczba genów", cex.lab=0.8, lwd=2)

  ### GENY RÓ¯NICUJ¥CE
  ## GENY O NAJNI¯SZEJ p-WARTOŒCI
  #thr=.2 # wartoœæ progu
  ind=which(pval<thr) # znalezienie genów o najni¿szej p-wartoœci (ni¿szej od progu)
  changes2=changes[ind[],] # najni¿sza p-wartoœæ
  
  #which(changes2[,1]>1) # sprawdzenie, czy dla któregoœ z istotnych statystycznie genów wzros³a ekspresja
  
  ## GENY O NAJMNIEJSZYM/NAJWIÊKSZYM FC
  ind2=which(abs(changes2[,1])>=thr2) # geny, dla których FC jest mniejszy od -(thr2) lub wiêkszy od (thr2)
  

  changes3=changes2[ind2[],] # wybrane geny ró¿nicuj¹ce (na podstawie p-wartoœci z testu i FC)
  ind3=ind[ind2] # indeksy genów ró¿nicuj¹cych
  
  # WARTOŒCI EKSPRESJI DLA GENÓW RÓ¯NICUJ¥CYCH
  exp_dif1=dane1[ind3[],] # grupa 1
  exp_dif2=dane2[ind3[],] # grupa 2


#   ### HISTOGRAMY (wartoœæ FC)
#   hist(changes2[,1],breaks=100,col='darkred',
#        main='Liczba genów o p-wartoœci mniejszej ni¿ próg w zale¿noœci od FC',
#        xlab='Wartoœæ FC',ylab='Liczba genów')
#   hist(changes3[,1],breaks=100,col='darkred',
#        main='Liczba genów wybranych jako ró¿nicuj¹ce w zale¿noœci od FC',
#        xlab='Wartoœæ FC',ylab='Liczba genów')
  
  # VOLCANO PLOT
  gene_list=data.frame(FC,pval)
  gene_list$threshold = as.factor((abs(gene_list$FC)>(thr2)) & (gene_list$pval < thr))
  
  ggplot(data=gene_list, aes(x=FC, y=-log10(pval),colour=threshold)) +
    geom_point(alpha=0.4, size=1.75) +
    theme(legend.position = "none") +
    xlim(c(-2, 2)) + ylim(c(0, 5)) +
    xlab("fold change") + ylab("-log10 p-value")+
    ggtitle("Volcano plot") + 
    theme(plot.title = element_text(lineheight=.8, face="bold"))


  ### SORTOWANIE W ZALE¯NOŒCI OD FC
  FCzm.sort=sort(changes3[,1],decreasing=FALSE,index.return=TRUE)
  ix=FCzm.sort$ix
  pvalzm.sort=changes3[ix,2]
  FCchanges.sort=FCzm.sort$x
  zm.sort=cbind(FCchanges.sort,pvalzm.sort)
  # Zapisanie do pliku
  write(row.names(zm.sort),file="Geny_roznicujace.txt")

  result<-list("expression1"=exp_dif1,"expression2"=exp_dif2)
  return (result)
}