### POR?WNANIE DANYCH METODAMI STATYSTYCZNYMI

gene_diff <- function(dane1,dane2,thr,thr2) { # dane1- dane ekspresji po normalizacja dla grupy 1; dane2- dane ekspresji po normalizacja dla grupy 2, thr - pr?g dla p-warto?ci, thr2 - pr?g dla FC (niezlogarytmowanego)
  ## SPRAWDZENIE NORMALNO?CI ROZK?ADU DANYCH (Test Lillieforsa)
  # Grupa 1:
  pLil1=matrix(data=NA,nrow=nrow(dane1),ncol=1)
  for (i in 1:nrow(dane1)) {
    temp1=lillie.test(dane1[i,]) # test
    pLil1[i,]=temp1$p.value # p-warto?ci 
  }
  # Obliczenie zawarto?ci danych o rozk?adzie normalnym:
  l_norm1=length(which(pLil1[]>=0.05))
  zaw_norm1=l_norm1/nrow(pLil1)
  
  # Grupa 2:
  pLil2=matrix(data=NA,nrow=nrow(dane2),ncol=1)
  for (i in 1:nrow(dane2)) {
    temp2=lillie.test(dane2[i,]) # test
    pLil2[i,]=temp2$p.value # p-warto?ci 
  }
  # Obliczenie zawarto?ci danych o rozk?adzie normalnym:
  l_norm2=length(which(pLil2[]>=0.05))
  zaw_norm2=l_norm2/nrow(pLil2)
  
  #   ## HISTOGRAMY P-WARTO?CI Z TESTU LILLEFORSA:
  #   hist(pLil1[],main='Grupa 1',xlab='p-value') # grupa 1
  #   hist(pLil2[],main='Grupa 2',xlab='p-value') # grupa 2
  
  ## SPRAWDZENIE JEDNORODNO?CI WARIANCJI (F test)
  pF=matrix(data=NA,nrow=nrow(dane1),ncol=1) # Macierz warto?ci logicznych
  pF2=matrix(data=NA,nrow=nrow(dane1),ncol=1) # Macierz p-warto?ci
  for (i in 1:nrow(dane1)) {
    temp1=var.test(dane1[i,],dane2[i,]) # test F 
    pF2[i]=temp1$p.value # p-warto?? 
    if ((temp1$statistic>=temp1$conf.int[1])&&(temp1$statistic<=temp1$conf.int[2]))pF[i,1]=1
    else pF[i,1]=0
  }
  rownames(pF2)=rownames(dane1)
  rownames(pF)=rownames(dane1)
  
  # Obliczenie zawarto?ci danych o jednorodnej wariancji:
  l_var=length(which(pF[]==1))
  zaw_var=l_var/nrow(pF)
  
  
  ## TEST t
  # Hipoteza zerowa: ?rednie s? sobie r?wne.
  # Hipoteza alternatywna: ?rednie nie s? sobie r?wne.
  # im mniejsza p-warto??, tym wi?ksze r??nice w ekspresji gen?w
  if ((zaw_norm1>=0.5)&& (zaw_norm2>=0.5)){ # je?li co najmniej po?owa ka?dej grupy ma rozk?ad normalny, to test t
    if (zaw_var>.5) {
      var.equal=T
      print("Test t")
    } else {
      var.equal=F
      print("Test Welcha") }
    
    pval_p=matrix(data=NA,nrow=nrow(dane1),ncol=1)
    for (i in 1:nrow(dane1)) {
      temp1=t.test(dane1[i,],dane2[i,],var.equal=var.equal)
      pval_p[i]=temp1$p.value # p-warto?ci
    }
    #pt.kor=mt.rawp2adjp(pt,proc=c("BH")) # korekcja Benjaminiego-Hochberga
    pval=p.adjust(pval_p,method="BH") # korekcja Benjaminiego-Hochberga
    
    
    ## TEST MANNA-WHITNEYA 
    # Hipoteza zerowa: mediany s? sobie r?wne.
    # Hipoteza alternatywna: mediany nie s? sobie r?wne.
    # im mniejsza p-warto??, tym wi?ksze r??nice w ekspresji gen?w
  } else { # test U Manna-Whitneya
    print('Test U Manna-Whitneya')
    pval_p=matrix(data=NA,nrow=nrow(dane1),ncol=1)
    for (i in 1:nrow(dane1)) {
      temp1=wilcox.test(dane1[i,],dane2[i,],adjust.method='BH')
      pval_p[i]=temp1$p.value # p-warto?ci
    }
    pval=p.adjust(pval_p,method="BH") # korekcja Benjaminiego-Hochberga
  }
  
  pval=cbind(pval)
  rownames(pval)=rownames(dane1)
  
  ### FOLD CHANGE
  ## U?REDNIENIE KOLUMN
  mn1=rowMeans(dane1)
  mn2=rowMeans(dane2)
  
  ## FC
  FC=mn2-mn1
  changes=cbind(FC)
  row.names(changes)<-rownames(dane1)
  
  
  
  ## ZESTAWIENIE WYNIK?W
  changes=cbind(changes,pval)
  il.pval=c()
  j=1;
  wek=seq(0,1,by=0.001)
  for (i in wek) {
    il.pval[j]=length(which(pval<=i)) # ile p-warto?ci jest mniejszych od i?
    j=j+1;
  }
  # Wykres dystrybuanty
  plot(wek,il.pval,type="l", col='violet', axes=T, ann=T,
       main='Wykres liczby gen?w o p-warto?ci mniejszej ni? poszczeg?lne warto?ci',
       xlab="p-warto??",
       ylab="Liczba gen?w", cex.lab=0.8, lwd=2)
  
  ### GENY RӯNICUJ?CE
  ## GENY O NAJNI?SZEJ p-WARTO?CI
  #thr=.2 # warto?? progu
  ind=which(pval<thr) # znalezienie gen?w o najni?szej p-warto?ci (ni?szej od progu)
  changes2=changes[ind[],] # najni?sza p-warto??
  
  #which(changes2[,1]>1) # sprawdzenie, czy dla kt?rego? z istotnych statystycznie gen?w wzros?a ekspresja
  
  ## GENY O NAJMNIEJSZYM/NAJWI?KSZYM FC
  ind2=which(abs(changes2[,1])>=thr2) # geny, dla kt?rych FC jest mniejszy od -(thr2) lub wi?kszy od (thr2)
  
  
  changes3=changes2[ind2[],] # wybrane geny r??nicuj?ce (na podstawie p-warto?ci z testu i FC)
  ind3=ind[ind2] # indeksy gen?w r??nicuj?cych
  
  # WARTO?CI EKSPRESJI DLA GEN?W RӯNICUJ?CYCH
  exp_dif1=dane1[ind3[],] # grupa 1
  exp_dif2=dane2[ind3[],] # grupa 2
  
  
  #   ### HISTOGRAMY (warto?? FC)
  #   hist(changes2[,1],breaks=100,col='darkred',
  #        main='Liczba gen?w o p-warto?ci mniejszej ni? pr?g w zale?no?ci od FC',
  #        xlab='Warto?? FC',ylab='Liczba gen?w')
  #   hist(changes3[,1],breaks=100,col='darkred',
  #        main='Liczba gen?w wybranych jako r??nicuj?ce w zale?no?ci od FC',
  #        xlab='Warto?? FC',ylab='Liczba gen?w')
  
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
  
  
  ### SORTOWANIE W ZALE?NO?CI OD FC
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