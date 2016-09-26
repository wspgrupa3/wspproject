get_expset = function(filename){#data file
  
#Marta Zarebska
  library('affy')
  library('Biobase')
  library('hgu95av2.db')
  library('gahgu95av2.db')
  
  ##### WCZYTANIE datasetA_scans #####
  dataScanFile <- filename #file.path(getwd(), filename)
  ##### STOWRZENIE phenoData #####
  pData <- read.table(filename, header=TRUE, sep="\t")
  summary(pData) # podsumowanie informacji o pr?bkach
  
  
  metadata <- data.frame(labelDescription=c("Anotacja","Tkanka","Pr?bka","Nazwa pliku CEL"),row.names=colnames(pData))
  
  phenoData=new("AnnotatedDataFrame",data=pData, varMetadata=metadata) # Annotated Data Frame z datainf
  #anot_datainf2=read.AnnotatedDataFrame("datasetA_scans.txt", header=TRUE, sep="\t")
  
  
  ##### WCZYTANIE DANYCH MIKROMACIERZOWYCH #####
  filePaths=paste(pData$scan,'.CEL',sep='') # nazwy plik?w CEL
  assign('varlookup', filePaths, .GlobalEnv);
  data_raw=ReadAffy(filenames=filePaths,phenoData=phenoData) # wczytanie, powstaje obiektpDataFile <- file.path(dataDirectory, "pData.txt")
  data_raw@annotation=paste("ga",data_raw@annotation,sep="") # zmiana anotacji
  data_raw@cdfName=paste("GA",data_raw@cdfName,sep="") # zmiana anotacji
  
  ##### NORMALIZACJA RMA #####
  RMA=expresso(data_raw,bgcorrect.method="rma",normalize.method="quantiles",
               pmcorrect.method="pmonly",summary.method="medianpolish")
  data=exprs(RMA)
  colnames(data)=pData$scan
  
  ##### FEATURE DATA ####
  # fData=matrix(,nrow(data),2)
  fData=data.frame(Gene=c(1:nrow(data)),Description=c(1:nrow(data)),row.names=rownames(data)) # przyk?adowe warto?ci
  # colnames(fData)=c("Gene", "Description")
  # rownames(fData)=rownames(data)
  
  metadata2 <- data.frame(labelDescription=c("Nazwa genu","Opis genu"),row.names=colnames(fData))
  
  featureData=new("AnnotatedDataFrame",data=fData, varMetadata=metadata2) # Annotated Data Frame z datainf
  
  
  
  ##### OPIS EKSPERYMENTU #####
  experimentData=new("MIAME",
                     name="Classification of Human Lung Carcinomas by mRNA Expression Profiling Reveals Distinct Adenocarcinoma Sub-classes",
                     lab="Harvard Medical School",
                     contact="Arindam_Bhattacharjee@dfci.harvard.edu, staunton@genome.wi.mit.edu, Matthew_Meyerson@dfci.harvard.edu, golub@genome.wi.mit.edu",
                     title="Classification of Human Lung Carcinomas by mRNA Expression Profiling Reveals Distinct Adenocarcinoma Sub-classes",
                     abstract="We have generated a molecular taxonomy of lung carcinoma, the leading cause of cancer death in the United States and worldwide. Using oligonucleotide microarrays, we analyzed mRNA expression levels corresponding to 12,600 transcript sequences in 186 lung tumor samples, including 139 adenocarcinomas resected from the lung. Hierarchical and probabilistic clustering of expression data defined distinct sub-classes of lung adenocarcinoma. Among these were tumors with high relative expression of neuroendocrine genes and of type II pneumocyte genes, respectively. Retrospective analysis revealed a less favorable outcome for the adenocarcinomas with neuroendocrine gene expression. The diagnostic potential of expression profiling is emphasized by its ability to discriminate primary lung adenocarcinomas from metastases of extra-pulmonary origin. These results suggest that integration of expression profile data with clinical parameters could aid in diagnosis of lung cancer patients",
                     url="http://www.broadinstitute.org/mpr/lung/")
  
  
  ##### STWORZENIE OBIEKTU ExpressionSet #####
  rownames(phenoData)=colnames(data)
  expset=ExpressionSet(assayData=data,phenoData=phenoData,
                       experimentData=experimentData,
                       annotation=data_raw@annotation,featureData=featureData)
  
  
 return(expset); 
}