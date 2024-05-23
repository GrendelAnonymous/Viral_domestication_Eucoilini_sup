library("msa")
library("texi2pdf")
#install.packages("texi2pdf")


setwd("/Users/bguinet/Desktop/Cynipoidea_paper/All_MSA_files/")

files <- list.files(path="/Users/bguinet/Desktop/Cynipoidea_paper/All_MSA_files", pattern="*.aa", full.names=TRUE, recursive=FALSE)
tot<-length(files)
n<-1
for (aln in files){

  ## call unified interface msa() for default method (ClustalW) and
  ## default parameters
  mySequences <- readAAStringSet(aln)
  Cluster_name<- gsub(".*\\/","",aln)
  Cluster_name<- gsub(".aa","",Cluster_name)
  Cluster_name<- gsub("_MSA_","",Cluster_name)
  Cluster_name<- gsub("_MSA","",Cluster_name)
  myAlignment <-  msaClustalOmega(mySequences,order="input")

  
  pdf_file = paste0("/Users/bguinet/Desktop/Cynipoidea_paper/All_MSA_files/",Cluster_name,"_MSA.pdf",sep="")
  
  ## create PDF file according to some custom settings
  #tmpFile <- tempfile(pattern="msa", tmpdir="/Users/bguinet/Desktop/Cynipoidea_paper", fileext=".pdf")
  #tmpFile
  msaPrettyPrint(myAlignment, file=pdf_file, output="pdf",
                 showNames="left", showNumbering="none", showLogo="top", 
                 logoColors="structure", shadingMode="functional",
                 shadingModeArg="hydropathy",
                 verbose=FALSE, askForOverwrite=FALSE , paperWidth = 12,paperHeight = 20)
  print(paste0(n,"/",tot))
  n=n+1
}

n=1
for (aln in files){
  mySequences <- readAAStringSet(aln)
  Cluster_name<- gsub(".*\\/","",aln)
  Cluster_name<- gsub(".aa","",Cluster_name)
  Cluster_name<- gsub("_MSA_","",Cluster_name)
  print(paste0("includepdf[pages={1-},scale=0.75]{",Cluster_name,"_MSA.pdf}"))
  #print(paste0("Figure S",n," Alignment of ",Cluster_name," and their homologs in Cynipoidea."))
  print("%")
  n=n+1
}


