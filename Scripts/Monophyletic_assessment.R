#!/usr/bin/env Rscript
library(ape)
library(phangorn)
library(seqinr)
library(tidytree)
library(ggplot2)
library(dplyr)
library(phytools)

args <- commandArgs(TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Please, add a file with one specie name / line and the path /path/Project/", call.=FALSE)
}

#Tree_file <- "/Users/bguinet/Desktop/Cynipoid_paper/Cynipoidea_phylogeny2.nwk" #The Species phylogenetic tree
#Cluster_directory <- "/Users/bguinet/Desktop/Cynipoid_paper/Cluster_phylogeny_dNdS_filamentous/" #The directory where are found all the cluster file with the extension .treefile
#monophylyprop <- 0.8 # The monophylyprop option determines how monophyletic a species group in the cluster tree should be in the species tree. 1 = 100% monoph. 0.6 = at least 60% monoph.
#Output_table_file <- "/Users/bguinet/Desktop/Cynipoid_paper/Cluster_phylogeny_dNdS_filamentous/Monophyletic_tab.tab" #The Outpute file with full directory
#Blast_table_scores <- "/Users/bguinet/Desktop/Cynipoid_paper/dsDNA_hmmer_clusters_cov_NR_ORF_euk_repeat_scores_synteny.tab" #The table with scores
#taxo<-read.table("/Users/bguinet/Desktop/Cynipoid_paper/Species_sub_sup_families_order.txt", header=T,sep=";")

Tree_file <-args[1]
Cluster_directory <- args[2]
monophylyprop <- args[3]
Blast_table_scores <- args[4]
taxo <- args[5]
Output_table_file <- args[6]



print (paste0(" Treefile : ",Tree_file))
print (paste0(" Cluster_directory : ",Cluster_directory))
print (paste0(" monophylyprop : ",monophylyprop))
print (paste0(" Blast_table_scores : ",Blast_table_scores))
print (paste0(" Taxo table : ",taxo))
print (paste0(" Output_table_file : ",Output_table_file))
Tree_file ="/beegfs/data/bguinet/Cynipoidea_project/Monophyletic_assessment/Cynipoidea_phylogeny.nwk"
Cluster_directory="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_phylogeny/"
monophylyprop=0.8
Blast_table_scores="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score.tab"
taxo="/beegfs/data/bguinet/Cynipoidea_project/Monophyletic_assessment/Species_sub_sup_families_order.txt"
Output_table_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Monophyletic_assessment/Monophyletic_table.tab"

# Example usage :
#/beegfs/data/soft/R-4.0.5/bin/Rscript /beegfs/data/bguinet/Cynipoidea_project/Monophyletic_assessment.R {input.treefile} {input.Cluster_directory} {input.monophylyprop} {input.Blast_table_scores} {input.taxo_table} {output.Output_table_file}


taxo<-read.table(taxo, header=T,sep=";")

sptr<-read.tree(Tree_file)
sptr$node.label[which(sptr$node.label=="")]<-"100/100"

#Function to test inclusion of one partition inside another

isincludedin<-function(arr, clan) {
  arr1<-clan[arr[1],]
  arr2<-clan[arr[2],]
  return(min(arr2-arr1)!="-1")
}

##transform species tree int a matrix
spmat<-t(do.call(cbind, lapply(lapply(subtrees(sptr), function(x) x$tip.label),function(a,b) is.element(b,a),b=sptr$tip.label)))+0
colnames(spmat)<-sptr$tip.label
rownames(spmat)<-sptr$node.label
mat2<-matrix(0,ncol=Ntip(sptr), nrow=Ntip(sptr))
diag(mat2)<-1
rownames(mat2)<-rep("100/100", Ntip(sptr))
spmat<-rbind(mat2, spmat)
##get all clustr names

clus_1<-list.files(path = Cluster_directory, pattern = "*.treefile")
clus<-c()
for (i in clus_1){
   if ( grepl('dNdS', i)){
    clus<- c(clus,i) 
  }
}

testinspmat<-function(sps,spmat) {
  max(apply(spmat[,sps],1,sum))
}

filtermonoph<-function(monoph,clans,spmat,thres,cluster) {
  sps<-unique(names(which(clans[monoph,]==1)))
  print(sps)
  sub<-which(apply(spmat[,sps,drop=F],1,sum)==length(sps))
  nbdesccontainingsps<-apply(spmat[sub,,drop=F],1,sum)
  #Allow a lower threshold when lonly 4 taxa within the monophyletic group
  small_nbdesccontainingsps <- nbdesccontainingsps[nbdesccontainingsps<5]
  prop<-length(sps)/nbdesccontainingsps
  small_prop <- length(sps)/small_nbdesccontainingsps
  if(length(small_prop[small_prop>0.60]>=0.60)>=1){
    if(length(small_prop[small_prop>0.80]>=0.60)==0){
      print(paste0(sps,cluster))
    }
  }
  if(length(small_prop[small_prop>0.60]>=0.60)>=1){
    return(max(small_prop)>=0.60)
  }else{
    return(max(prop)>=thres)
  }
}

Blast_table<-read.csv(Blast_table_scores ,sep=";",h=T)
Blast_table$label<- gsub("):","__",Blast_table$Names3)
Blast_table$label<- gsub("\\(","_",Blast_table$label)
Blast_table$label<- gsub(":","_",Blast_table$label)
Blast_table$label<- gsub("\\+","_",Blast_table$label)
Blast_table$label<- gsub("____","_+__",Blast_table$label)
Blast_table$label<- gsub(" \\[ORF\\]","",Blast_table$label)
Blast_table$label<- gsub("\\.1\\.1",".1",Blast_table$label)
#Blast_table$label<- gsub("GCA_009025955","GCA_009025955.1",Blast_table$label)
Blast_table$Scaffold_score[Blast_table$label=="C2179507_235-390_-__Trybliographa"]<-"C"
Blast_table$Scaffold_score["Leptopilina_boulardi_GCA_015476485.1" %in% Blast_table$Species_name]<-"C"
#Blast_table$Scaffold_score["Leptopilina_heterotoma_GCA_015476485.1" %in% Blast_table$Species_name]<-"C"
#Remove duplicate keep first 
Blast_table<-Blast_table %>% distinct(label, .keep_all = TRUE)

#Manually stuff 
#Blast_table<-Blast_table[!grepl("5480-5854",Blast_table$label),]
#Blast_table<- Blast_table[!Blast_table$label %in% c("5480-5854"),]

cluster_rerun<-c('')
level<-"Family"
GetClades<-function(cluster="Cluster", level="Family",monophylyprop=monophylyprop,otherwise = NA) {
  print(cluster)
  #root the tree and save it 
  gntr<-read.tree(paste0(Cluster_directory,cluster))
  #gntr$label<-gsub("GCA_015476485","GCA_015476485.1",gntr$label)
  #gntr$label<-gsub("GCA_009025955","GCA_009025955.1",gntr$label)
  gntr<-midpoint.root(gntr)
  write.tree(gntr,file=paste0(Cluster_directory,cluster))
  #print(monophylyprop)
  #print(level)
  #ggtree(gntr) + geom_tiplab()  +geom_nodelab(size = 7, col= "red") 
  #read gene tree
  #gntr<-read.tree(paste0(Cluster_directory,cluster))
  #gntr<- unroot( gntr)
  gntr2<-as_tibble(gntr)
  #gntr2$label<-gsub("GCA_015476485","GCA_015476485.1",gntr2$label)
  #gntr2$label<-gsub("GCA_009025955","GCA_009025955.1",gntr2$label)
  #Merge with scaffold score
  gntr2<-left_join(fortify(gntr2), fortify(select(Blast_table,'label','Scaffold_score')), by=c('label'))
  gntr2$label<-gsub("___PoEFV","_+__PoEFV",gntr2$label)
  gntr2$label<-gsub("_-_PoEFV","_-__PoEFV",gntr2$label)
  toMatch <- c("A","B","C","D")
  if (length(unique (grep(paste(toMatch,collapse="|"), gntr2$Scaffold_score, value=TRUE)))>=1) {    
    gntr2$label<-paste0(gntr2$Scaffold_score,"_",gntr2$label)
    gntr2$label<-gsub("NA_","",gntr2$label)
    gntr2<-left_join(as_tibble(gntr, ladderize = FALSE),select(gntr2,-c(,"parent",'branch.length','Scaffold_score')) , by = c("node" = "node"))
    gntr2$label.x<-NULL
    names(gntr2)[names(gntr2)=="label.y"] <- "label"
    #If X to E, then remove __ to _ in order to not pass the filter
    nb <- 1
    label_list<-c()
    for (i in gntr2$label){
      if (nchar(i) > 4) {
        label_list[[nb]] <-i
        nb <- nb + 1
      }
    }
    i <- grepl( "(X_|F_|E_)", label_list )
    label_list <- gsub( "____", "_+__", label_list )
    label_list[i] <- gsub( "__", "_", label_list[i] )
    gntr$tip.label <-label_list 
    #Remove letters 
    gntr$tip.label<-sub("^[A-X]_", "",gntr$tip.label)
    #prepare clans
    cl2<-getClans(gntr)
    #add rownames if missing
    if (is.null(rownames(cl2))) rownames(cl2)<-rep("-1", nrow(cl2))
    #	sizebipar<-apply(cl2,1,sum)
    nam<-colnames(cl2)
    spnam<-unlist(lapply(strsplit(colnames(cl2), "__"), function(x) x[2]))
    if (sum(!is.na(spnam))>0){
      #rename clans based on the taxonomic level chosen
      nam2<-taxo[,level][match(spnam,taxo$Species)]
      nam2<-as.character(nam2)
      nam2[is.na(nam2)]<-"Virus"
      cl3<-cl2
      cl4<-cl2
      colnames(cl3)<-nam2
      colnames(cl4)<-spnam
      testMonophyly<-apply(cl3,1, function(x) unique(names(x[x==1])))
      monophgroups<-which(unlist(lapply(testMonophyly, function(x) (length(x)==1)&x[1]!="Virus")))
      ##We only keep monoph groups that are monophyletic (or almost, depending on monophylyprop)
      monophgroups<-monophgroups[sapply(monophgroups, filtermonoph, clans=cl4,spmat=spmat,thres=monophylyprop,cluster=cluster)] #issue 
      if (length(monophgroups)==1) clades<-monophgroups 
      ##now we remove groups included in others, to only keep the largest ones.
      else {
        allcompare<-cbind(rep(monophgroups, each=length(monophgroups)), rep(monophgroups, length(monophgroups)))
        allcompare<-allcompare[(allcompare[,1]-allcompare[,2])!=0,] ##remùove self comparison
        getinclusions<-apply(allcompare,1,isincludedin,clan=cl2) #is group of column 1 included in group of colmumn two
        toremove<-unique(allcompare[getinclusions,1])
        tokeep<-setdiff(monophgroups, toremove)
        clades<-monophgroups[match(tokeep,monophgroups)]
      }
      
      #PREPARE DATA FOR OUTPUT
      Cluster<-rep(strsplit(cluster, ".fa")[[1]][1], length(clades))
      Cluster<-gsub("\\..*","",Cluster)
      Cluster_name<-strsplit(cluster, ".fa")[[1]][1]
      Cluster_name<-gsub("\\..*","",Cluster_name)
      #Cluster_name<-sub('\\..*', "", Cluster_name)
      Event<-1:length(clades) #as YOU did
      Nloc<--array() #nb of species in each monoph clade
      Nsp<-array()
      Nsp_MRCA<-array() #size of the smallest clade containing these species in the species tree
      Tips4dNdS<-apply(cl2[clades,,drop=FALSE],1,function(x) names(x[x==1])) #C'est une liste, plus pratique pour caluler ensuite le dNdS?
      
      ##CALCUL DNDS.
      #alignment_file <- paste("/beegfs/data/bguinet/M2/dNdS_calculation/",Cluster_name,".codon.fa.aln",sep="")
      #alignment<-read.alignment(alignment_file,'fasta')
      #kaks_matrix<-kaks(alignment, verbose =F, debug = FALSE, forceUpperCase = TRUE)
      #omega<-kaks_matrix$ka/kaks_matrix$ks
      #omega<-as.matrix(omega)
      #write.table(omega, file=paste0("/beegfs/data/bguinet/M2/dNdS_calculation/",Cluster_name,"_kaks_matrix"), sep = "\t")
      
      Tips4dNdS.array<-apply(cl2[clades,,drop=FALSE],1,function(x) paste("[",paste(paste("'",names(x[x==1]),"'",sep=""),collapse=","),"]",sep="")) #idem en mode caractères
      Boot<-names(clades) #tout simplement.
      Family<-unlist(testMonophyly[clades])
      for (i in 1:length(clades)) {
        speciesinclade<-spnam[cl2[clades[i],]==1]
        Nsp[i]<-length(unique(speciesinclade))
        Nloc[i]<-length(speciesinclade)
        if (Nsp[i]>1) {
          Nsp_MRCA[i]<-length(Descendants(sptr, getMRCA(sptr,unique(speciesinclade)), type="tips")[[1]])
        }
        else {
          Nsp_MRCA[i]<-1 #should be 1?
        }
      }	
      RES<-cbind(cluster, Event, Nloc,Nsp, Nsp_MRCA,Tips4dNdS.array, Boot)
      rownames(RES)<-NULL
      return(RES)
    }
    else {
      RES<-NULL
      
      RES<-cbind(cluster, NA, NA,NA, NA,NA, NA)
      return(RES)
    }
  }
}

#setwd(Cluster_directory)
#clus<-list.files(pattern = "\\.treefile")


clus<-clus[!grepl(".bak",clus)]
clus<-clus[!grepl(".pruned",clus)]


#clus

RESULT<-NULL
for (i in 1:length(clus)) {
  #print(paste(i,clus[i],sep=" - "))
  RESULT<-rbind(RESULT, GetClades(clus[i],"Family",monophylyprop))
}

RESULT<- as.data.frame(RESULT)

write.table(RESULT,Output_table_file,sep=";")

