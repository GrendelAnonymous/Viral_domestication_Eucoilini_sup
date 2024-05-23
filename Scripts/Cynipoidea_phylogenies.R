# Midpoint root a tree 
library("phytools")
library(ggplot2)
library(ggtree)





files <- list.files(path="/Users/bguinet/Desktop/Cynipoidea_paper/All_plot_trees", pattern="*.treefile", full.names=TRUE, recursive=FALSE)

Figure_number=0
#cluster<-"/Users/bguinet/Desktop/Cynipoid_paper/Cluster138_LbFVorf5.faa.ali.treefile"
#cluster<-"/Users/bguinet/Desktop/Cynipoid_paper/ALL_trees/Cluster98_redefined-Cluster795_redefined-Cluster1743_redefined_LbFVorf106.faa.ali.treefile"
cluster<-"/Users/bguinet/Desktop/Cynipoidea_paper/All_plot_trees/Cluster279_plot.aa.ali.treefile"  

#Open table and change label to fit label trees
tab<-read.table("/Users/bguinet/Desktop/Cynipoidea_paper/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event_dNdS.tab",sep=";",h=T)

tab<-tab[c('Cluster_hmmer','Names3','Mean_dNdS','Pvalue_dNdS','SE_dNdS','Prot_name','Best_hit_ORF_perc','Scaffold_score','Pseudogenized')]

tab$Names4<- gsub(" \\[ORF\\]","",tab$Names3)
tab$Names4<- gsub("\\(\\+\\):","_+__",tab$Names4)
tab$Names4<- gsub("\\(-\\):","_-__",tab$Names4)
tab$Names4<- gsub(":","_",tab$Names4)

for (cluster in files){

  tree_file<-cluster
  tree_name<-gsub(".*/","",tree_file)
  tree_name<-gsub("\\..*","",tree_name)
  tree_name<-gsub("_plot","",tree_name)
  
  subtab<-tab[tab$Cluster_hmmer==tree_name,]
  subtab$Best_hit_ORF_perc[is.na(subtab$Best_hit_ORF_perc)]<-0
  Prot_name<-unique(subtab$Prot_name)
  tree<-read.tree(tree_file)
  
  nb_taxa<-length(tree$tip.label)
  
  tree$tip.label<- gsub("-_Purifying_"," [Purifying]",tree$tip.label)
  tree$tip.label<- gsub("-_ORF_"," [ORF]",tree$tip.label)
  tree$tip.label<- gsub("-_VIRUS_"," [VIRUS]",tree$tip.label)
  tree$tip.label<- gsub("-_BACTERIA_"," [BACTERIA]",tree$tip.label)
  tree$tip.label<- gsub("-_EUKARYOTA_"," [EUKARYOTA]",tree$tip.label)
  tree$tip.label<- gsub("−_","-",tree$tip.label)

  
  #Change tiplabel plus strand
for (label in tree$tip.label){
    if (grepl('Rhoptromeris|Trybliographa|Thrichoplasta|Leptopilina',label)){
     label_save<-label
     label <- gsub(" \\[.*","",label)
     label<- gsub("\\(\\+\\):","_+__",label)
     label<- gsub("\\(-\\):","_-__",label)
     label<- gsub(":","_",label)
     label_score<-subtab$Scaffold_score[subtab$Names4==label]
     ORF_perc<-subtab$Best_hit_ORF_perc[subtab$Names4==label]
     Mean_dNdS<-subtab$Mean_dNdS[subtab$Names4==label]
     SE_dNdS<-subtab$SE_dNdS[subtab$Names4==label]
     Pvalue_dNdS<-subtab$Pvalue_dNdS[subtab$Names4==label]
     new_name<-(paste0(label," [",label_score,']'))
     Pseudo<-subtab$Pseudogenized[subtab$Names4==label]
     if (ORF_perc>0.70){
       new_name<-(paste0(new_name," [ORF]"))
     }else{
       if (Pseudo=="yes"){
         new_name<-(paste0(new_name," [*]"))
         }else{
       new_name<-new_name
       }
     }
     if (is.na(Mean_dNdS)==FALSE){
      if(Pvalue_dNdS<0.05 & Mean_dNdS+SE_dNdS<0.8 ){
        new_name<-(paste0(new_name," [Purifying]"))
      }else{
        new_name<-new_name
      }
     }else{
        new_name<-new_name
      }
     new_name<-gsub("_-__",'(-):',new_name)
     new_name<-gsub("_\\+__",'(+):',new_name)
     tree$tip.label[tree$tip.label==label_save]<-new_name
    }else{
  new_name<-gsub("_-__",'(-):',label)
  new_name<-gsub("_\\+__",'(+):',new_name)
  new_name<-gsub("___",'_',new_name)
  tree$tip.label[tree$tip.label==label]<-new_name
 }
}

  if (length(tree$tip.label[grepl("WSSV", tree$tip.label)]) >0){
    tryCatch(
      {
        tree<- root(tree, outgroup = tree$tip.label[grepl("WSSV", tree$tip.label)], resolve.root = TRUE)
      },
      error = function(e){
        print("no root")
      }
    )
    
  }else{
    tree<-midpoint.root(tree)
  }
  
  if (nb_taxa < 150){
    taxa_size=7
  }else{
    taxa_size=4
  }
  
  if (nb_taxa < 30){
    pdf_size_heigth=15
  }else{
    pdf_size_heigth=30
  }
  
  if (grepl("DNApol|lef-8|Cluster128",tree_name)){
    lim_right=30
  }else{
    lim_right=20
  }
  
  if (nb_taxa > 130){
    lim_right=40
  }
  
  if (nb_taxa > 150){
    pdf_size_heigth=70
  }
  
  p1<- ggtree(tree)
  
  p1$data$phylum <- "Virus"
  p1$data$phylum[grepl('EFV|boulardi|heterotoma|clavipes|Rhoptromeris|Dolichomitus|Chelonus|CiBV|MdBV|CcBV|Eurytoma|Melanaphis|FaBV|BtBV|VcBV|Trybliographa|Thrichoplasta|Phanerotoma',p1$data$label)]<-"Other EVEs"
  p1$data$phylum[grepl('boulardi|heterotoma|clavipes|Rhoptromeris|Trybliographa|Thrichoplasta',p1$data$label)]<-"Cynipids EVEs"
  p1$data$phylum[grepl('BACTERIA',p1$data$label)]<-"Bacteria"
  p1$data$phylum[grepl('EUKARYOTA',p1$data$label)]<-"Eukaryota"
  tree$tip.label[tree$tip.label %in% p1$data$label[p1$data$phylum=="Other EVEs" & grepl("VIRUS",p1$data$label)]]<-gsub(" \\[VIRUS\\]","",tree$tip.label[tree$tip.label %in% p1$data$label[p1$data$phylum=="Other EVEs" & grepl("VIRUS",p1$data$label)]])
  p1$data$label[p1$data$phylum=="Other EVEs" & grepl("VIRUS",p1$data$label)]<-gsub(" \\[VIRUS\\]","",p1$data$label[p1$data$phylum=="Other EVEs" & grepl("VIRUS",p1$data$label)])

  #p1$data$phylum[grepl('Bacteria',p1$data$label)]<-"Bacteria"
  #p1$data$phylum[grepl('nan−nan−nan',p1$data$label)]<-"NaN"
  #p1$data$phylum[grepl('Euk',p1$data$label)]<-"Eukaryota"
  
  # Continue only if filamentous viruses 
  
  mycolors=c("Virus"="#CC3636","Other EVEs"="#88BDD4","Cynipids EVEs"="#0CAEED","Bacteria"="#5B4B49","Eukaryota"="#42855B")
  
  title<-paste0(tree_name," - ",Prot_name)
    
  if (cluster=="/Users/bguinet/Desktop/Cynipoidea_paper/All_plot_trees/Cluster559-Cluster482-Cluster457-Cluster353-Cluster2746-Cluster2234-Cluster1338-Cluster1713-Cluster1176-Cluster1383-Cluster2222-Cluster2543-Cluster735-Cluster167_plot.aa.ali.treefile"){
    title="Lef-5"
  }
  
  p2 <- ggtree(tree, layout='rectangular') %<+% p1$data +
    geom_tippoint(
      mapping = aes(color = phylum),          # tip color by phyla. 
      size = 2,
      show.legend = FALSE) +
    #geom_nodelab(size = 4, col= "red",vjust=-.5)  +
    geom_nodelab(size = 5,nudge_x = 0.1)  +
    geom_tiplab(    mapping = aes(color = phylum),                      # adds name of phyla to tip of its branch
                    offset = 0.5,
                    size = taxa_size,
                    #geom = "label",
                    align = TRUE,
                    face = "bold",
                    #label.size = 5,
                    #label.padding = unit(padding, "lines"), # amount of padding around the labels
                    linetype = "dashed") +
    theme(
      axis.title.x = element_blank(), # removes x-axis title
      axis.title.y = element_blank(), # removes y-axis title
      legend.title = element_text(    # defines font size and format of the legend title
        face = "bold",
        size = 12),   
      legend.text=element_text(       # defines font size and format of the legend text
        face = "bold",
        size = 10),  
      plot.title = element_text(      # defines font size and format of the plot title
        size = 12,
        face = "bold"),  
      legend.position = "bottom",     # defines placement of the legend
      legend.box = "vertical",        # defines placement of the legend
      legend.margin = margin()) +
    scale_color_manual(values = mycolors)+
    scale_fill_manual(values = mycolors) + xlim(0,lim_right)  + ggtitle(title)
  
  p2
  ggsave(paste0("/Users/bguinet/Desktop/Cynipoidea_paper/Phylogenies_pdf/",paste0(tree_name," - ",Prot_name),".pdf"), width = 20, height = pdf_size_heigth, units = "in", limitsize = FALSE)
}

