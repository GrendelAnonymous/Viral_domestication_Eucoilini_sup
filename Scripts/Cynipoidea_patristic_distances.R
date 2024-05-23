

tree<-read.tree("/Users/bguinet/Desktop/Cynipoidea_paper/Concatenated_sequences_main_datation.treefile")
distance_matrix<-cophenetic.phylo(tree)
patristic_tab<-as.data.frame(as.table(distance_matrix))
colnames(patristic_tab)<-c("C1","C2","Distance")

patristic_tab<- patristic_tab[!patristic_tab$C1 == patristic_tab$C2,]
patristic_tab<-patristic_tab[!duplicated(patristic_tab$Distance),]


Microgastroid<-c("MdENV","CcENV","TnENV","Chelonus_insularis","CiENV")
Cynipids<-c("LbEFV","LcEFV","LhEFV","ThEFV","TrEFV","RhEFV")
Viruse<-c("CoBV","AcMNPV","LdMNPV","CpV","NeseNPV","CuniNPV","DhNV","PmNV","HzNV-1","PhENV","ToNV","EbENV","FaENV","OrNV","tom",    
          "NlENV","VcENV","BtENV","GbNV","MsENV","ApENV","DnENV","PoEFV","PoFV","CcFV1","CcFV2","DFV","LhFV","LbFV","EfFV",
          "PcFV","GpSGHV","MdSGHV")


Microgastroid_tab<-patristic_tab[patristic_tab$C1 %in% Microgastroid & patristic_tab$C2 %in% Microgastroid,]
Microgastroid_tab$Type <-"Microgastroids"
Cynipids_tab<-patristic_tab[patristic_tab$C1 %in% Cynipids & patristic_tab$C2  %in% Cynipids ,]
Cynipids_tab$Type <-"Cynipids"
Viruses_tab <- patristic_tab[patristic_tab$C1 %in% Viruse & patristic_tab$C2 %in% Viruse,]
Viruses_tab$Type <- "Viruses"

All_tab<-rbind(Microgastroid_tab,Cynipids_tab,Viruses_tab)

library("ggpubr")
my_comparisons <- list( c("Microgastroids", "Cynipids"), c("Cynipids", "Viruses"), c("Microgastroids", "Viruses") )
plot1<-ggplot(All_tab, aes(x = Type, y = Distance)) + 
  geom_boxplot(aes(fill = Type))+ theme_bw()+
  ylab("Patristic distances") + 
  stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6)   +
  theme(axis.title.x=element_blank())

plot1
ggsave(paste0("/Users/bguinet/Desktop/Cynipoidea_paper/Patristic_distance_boxplot.pdf"), width = 6, height =6, units = "in", limitsize = FALSE)


#####
## 

# Load all 16 gene phylogenies 
library(ape)

files <- list.files(path="/Users/bguinet/Desktop/Cynipoidea_paper/All_plot_trees", pattern="*.treefile", full.names=TRUE, recursive=FALSE)

dat<-as.data.frame(matrix(ncol=3, nrow=0))

files<-files[!grepl("Cluster29",files)]
i<-1
All_relative_distances<-c()
for (cluster in files){

  tree<-read.tree(cluster)

  if (length(tree$tip.label[grepl( c("Leptopilina"),tree$tip.label) | 
                            grepl( c("Thrichoplasta"),tree$tip.label ) | 
                            grepl( c("Rhoptromeris"),tree$tip.label ) | 
                            grepl( c("Trybliographa"),tree$tip.label )]) >= 6 ){
      print(cluster)
      distance_matrix<-cophenetic.phylo(tree)
      patristic_tab<-as.data.frame(as.table(distance_matrix))
      colnames(patristic_tab)<-c("C1","C2","Distance")
      
      patristic_tab<- patristic_tab[!patristic_tab$C1 == patristic_tab$C2,]
      #patristic_tab<-patristic_tab[!duplicated(patristic_tab$Distance),]
      
      # Remove second event scaffolds
      
      patristic_tab<-patristic_tab[!grepl("scaffold39978|scaffold10766|scaffold7282|scaffold49329|scaffold98928|scaffold52008|scaffold51610",patristic_tab$C1),]
      patristic_tab<-patristic_tab[!grepl("scaffold39978|scaffold10766|scaffold7282|scaffold49329|scaffold98928|scaffold52008|scaffold51610",patristic_tab$C2),]
      
      a<-patristic_tab$Distance[grepl("boulardi",patristic_tab$C1) & grepl("Trybliographa",patristic_tab$C2) | 
                               grepl("boulardi",patristic_tab$C2) & grepl("Trybliographa",patristic_tab$C1)  ]
 
      Mean_Intra_distances <- unique(a)
      
  
      
      a<-unique(patristic_tab$Distance[grepl("LhFV",patristic_tab$C1) & grepl("Trybliographa",patristic_tab$C2) | 
                                  grepl("LhFV",patristic_tab$C2) & grepl("Trybliographa",patristic_tab$C1)])
      
      b<-unique(patristic_tab$Distance[grepl("LhFV",patristic_tab$C1) & grepl("boulardi",patristic_tab$C2) | 
                                  grepl("LhFV",patristic_tab$C1) & grepl("boulardi",patristic_tab$C2)  ])
      
      Mean_Extra_distances <- mean(c(a,b))
      
      relative_distance=Mean_Intra_distances/Mean_Extra_distances
      
      dat[i,1] = relative_distance
      cluster_name<-gsub("_plot.aa.ali.treefile","",gsub(".*\\/",'',cluster))
      dat[i,2] = gsub("-.*","",cluster_name)
      dat[i,3] = "Event 1"
      i<-i+1

      
      
      #All_relative_distances<- append(All_relative_distances,relative_distance)

  }
}

library(ggrepel)
dat %>% 
  ggplot(aes(x=V3,y=V1, label = V2))+
  geom_boxplot(width=.5)+
  geom_text_repel(arrow = arrow(length = unit(0.02, "npc")),box.padding = 1)+
  theme(legend.position="none")+  theme_classic()

