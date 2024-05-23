


>NW_025111203.1:2150605-2150865(+):Leptopilina_heterotoma_GCA_015476425.1
MKAKGSYISFHDEFEYIAGGQSVRDSLSNFLFFDYGGTCTKSQCPGKNGGSFTTYRFSTSNIIIYRQEIKKKKLHCKSNILQKYFYY

>JADEYJ010000123.1:877336-878082(-):Leptopilina_boulardi_GCA_019393585.1
MTNESSAISQQPDHKVYQLREFSLIYTIKWLGTFDEENEIFTINPKILENCQIVVDKSTLEFETIFLRIFQLINPSLLCGIKLDKNVIIVPTLSMRAKGSFISFHDEFEYISGGQTVRDSLSNFVFFDYGGTFTKSQCPGKNGGSFTTYRFLTSNMIVYRQEIKKNNCTVIQTYNKNIFNINNRRTAKLFLFRNNQFFNIVTQSEISYELYQPFFTVFKIVWSNYKRSEDGRMEIWPTIEIIAACVQLI

>JABAIF010000692.1:288490-289236(+):Leptopilina_boulardi_GCA_015476485.1
MTNESSAISQQPDHKVYQLREFSLIYTIKWLGTFDEENEIFTINPKILENCQIVVDKSTLEFETIFLRIFQLINPSLLCGIKLDKNVIIVPTLSMRAKGSFISFHDEFEYISGGQTVRDSLSNFVFFDYGGTFTKSQCPGKNGGSFTTYRFLTSNMIVYRQEIKKNNCTVIQTYNKNIFNINNRRTAKLFLFRNNQFFNIVTQSEISYELYQPFFTVFKIVWSNYKRSEDGRMEIWPTIEIIAACVQLI


>JABAIF010007716.1:423933-424160(-):Leptopilina_boulardi_GCA_015476485.1
IGFVKSEIRKGVLPRMLSEILNTRLMVKKAMKDHSKDDRSLQQVLHNRQLGLKLIANVTYGYTSANFSGRMSCIEV---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

>JABAIF010000692.1:288490-289236(+):Leptopilina_boulardi_GCA_015476485.1
ATGACAAATGAAAGTAGTGCTATTTCTCAGCAACCCGATCACAAAGTATATCAATTGCGTGAATTTAGTCTTATTTACACAATAAAATGGTTAGGTACATTCGATGAGGAGAATGAAATATTTACCATAAATCCGAAGATATTGGAAAATTGTCAAATTGTCGTGGACAAATCAACTTTAGAATTTGAGACAATCTTTTTGCGGATCTTTCAGCTAATTAATCCATCGTTATTGTGTGGAATAAAGTTAGATAAGAATGTTATTATTGTACCGACATTGTCAATGAGAGCGAAAGGTTCATTTATCAGTTTTCACGATGAGTTTGAATACATTTCCGGTGGACAAACAGTTCGTGATTCTTTAAGTAATTTTGTTTTTTTCGATTATGGTGGCACATTTACAAAATCACAGTGTCCAGGCAAGAATGGAGGTAGTTTTACAACTTATCGTTTCTTAACTTCCAACATGATTGTTTATCGTCAAGAAATTAAAAAAAATAATTGCACTGTTATTCAAACTTATAATAAAAATATTTTTAACATTAATAATCGCAGAACAGCTAAGCTGTTTCTGTTTAGAAATAATCAATTTTTTAATATAGTCACCCAAAGTGAAATTTCATATGAATTGTATCAACCATTTTTTACCGTTTTTAAAATAGTATGGTCTAATTACAAACGCAGCGAAGATGGTCGCATGGAAATTTGGCCGACAATAGAAATAATTGCTGCGTGTGTTCAGCTAATT
>JABAIF010007716.1:423933-424160(-):Leptopilina_boulardi_GCA_015476485.1
ATTGGTTTTGTAAAATCGGAAATTCGTAAAGGTGTTCTACCTCGAATGTTGAGTGAAATTTTAAATACTCGTTTAATGGTGAAAAAAGCAATGAAGGATCATTCAAAGGATGATCGATCTTTACAACAGGTTCTTCATAATAGACAATTGGGATTAAAATTAATCGCTAATGTAACTTATGGTTATACATCGGCAAATTTCAGTGGAAGAATGTCATGTATTGAGGTA
#Usage example on cluster 

#nohup snakemake -j 8000  -s Snakemake_Alignment_main_Phylogeny  --cluster "sbatch -J {params.name} -p normal -N 1 --cpus-per-task  {params.threads}  --exclude='pbil-deb33' -o {params.out} -e {params.err}  " &> nohup_Alignment_main_Phylogeny_snakemake.out &

the_list=['p47_VcBV',
'p47_MdBV',
'CAR31573.1_CcBV',
'B7SV41_OrNV',
'AHW98279.1_P47_NlBV',
'ATY70215_DmNV_tom',
'p47_FaBV',
'p47-like_BtBV',
'A0A0B4VFM8_ToNV',
'A4L232_GbNV',
'A0A7D5UN34_DhNV',
'A0A076FC45_PmNV',
'AAN04369_HzNV-1',
'Q9YMS6_LdMNPV',
'P34051_AcMNPV',
'A0A097P1H3_CpV',
'Q91EY7_CpV']
  
import re 
from Bio import SeqIO 
import os 
#Define your paths :

Cluster_seqs_path="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final" #The path where all cluster seq not aligned will be present
Cluster_alignment_path="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_alignment" #The path where all aligned cluster seq will be present
Cluster_phylogeny_path="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_phylogeny"  #The path where all cluster phylognies seq will be present

mmseqs="/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs"


#If more than 2 taxa : 
def file_size(filePath):
    with open(filePath, 'r') as myfile:
	    data=myfile.read().replace('\n', ' ')
	    nb=data.count(">")
	    return nb

#For alignment with more than 1 taxa
SAMPLES=[]
for file in os.listdir(Cluster_seqs_path):
    taxa_list=[]
    if file.endswith("_phylogeny.aa"):
      taxa_list=[]
      for record in SeqIO.parse(Cluster_seqs_path+"/"+file, "fasta"):
         taxa_list.append(re.sub(".*_","",record.id))
      if len(list(dict.fromkeys(taxa_list)))>1:
         Clustername=re.sub("_phylogeny.aa","",file)
         SAMPLES.append(Clustername)	

print(SAMPLES)

#For phylogeny with more than 2 taxas 
SAMPLES2=[]
for file in os.listdir(Cluster_seqs_path):
    if file.endswith("_phylogeny.aa"):
       if file_size(Cluster_seqs_path+"/"+file)>2:
         Clustername=re.sub("_phylogeny.aa","",file)
         SAMPLES2.append(Clustername)	


rule all:
  input:
        expand(Cluster_alignment_path+"/{cluster_number}_main_phylogeny.aa.ali", cluster_number = SAMPLES),
	      expand(Cluster_alignment_path+"/{cluster_number}_main_phylogeny.aa.ali.trimmed", cluster_number = SAMPLES),
	      Cluster_phylogeny_path+"/Concatenated_sequences_main_phylogeny.treefile"


#Run Clustal analysis 
rule Clustal_cluster_alignment:
  log: Cluster_alignment_path+"/{cluster_number}.log"
  params:
     threads="3",
     time="2:00:00",
     name="Cluster_alignment_{cluster_number}",
     out=Cluster_alignment_path+"/Clustal_run_{cluster_number}.out",
     err=Cluster_alignment_path+"/Clustal_run_{cluster_number}.error"
  input:
     Cluster_file=Cluster_seqs_path+"/{cluster_number}_phylogeny.aa"
  output: 
     Alignment_cluster_file=Cluster_alignment_path+"/{cluster_number}_main_phylogeny.aa.ali",
  shell:
     """
     hostname
     cd {Cluster_seqs_path}
     /beegfs/data/bguinet/TOOLS/clustalo -i {input.Cluster_file} -o {output.Alignment_cluster_file} --threads {params.threads}
     """

#Run trimal analysis 
rule Trimal_cluster:
  log: Cluster_alignment_path+"/{cluster_number}.log"
  params:
     threads="1",
     time="10:00",
     name="Cluster_alignment_{cluster_number}",
     out=Cluster_alignment_path+"/Trimal_run_{cluster_number}.out",
     err=Cluster_alignment_path+"/Trimal_run_{cluster_number}.error"
  input:
     Alignment_cluster_file=Cluster_alignment_path+"/{cluster_number}_main_phylogeny.aa.ali"
  output:
     Alignment_cluster_file_trimmed=Cluster_alignment_path+"/{cluster_number}_main_phylogeny.aa.ali.trimmed"
  shell:
     """
     hostname
     cd {Cluster_seqs_path}
     /beegfs/data/bguinet/TOOLS/trimal/source/trimal -in {input.Alignment_cluster_file} -out {output.Alignment_cluster_file_trimmed} -automated1 -resoverlap 0.30 -seqoverlap 30 -fasta 
     """

#all dsDNA phylogeny 
rule dsDNA_phylogeny:
  params:
     threads="20",
     time="6:00:00",
     name="dsDNA_phylogeny",
     out=Cluster_phylogeny_path+"/Concatenated_sequences_main_phylogeny.out",
     err=Cluster_phylogeny_path+"/Concatenated_sequences_main_phylogeny.error"
  output:
     dsDNA_partition=Cluster_phylogeny_path+"/partition_main_phylogeny.txt",
     dsDNA_partition2=Cluster_phylogeny_path+"/partition_main_phylogeny.tab",
     Concatenated_sequence_alignment=Cluster_phylogeny_path+"/Concatenated_sequences_main_phylogeny.aln",
     Treefile_output=Cluster_phylogeny_path+"/Concatenated_sequences_main_phylogeny.treefile"
  shell:
     """
     hostname
     cd {Cluster_phylogeny_path}
     #Concatenate cluster trimmed alignment 
     perl /beegfs/home/bguinet/these_scripts_2/catfasta2phyml.pl -f --concatenate --verbose {Cluster_alignment_path}/*_main_phylogeny.aa.ali.trimmed > {output.Concatenated_sequence_alignment}  2> {output.dsDNA_partition} || true
     #Transform to partition.tab 
     grep '_main_phylogeny.aa.ali.trimmed =' {output.dsDNA_partition} >> {output.dsDNA_partition2} || true
     sed -i "s@Cluster@AA, Cluster@g" {output.dsDNA_partition2}
     cat {output.dsDNA_partition2}
     #Inferring species tree with ultraboostap  and alrt 
     /beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s {output.Concatenated_sequence_alignment} -spp {output.dsDNA_partition2} --prefix Concatenated_sequences_main_phylogeny -m MFP -alrt 5000  -bb 5000 -bnni --symtest-remove-bad  -nt {params.threads}
     #Inferring locus trees
     #/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s {output.Concatenated_sequence_alignment} -S concat.best_model.nex --prefix loci  -nt 10 
     #compute concordance factors
     #/beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -t {output.Treefile_output} --gcf loci.treefile -s {output.Concatenated_sequence_alignment}  --scf 100 --prefix concord  -nt 10
     """


#
Blast_table['strand'] = np.where(Blast_table['Names'].str.contains("\\+"), "+", "-")
Blast_table['start-end']=Blast_table['Names'].str.replace(".*:","")
Blast_table['start-end']=Blast_table['start-end'].str.replace("\\(.*","")
Blast_table['start']=Blast_table['start-end'].str.replace("-.*","")                  
Blast_table['end']=Blast_table['start-end'].str.replace(".*-","")


for group_name, df_group in Blast_table_grouped:
  # if there are viruses and more than two sequence within a cluster, create it !
  if df_group['Virus_count'].iloc[0]>0:
...       if len(df_group['Names'].unique())>=2:
...         #print("writting cluster : ", group_name)
...         # EVE protein loci 
...         # ALL sequences
...         with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_plot.aa","w") as output_aa:
...                 for row_index, row in df_group.iterrows():
...                  if row['Names2']in row['Names2'] in list_cynipoidea:
...                   if row['Names'] not in Added_sequences1:
...                     if row['Best_hit_ORF_perc']> 0.70:
...                       New_names= row['Scaffold_name']+":"+str(int(row['ORF_start']))+"-"+str(int(row['ORF_end']))+"("+str(row['ORF_strand'])+"):"+row['Species_name']+' [ORF]'
...                       if row['Mean_dNdS']+row['SE_dNdS']< 0.8 and row['Pvalue_dNdS']<0.05 :
...                           New_names=New_names+' [Purifying]'
...                       print('>',New_names,sep="",file=output_aa)
...                       print(AA_records_ORFs[row['ORF_target']].seq,file=output_aa)
...                       Added_sequences1.append(row['Names'])
...                     else:
...                       New_names= row['Scaffold_name']+":"+str(int(row['start']))+"-"+str(int(row['end']))+"("+str(row['strand'])+"):"+row['Species_name']
...                       if row['Mean_dNdS']+row['SE_dNdS']< 0.8 and  row['Pvalue_dNdS']<0.05 :
...                         New_names=New_names+' [Purifying]'
...                       else:
...                         New_names=New_names
...                       print('>',New_names,sep="",file=output_aa)
...                       print(AA_records_origin[row['Names']].seq,file=output_aa)
...                       Added_sequences1.append(row['Names'])
...                   else:
...                     print("")
...                  else:
...                     print('>',row['Names'],sep="",file=output_aa)
...                     print(AA_records_origin[row['Names']].seq,file=output_aa)
...                     Added_sequences1.append(row['Names'])














#To decide of the threshold we need to first run the /beegfs/data/bguinet/Cynipoidea_project/Synteny_simulation.py script 
/beegfs/data/bguinet/Cynipoidea_project/Synteny_analysis.py --main_file /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/dsDNA_hmmer_clusters_cov_NR_ORF_euk_repeat_scores.tab -o {output.Synteny_output}



>JABAIF010007716.1_423933-424160_-__Leptopilina_boulardi_GCA_015476485
IGFVKSEIRKGVLPRMLSEILNTRLMVKKAMKDHSKDDRSLQQVLHNRQLGLKLIANVTYGYTSANFSGRMSCIEV

['Cluster89_dNdS_LbFVorf94', 'Cluster66-Cluster338_dNdS_LbFVorf58-DNApol', 'Cluster183-Cluster1137-Cluster172-Cluster1324_dNdS_LbFVorf85-Ac81', 'Cluster66-Cluster338_dNdS_LbFVorf58-DNApol', 'Cluster198-Cluster2434_dNdS_LbFVorf108', 'Cluster2376-Cluster90_dNdS_LbFVorf96-lef-8', 'Cluster71_dNdS_LbFVorf13', 'Cluster20_dNdS_LbFVorf5', 'Cluster945_dNdS_LbFVorf10', 'Cluster52_dNdS_LbFVorf78-lef-9', 'Cluster60-Cluster21-Cluster259_dNdS_LbFVorf107-lef-4', 'Cluster60-Cluster21-Cluster259_dNdS_LbFVorf107-lef-4', 'Cluster2141-Cluster80-Cluster2173-Cluster1262-Cluster40_dNdS_LbFVorf68-helicase', 'Cluster2141-Cluster80-Cluster2173-Cluster1262-Cluster40_dNdS_LbFVorf68-helicase', 'Cluster183-Cluster1137-Cluster172-Cluster1324_dNdS_LbFVorf85-Ac81', 'Cluster222_dNdS_LbFVorf87', 'Cluster34_dNdS_LbFVorf72', 'Cluster71_dNdS_LbFVorf13', 'Cluster89_dNdS_LbFVorf94', 'Cluster2376-Cluster90_dNdS_LbFVorf96-lef-8', 'Cluster945_dNdS_LbFVorf10', 'Cluster843-Cluster250_dNdS_LbFVorf83', 'Cluster3_dNdS_LbFVorf60-lcat', 'Cluster3_dNdS_LbFVorf60-lcat', 'Cluster843-Cluster250_dNdS_LbFVorf83', 'Cluster34_dNdS_LbFVorf72', 'Cluster222_dNdS_LbFVorf87', 'Cluster130-Cluster2564_dNdS_LbFVorf92', 'Cluster198-Cluster2434_dNdS_LbFVorf108', 'Cluster20_dNdS_LbFVorf5', 'Cluster52_dNdS_LbFVorf78-lef-9'

#Usage example on cluster 

#nohup snakemake -j 8000  -s Snakemake_alignment_and_phylogenies_step4   --cluster "sbatch -J {params.name} -p normal -N 1 --cpus-per-task  {params.threads}  --exclude='pbil-deb33' -o {params.out}  " &> nohup_Snakemake_alignment_and_phylogenies_step4.out &


import re 
import os 
from Bio import SeqIO
#Define your paths :
Cluster_seq_path="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/" #The path where all cluster seq not aligned will be present
Cluster_alignment_path="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_alignment/" #The path where all aligned cluster seq will be present
Cluster_phylogeny_path="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_phylogeny/"  #The path where all cluster phylognies seq will be present


#For alignment with more than Nb_taxa_in_phylo taxa
SAMPLES=[]
for file in os.listdir(Cluster_seq_path):
    if "_plot.aa" in file:
      Clustername=re.sub("_plot.aa","",file)
      SAMPLES.append(Clustername)


#SAMPLES=['Cluster1677-Cluster463-Cluster352-Cluster1676-Cluster648_LbFVorf108']      
rule all:
  input:
        expand(Cluster_alignment_path+"{cluster_number}_plot.aa.ali", cluster_number = SAMPLES),
        expand(Cluster_phylogeny_path+"{cluster_number}_plot.aa.ali.treefile", cluster_number = SAMPLES),



#Run Clustal analysis 
rule Clustal_cluster_alignment:
  log: Cluster_alignment_path+"/{cluster_number}.log"
  params:
     threads="3",
     time="2:00:00",
     name="Cluster_alignment_{cluster_number}",
     out=Cluster_alignment_path+"/Clustal_run_{cluster_number}.out",
     err=Cluster_alignment_path+"/Clustal_run_{cluster_number}.error"
  input:
     Cluster_file=Cluster_seqs_path+"/{cluster_number}_plot.aa"
  output:
     Alignment_cluster_file=Cluster_alignment_path+"/{cluster_number}_plot.aa.ali",
  shell:
     """
     hostname
     cd {Cluster_seqs_path}
     /beegfs/data/bguinet/TOOLS/clustalo -i {input.Cluster_file} -o {output.Alignment_cluster_file} --threads {params.threads}
     """

#Run IQTREE analysis 
rule Cluster_phylogeny:
  params:
    threads = "5",
    time = "5:00:00",
    name = "Cluster_phylogeny_plot_{cluster_number}",
    out = Cluster_phylogeny_path+"IQTREE_run_plot_{cluster_number}.out",
    err = Cluster_phylogeny_path+"IQTREE_run_plot_{cluster_number}.error"
  input:
    Cluster_alignment_file=Cluster_alignment_path+"{cluster_number}_plot.aa.ali"
  output:
    Cluster_alignment_output=Cluster_phylogeny_path+"{cluster_number}_plot.aa.ali.treefile"
  shell:
    """
    sed -i 's@!@-@g' {input.Cluster_alignment_file}
    /beegfs/data/bguinet/TOOLS/iqtree-2.1.2-Linux/bin/iqtree2 -s {input.Cluster_alignment_file} -m MFP -alrt 5000 -bb 5000  -nt {params.threads} -redo
    sed -i 's@____@_+__@g' {Cluster_alignment_path}{wildcards.cluster_number2}_plot.aa.ali.treefile
    mv {Cluster_alignment_path}{wildcards.cluster_number2}_plot.aa.* {Cluster_phylogeny_path}
    """


grouped = Blast_table.groupby('Cluster_hmmer')
for group_name, group in grouped:
  if len(group.loc[group['Cynip_count'].gt(0) & group['Filamentous_count'].gt(0)])>0:
      print('Name : ',group['Prot_name'].iloc[0])
      try:
        cluster=group.loc[group['Names2'].isin(['RhEFV','ThEFV','TrEFV','LbEFV','LbEFV','LcEFV','LhEFV','CcFV1','CcFV2','EfFV','PoFV','PcFV','LhFV']) & group['Best_hit_ORF_perc'].lt(0.70)][['Names','Cluster_hmmer','Names2','Prot_name','Pseudogenized','Mean_dNdS','Pvalue_dNdS','Best_hit_ORF_perc']]['Cluster_hmmer'].iloc[0]
        Blast_table.loc[Blast_table['Cluster_hmmer'].str.contains(cluster,na=False) & group['Names2'].isin(['RhEFV','ThEFV','TrEFV','LbEFV','LbEFV','LcEFV','LhEFV','CcFV1','CcFV2','EfFV','PoFV','PcFV','LhFV'])][['Names','Cluster_hmmer','Names2','Prot_name','Pseudogenized','Mean_dNdS','Pvalue_dNdS','Best_hit_ORF_perc']]
        print("\n")
      except:
        continue
###
#Run Event analysis 
rule Asses_events_along_phylogeny:
  params:
    threads = "1",
    time = "15:00",
    name = "Cynipoidea_monophyletique_assessment",
    out = Monophyletic_assessment_path+"Cynipoidea_monophyletique_assessment.out",
    err = Monophyletic_assessment_path+"Cynipoidea_monophyletique_assessment.error"
  input:
    treefile = Monophyletic_assessment_path+"Cynipoidea_phylogeny.nwk",
    Cluster_directory = Cluster_phylogeny_path,
    Blast_table_scores = Clustering_path+"dsDNA_hmmer_clusters_cov_NR_ORF_euk_repeat_scores_synteny.tab",
    taxo_table = Monophyletic_assessment_path+"Species_sub_sup_families_order.txt",
  output:
    Output_table_file = Monophyletic_assessment_path+"Monophyletic_table.tab",
    output_monop = Clustering_path+"dsDNA_hmmer_clusters_cov_NR_ORF_euk_repeat_scores_synteny_event.tab"
  shell:
    """
    /beegfs/data/soft/R-4.0.5/bin/Rscript /beegfs/data/bguinet/Cynipoidea_project/Monophyletic_assessment.R {input.treefile} {input.Cluster_directory} 0.8 {input.Blast_table_scores} {input.taxo_table} {output.Output_table_file}
    python3 /beegfs/data/bguinet/Cynipoidea_project/Add_Monophyletic_table.py -M {output.Output_table_file} -b {input.Blast_table_scores} -o {output.output_monop}
    """


































#



import pandas as pd 
import sys 
import argparse
import os
from collections import defaultdict
import pandas as pd
from Bio import SearchIO
import re 
from Bio import SeqIO
import subprocess
import numpy as np 

# Print out a message when the program is initiated.
print('------------------------------------------------------------------------------\n')
print('                        Create cluster files                                   .\n')
print('------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Gather HMMEr files')
parser.add_argument("-b", "--blast_cluster_file", help="The blast cluster file")

blast_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event_dNdS.tab"

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}



#print ("\n")

print("Writting all the cluster files...")

Blast_table=pd.read_csv(blast_file,sep=";")

Blast_table=Blast_table.loc[~Blast_table['Species_name'].isin(['Leptopilina_boulardi_GCA_015476485.1','Leptopilina_boulardi_GCA_003121605.1',
       'Leptopilina_boulardi_GCA_011634795.1','Leptopilina_heterotoma_GCA_009602685.1',
       'Leptopilina_heterotoma_GCA_010016045.1',
       'Leptopilina_heterotoma_GCA_009026005.1',
       'Leptopilina_heterotoma_GCA_009025955.1','Leptopilina_boulardi_GCA_015476485','Leptopilina_boulardi_GCA_003121605',
              'Leptopilina_boulardi_GCA_011634795','Leptopilina_heterotoma_GCA_009602685',
              'Leptopilina_heterotoma_GCA_010016045',
              'Leptopilina_heterotoma_GCA_009026005',
              'Leptopilina_heterotoma_GCA_009025955'])]

Blast_table.loc[Blast_table['Names'].str.contains("MdBV"),"Names2"]="MdENV"
Blast_table.loc[Blast_table['Names'].str.contains("CcBV"),"Names2"]="CcENV"
#############################################
# Write the AA clusters                    ##
#############################################

print("Writting AA cluster files...")
# iterate over each group to get EVE loci sequences 

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}
    
AA_records_origin=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses.faa","fasta"))
AA_records_ORFs=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.aa","fasta"))


list_cynipoidea = ['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV']
list_filamentous = ['LbFV','EfFV','PoFV','PcFV','CcFV1','CcFV2','LhFV']

    

Blast_table=Blast_table.loc[Blast_table['Cynip_count'].gt(0) &  Blast_table['Filamentous_count'].gt(0)]

Blast_table['strand'] = np.where(Blast_table['Names'].str.contains("\\+"), "+", "-")
Blast_table['start-end']=Blast_table['Names'].str.replace("\\(.*","")
Blast_table['start-end']=Blast_table['start-end'].str.replace(".*:","")
Blast_table['start']=Blast_table['start-end'].str.replace("-.*","")                  
Blast_table['end']=Blast_table['start-end'].str.replace(".*-","")

#Blast_table=Blast_table.loc[Blast_table['Cluster_hmmer'].str.contains("279")]

Added_sequences1=[]
Added_sequences2=[]
Blast_table_grouped = Blast_table.groupby('Cluster_hmmer')
for group_name, df_group in Blast_table_grouped:
  # if there are viruses and more than two sequence within a cluster, create it !
  if df_group['Virus_count'].iloc[0]>0:
          if len(df_group['Names'].unique())>=2:
            #print("writting cluster : ", group_name)
            # EVE protein loci 
            # ALL sequences
            with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_plot.aa","w") as output_aa:
                    for row_index, row in df_group.iterrows():
                     if row['Names2']in row['Names2'] in list_cynipoidea:
                       if row['Names'] not in Added_sequences1:
                        if row['Best_hit_ORF_perc']> 0.70:
                          New_names= row['Scaffold_name']+":"+str(int(row['ORF_start']))+"-"+str(int(row['ORF_end']))+"("+str(row['ORF_strand'])+"):"+row['Species_name']+' [ORF]'
                          if row['Mean_dNdS']+row['SE_dNdS']< 0.8 and row['Pvalue_dNdS']<0.05 :
                             New_names=New_names+' [Purifying]'
                          print('>',New_names,sep="",file=output_aa)
                          print(AA_records_ORFs[row['ORF_target']].seq,file=output_aa)
                          Added_sequences1.append(row['Names'])
                        else:
                          New_names= row['Scaffold_name']+":"+str(int(row['start']))+"-"+str(int(row['end']))+"("+str(row['strand'])+"):"+row['Species_name']
                          if row['Mean_dNdS']+row['SE_dNdS']< 0.8 and  row['Pvalue_dNdS']<0.05 :
                            New_names=New_names+' [Purifying]'
                          else:
                            New_names=New_names
                          print('>',New_names,sep="",file=output_aa)
                          print(AA_records_origin[row['Names']].seq,file=output_aa)
                          Added_sequences1.append(row['Names'])
                     else:
                        print('>',row['Names'],sep="",file=output_aa)
                        print(AA_records_origin[row['Names']].seq,file=output_aa)
                        Added_sequences1.append(row['Names'])

print ("ALL AA files with all sequences written to : /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/ with extension plot.aa")



row['Pvalue_dNdS']=row['Pvalue_dNdS'].astype(int)
for group_name, df_group in Blast_table_grouped:
    # if there are viruses and more than two sequence within a cluster, create it !
    if df_group['Virus_count'].iloc[0]>0:
      if len(df_group['Names'].unique())>=2:
        #print("writting cluster : ", group_name)
        # EVE protein loci 
        # ALL sequences
        with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_plot.aa","w") as output_aa:
                for row_index, row in df_group.iterrows():
                  if row['Names'] not in Added_sequences1:
                    if row['Best_hit_ORF_perc']> 0.70:
                      New_names= row['Scaffold_name']+":"+str(int(row['ORF_start']))+"-"+str(int(row['ORF_end']))+"("+str(row['ORF_strand'])+"):"+row['Species_name']+' [ORF]'
                      if row['Mean_dNdS']+row['SE_dNdS']< 0.8 and row['Pvalue_dNdS']<0.05 :
                          New_names=New_names+' [Purifying]'
                      print('>',New_names,sep="",file=output_aa)
                      print(AA_records_ORFs[row['ORF_target']].seq,file=output_aa)
                      Added_sequences1.append(row['Names'])
                    else:
                      if row['Mean_dNdS']+row['SE_dNdS']< 0.8 and  row['Pvalue_dNdS']<0.05 :
                        New_names=row['Names']+' [Purifying]'
                      else:
                        New_names=row['Names']
                      print('>',New_names,sep="",file=output_aa)
                      print(AA_records_origin[row['Names']].seq,file=output_aa)
                      Added_sequences1.append(row['Names'])
                  else:
                    print("")
##
print ("ALL AA files with all sequences written to : /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/ with extension .aa")




#SBATCH --mem 10G
#SBATCH -t 5:00:00
#SBATCH --cpus-per-task=1
#SBATCH -e /beegfs/data/bguinet/Cynipoidea_project/Scripts/Synteny_job.error
#SBATCH -o /beegfs/data/bguinet/Cynipoidea_project/Scripts/Synteny_job.out
#SBATCH -J Synteny job

python3  /beegfs/data/bguinet/Cynipoidea_project/Scripts/Synteny.py
    
# Create clusters 
import pandas as pd 
import re 
from Bio import SeqIO
import argparse
import os 
import subprocess
import numpy as np

# Print out a message when the program is initiated.
print('--------------------------------------------------------------------------------------------------\n')
print('              Run blast synteny analysis                                                         .\n')
print('--------------------------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Run blast synteny analysis')
parser.add_argument("-main", "--main_file", help="The name of the main file with all the other informations")
parser.add_argument("-o", "--output_file", help="The output desired file name ")
args = parser.parse_args()

# example usage
"""
python3 Add_metaeuk_repeat_and_scores.py 
 --main_file /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/dsDNA_hmmer_clusters_cov_NR_ORF_euk_repeat_scores.tab \
 -o /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/dsDNA_hmmer_clusters_cov_NR_ORF_euk_repeat_scores_synteny.tab

""" 

Output_file=args.output_file
Main_file=args.main_file

Output_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Syntenie.tab"
Main_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score.tab"
Main_table=pd.read_csv(Main_file,sep=";")

# keep only scaffolds of interest

Main_table=Main_table.loc[Main_table['Scaffold_score'].isin(["A","B","C","D"])]
Main_table=Main_table.loc[Main_table['Names2'].isin(['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV'])]


if os.path.exists("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result_reduced.m8"):
  Synteny_table=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result_reduced.m8",sep=";")
else:

  # run blast all vs all on scaffold containing EVEs 

  # Create database with all candidate scaffolds 

  def to_dict_remove_dups(sequences):
      return {record.id: record for record in sequences}
      
  record_dict = to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/All_candidate_scaffolds.dna","fasta"))

  with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred.dna","w") as output :
    for loci in Main_table['Scaffold_and_Species_name'].unique():
      for sp in ['Leptopilina_heterotoma_GCA_015476425.1','Leptopilina_boulardi_GCA_019393585.1', 'Trybliographa','Rhoptromeris', 'Leptopilina_clavipes', 'Thrichoplasta']:
        if sp in loci:
          print(loci)
          print(">",record_dict[loci].id,sep="",file=output)
          print(record_dict[loci].seq,file=output)
    #Add also LbFV and LhFV genomes 
    for loci in SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Genomes/LbFV/LbFV_free.dna","fasta"):
      print(">LbFV:",loci.id,sep="",file=output)
      print(loci.seq,file=output)
    for loci in SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Genomes/LhFV/LhFV_free.dna","fasta"):
      print(">LhFV:",loci.id,sep="",file=output)
      print(loci.seq,file=output)
    #Add also EVEs in dna format 
    for loci in SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses.dna","fasta"):
      for sp in ['Leptopilina_heterotoma_GCA_015476425.1','Leptopilina_boulardi_GCA_019393585.1', 'Trybliographa','Rhoptromeris', 'Leptopilina_clavipes', 'Thrichoplasta']:
        if sp in loci.id:
          #print(">",record_dict[loci].id,sep="",file=output)
          #print(record_dict[loci].seq,file=output)
          print(">",loci.id,sep="",file=output)
          print(loci.seq,file=output)
    

  if os.path.exists("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result.m8"):
    print("Blast file already exist")
  else:
    subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred.dna /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_db", shell=True)
    subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs search /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_db /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_db /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_tpm --threads 20 --search-type 4", shell=True)
    subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs convertalis --search-type 4  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tlen,qlen,qcov,tcov' /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_db /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_db /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result.m8" ,shell=True)

  # Open blast synteny 
  Synteny_table=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result.m8",sep="\t",header=None)
  Synteny_table.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']

  Synteny_table=Synteny_table.loc[~(Synteny_table['query'].str.contains("_orf") | Synteny_table['target'].str.contains("_orf"))]
  Synteny_table=Synteny_table.loc[~(Synteny_table['query'].str.contains("AQQ") | Synteny_table['target'].str.contains("AQQ"))]

  #Generate the comparison format table from GenoplotR 

  #Synteny_table.loc[Synteny_table['query'].str.contains("LhFV") & Synteny_table['target'].str.contains("LbFV") | Synteny_table['query'].str.contains("LbFV") & Synteny_table['target'].str.contains("LhFV") ]

  #Keep only convincing hits 
  Synteny_table=Synteny_table.loc[Synteny_table['bits'].ge(50)]
  #Synteny_table=Synteny_table.loc[Synteny_table['evalue'].le(0.000001)]

  Synteny_table['query_SP']="NA"
  Synteny_table['query_SP'].mask(Synteny_table["query"].str.contains("LhFV"), 'LhFV', inplace=True)
  Synteny_table['query_SP'].mask(Synteny_table["query"].str.contains("LbFV"), 'LbFV', inplace=True)
  Synteny_table['query_SP'].mask(Synteny_table["query"].str.contains("Thrichoplasta"), 'ThEFV', inplace=True)
  Synteny_table['query_SP'].mask(Synteny_table["query"].str.contains("Trybliographa"), 'TrEFV', inplace=True)
  Synteny_table['query_SP'].mask(Synteny_table["query"].str.contains("Rhoptromeris"), 'RhEFV', inplace=True)
  Synteny_table['query_SP'].mask(Synteny_table["query"].str.contains("clavipes"), 'LcEFV', inplace=True)
  Synteny_table['query_SP'].mask(Synteny_table["query"].str.contains("boulardi"), 'LbEFV', inplace=True)
  Synteny_table['query_SP'].mask(Synteny_table["query"].str.contains("heterotoma"), 'LhEFV', inplace=True)
  Synteny_table['target_SP']="NA"
  Synteny_table['target_SP'].mask(Synteny_table["target"].str.contains("LhFV"), 'LhFV', inplace=True)
  Synteny_table['target_SP'].mask(Synteny_table["target"].str.contains("LbFV"), 'LbFV', inplace=True)
  Synteny_table['target_SP'].mask(Synteny_table["target"].str.contains("Thrichoplasta"), 'ThEFV', inplace=True)
  Synteny_table['target_SP'].mask(Synteny_table["target"].str.contains("Trybliographa"), 'TrEFV', inplace=True)
  Synteny_table['target_SP'].mask(Synteny_table["target"].str.contains("Rhoptromeris"), 'RhEFV', inplace=True)
  Synteny_table['target_SP'].mask(Synteny_table["target"].str.contains("clavipes"), 'LcEFV', inplace=True)
  Synteny_table['target_SP'].mask(Synteny_table["target"].str.contains("boulardi"), 'LbEFV', inplace=True)
  Synteny_table['target_SP'].mask(Synteny_table["target"].str.contains("heterotoma"), 'LhEFV', inplace=True)

  Synteny_table=Synteny_table.loc[~(Synteny_table['query_SP'].eq("NA") | Synteny_table['target_SP'].eq("NA"))]
  Synteny_table=Synteny_table.loc[~(Synteny_table['query_SP'] == Synteny_table['target_SP'])]

  """
  Synteny_table[['query_SP', 'query']] = Synteny_table['query'].str.split(':', 1, expand=True)
  Synteny_table[['target_SP', 'target']] = Synteny_table['target'].str.split(':', 1, expand=True)
  Synteny_table=Synteny_table.loc[~(Synteny_table['query_SP'] == Synteny_table['target_SP'])]
  Synteny_table['query'] = Synteny_table['query']+':'+Synteny_table['query_SP']
  Synteny_table['target'] = Synteny_table['target']+':'+Synteny_table['target_SP']
  """

  #Remove overlapping hits and keep the best one 


  #Change strand 

  import numpy as np
  Synteny_table['qstrand']=np.where(Synteny_table["qstart"]>Synteny_table["qend"],'-','+')
  m = Synteny_table['qstrand'].eq('-')
  Synteny_table.loc[m, ['qstart','qend']] = Synteny_table.loc[m, ['qend','qstart']].to_numpy()

  Synteny_table['tstrand']=np.where(Synteny_table["tstart"]>Synteny_table["tend"],'-','+')
  m = Synteny_table['qstrand'].eq('-')
  Synteny_table.loc[m, ['tstart','tend']] = Synteny_table.loc[m, ['tend','tstart']].to_numpy()

  is_overlapped = lambda x: x['qstart'] >= x['qend'].shift(fill_value=-1)
  Synteny_table['group'] = Synteny_table.sort_values(['query', 'qstart', 'qend']) \
      .groupby(['query','target'],as_index=False).apply(is_overlapped).droplevel(0).cumsum()

  Synteny_table = Synteny_table.sort_values(['group', 'evalue', 'bits'], ascending=[True, True, False]) \
      .groupby(Synteny_table['group']).head(1)

  is_overlapped = lambda x: x['tstart'] >= x['tend'].shift(fill_value=-1)
  Synteny_table['group2'] = Synteny_table.sort_values(['target', 'tstart', 'tend']) \
      .groupby(['target','query'],as_index=False).apply(is_overlapped).droplevel(0).cumsum()

  Synteny_table = Synteny_table.sort_values(['group2', 'evalue', 'bits'], ascending=[True, True, False]) \
      .groupby(Synteny_table['group2']).head(1)

  Synteny_table.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result_reduced.m8",sep=";",index=False)
  #

subtab=Main_table # <- voir avec l'ancien
#subtab['query2'] = subtab['query'].str.replace(r'\(.*', '')
#subtab['query2'] = subtab['query2'].str.replace(r'.*:', '')
#subtab[['start', 'end']] = subtab['query2'].str.split('-', 1, expand=True)

subtab=subtab.drop_duplicates(subset = 'Names', keep = 'first')

subtab_copy= subtab.copy()
subtab_copy=subtab_copy[['Names','start','end','Scaffold_name','Species_name']]
subtab_copy['query'] = subtab_copy['Species_name']+':'+subtab_copy['Scaffold_name']

#First we will remove every overlaping coordinates with EVEs 

Synteny_table=Synteny_table.loc[~Synteny_table['query'].str.contains("\\)")]
Synteny_table=Synteny_table.loc[~Synteny_table['target'].str.contains("\\)")]
Synteny_table1=Synteny_table[['query','qstart','qend']]
Synteny_table1.columns=['query','start','end']
subtab_copy['key'] = subtab_copy['query'].str.upper()
Synteny_table1['key']=Synteny_table1['query'].str.upper()

Synteny_table_subtab_copy1 = Synteny_table1.reset_index().merge(subtab_copy, on='key', how='left', suffixes=['', '_y'])
Synteny_table_subtab_copy1['start'] = Synteny_table_subtab_copy1['start'].astype(float)
Synteny_table_subtab_copy1['end'] = Synteny_table_subtab_copy1['end'].astype(float)
Synteny_table_subtab_copy1['start_y'] = Synteny_table_subtab_copy1['start_y'].astype(float)
Synteny_table_subtab_copy1['end_y'] = Synteny_table_subtab_copy1['end_y'].astype(float)

Synteny_table_subtab_copy1_1=Synteny_table_subtab_copy1.loc[(Synteny_table_subtab_copy1['start_y']>= Synteny_table_subtab_copy1['start'])  & (Synteny_table_subtab_copy1['end']>= Synteny_table_subtab_copy1['start_y'])]
Synteny_table_subtab_copy1_2=Synteny_table_subtab_copy1.loc[(Synteny_table_subtab_copy1['start']>= Synteny_table_subtab_copy1['start_y'])  & (Synteny_table_subtab_copy1['end_y']>= Synteny_table_subtab_copy1['start'])]
Synteny_table_subtab_copy1_3=Synteny_table_subtab_copy1.loc[(Synteny_table_subtab_copy1['start_y']>= Synteny_table_subtab_copy1['start'])  & (Synteny_table_subtab_copy1['end_y']<= Synteny_table_subtab_copy1['end'])]
Synteny_table_subtab_copy1_4=Synteny_table_subtab_copy1.loc[(Synteny_table_subtab_copy1['start']>= Synteny_table_subtab_copy1['start_y'])  & (Synteny_table_subtab_copy1['end']<= Synteny_table_subtab_copy1['end_y'])]

Synteny_table_subtab_copy1=pd.concat([Synteny_table_subtab_copy1_1,Synteny_table_subtab_copy1_2,Synteny_table_subtab_copy1_3,Synteny_table_subtab_copy1_4])
Synteny_table_subtab_copy1['full_name']=Synteny_table_subtab_copy1['query']+"-"+Synteny_table_subtab_copy1['start'].astype(int).astype(str)+"-"+Synteny_table_subtab_copy1['end'].astype(int).astype(str)
Synteny_table_subtab_copy1['full_name']=Synteny_table_subtab_copy1['full_name'].str.upper()


Synteny_table2=Synteny_table[['query','tstart','tend']]
Synteny_table2.columns=['query','start','end']
subtab_copy['key'] = subtab_copy['query'].str.upper()
Synteny_table2=Synteny_table[['query','tstart','tend']]
Synteny_table2.columns=['query','start','end']
subtab_copy['key'] = subtab_copy['query'].str.upper()
Synteny_table2['key']=Synteny_table2['query'].str.upper()

Synteny_table_subtab_copy2 = Synteny_table2.reset_index().merge(subtab_copy, on='key', how='left', suffixes=['', '_y'])
Synteny_table_subtab_copy2['start'] = Synteny_table_subtab_copy2['start'].astype(float)
Synteny_table_subtab_copy2['end'] = Synteny_table_subtab_copy2['end'].astype(float)
Synteny_table_subtab_copy2['start_y'] = Synteny_table_subtab_copy2['start_y'].astype(float)
Synteny_table_subtab_copy2['end_y'] = Synteny_table_subtab_copy2['end_y'].astype(float)

Synteny_table_subtab_copy2_1=Synteny_table_subtab_copy2.loc[(Synteny_table_subtab_copy2['start_y']>= Synteny_table_subtab_copy2['start'])  & (Synteny_table_subtab_copy2['end']>= Synteny_table_subtab_copy2['start_y'])]
Synteny_table_subtab_copy2_2=Synteny_table_subtab_copy2.loc[(Synteny_table_subtab_copy2['start']>= Synteny_table_subtab_copy2['start_y'])  & (Synteny_table_subtab_copy2['end_y']>= Synteny_table_subtab_copy2['start'])]
Synteny_table_subtab_copy2_3=Synteny_table_subtab_copy2.loc[(Synteny_table_subtab_copy2['start_y']>= Synteny_table_subtab_copy2['start'])  & (Synteny_table_subtab_copy2['end_y']<= Synteny_table_subtab_copy2['end'])]
Synteny_table_subtab_copy2_4=Synteny_table_subtab_copy2.loc[(Synteny_table_subtab_copy2['start']>= Synteny_table_subtab_copy2['start_y'])  & (Synteny_table_subtab_copy2['end']<= Synteny_table_subtab_copy2['end_y'])]

Synteny_table_subtab_copy2=pd.concat([Synteny_table_subtab_copy2_1,Synteny_table_subtab_copy2_2,Synteny_table_subtab_copy2_3,Synteny_table_subtab_copy2_4])
Synteny_table_subtab_copy2['full_name']=Synteny_table_subtab_copy2['query']+"-"+Synteny_table_subtab_copy2['start'].astype(str)+"-"+Synteny_table_subtab_copy2['end'].astype(str)


Synteny_table['full_name_query']=Synteny_table['query']+'-'+Synteny_table['qstart'].astype(str)+'-'+Synteny_table['qend'].astype(str)
Synteny_table['full_name_target']=Synteny_table['target']+'-'+Synteny_table['tstart'].astype(str)+'-'+Synteny_table['tend'].astype(str)
Synteny_table['full_name_query']=Synteny_table['full_name_query'].str.upper()
Synteny_table['full_name_target']=Synteny_table['full_name_target'].str.upper()

Synteny_table2=Synteny_table.loc[~Synteny_table['full_name_query'].isin(Synteny_table_subtab_copy2['full_name'])] 
Synteny_table2=Synteny_table2.loc[~Synteny_table2['full_name_target'].isin(Synteny_table_subtab_copy2['full_name'])] 
Synteny_table2=Synteny_table2.loc[~Synteny_table2['full_name_query'].isin(Synteny_table_subtab_copy1['full_name'])] 
Synteny_table2=Synteny_table2.loc[~Synteny_table2['full_name_target'].isin(Synteny_table_subtab_copy1['full_name'])] 
##

count_success=0
count_success_high=0



def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences} 


              
if os.path.exists("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds_length.tab") :
  Length_table = pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds_length.tab",sep=";")
else:
  Record_dict=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.dna","fasta"))
  Length_table = pd.DataFrame(columns=("Scaffold_species_name","Scaffold_length"))
  for index, row in subtab.iterrows():
    Length_table=Length_table.append({"Scaffold_species_name":Record_dict[re.sub(".*:","",row['Names'])+":"+row['Scaffold_name']].id,"Scaffold_length":len(Record_dict[re.sub(".*:","",row['Names'])+":"+row['Scaffold_name']].seq)}, ignore_index=True)
  Length_table = Length_table.drop_duplicates(subset = ["Scaffold_species_name"])
  Length_table.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds_length.tab",sep=";",index=False)


subtab['Scaffold_species_name']=subtab['Names'].str.replace(".*:","")+":"+subtab['Scaffold_name']
subtab=subtab.merge(Length_table,on="Scaffold_species_name")
subtab=subtab.loc[subtab['Cynip_count'].gt(0) &  subtab['Filamentous_count'].gt(0)]


grouped = subtab.groupby(['Cluster_hmmer'])
#grouped = subsubtab.groupby(['Cluster_hmmer'])
  

from  statistics import * 
i=0


Synteny_table2_2=Synteny_table2[['target', 'query', 'pident', 'alnlen', 'mismatch', 'gapopen', 'tstart', 'tend', 'qstart', 'qend', 'evalue', 'bits', 'qlen', 'tlen', 'tcov', 'qcov', 'target_SP', 'query_SP', 'tstrand', 'qstrand', 'group', 'group2', 'full_name_target', 'full_name_query']]
Synteny_table2_2.columns=['query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tlen', 'qlen', 'qcov', 'tcov', 'query_SP', 'target_SP', 'qstrand', 'tstrand', 'group', 'group2', 'full_name_query', 'full_name_target']
Synteny_table2=pd.concat([Synteny_table2,Synteny_table2_2])
Synteny_table2 = Synteny_table2.drop_duplicates()
Synteny_table2 = Synteny_table2[~Synteny_table2[['full_name_query', 'full_name_target']].agg(frozenset, axis=1).duplicated()]

Result_table = pd.DataFrame(columns=("Windows","Cluster_hmmer","Mean_homologous","Median_homologous","Tot_alnlen","Mean_dNdS"))
subResult_table = pd.DataFrame(columns=("Windows","Cluster_hmmer","Scaffold_name","Species","Nb_HSPs","Tot_alnlen","PASS"))


for group_name, group in grouped:
        print(group_name)
        Cluster_hmmer=group['Cluster_hmmer'].iloc[0]
        Prot_name=group['Prot_name'].iloc[0]
        len_homologous=[]
        #print(group)
        #For each  scaffold containing EVE within a Cluster, look at the other scaffolds containing also EVEs within that cluster 
        for index1, row1 in group.iterrows():
            query=row1['query']
            EVE_start1 = int(row1['start'])
            EVE_end1 = int(row1['end'])
            EVE_scaffold1 = row1['Scaffold_name']
            EVE_species1= row1['Species_name']
            subsubtab=group.loc[~(group['query'].str.contains(EVE_species1,na=False))]
            #print(subtab)
            for index2, row2 in subsubtab.iterrows():
                target=row2['query']
                EVE_start2 = int(row2['start'])
                EVE_end2 = int(row2['end'])
                EVE_scaffold2 = row2['Scaffold_name']
                EVE_species2= row2['Species_name']
                # Look for HSPs in common between the scaffolds containing EVEs 
                subSynteny_table=Synteny_table2.loc[(Synteny_table2['query']==  EVE_species1+":"+EVE_scaffold1) & (Synteny_table2['target']== EVE_species2+":"+EVE_scaffold2)]
                #subSynteny_table[['qstart','qend','tstart','tend','query','target']]
                #print(EVE_scaffold1,":",EVE_species1," vs ",EVE_scaffold2,":",EVE_species2)
                #print(query, " - ", target)
                #Define the windows threshold size based on the min lenght of one of the two scaffolds 
                if min(subtab.loc[subtab['Scaffold_name']==EVE_scaffold2]['Scaffold_length'].iloc[0],subtab.loc[subtab['Scaffold_name']==EVE_scaffold1]['Scaffold_length'].iloc[0]) <= 10000 :
                    windows=10000
                elif min(subtab.loc[subtab['Scaffold_name']==EVE_scaffold2]['Scaffold_length'].iloc[0],subtab.loc[subtab['Scaffold_name']==EVE_scaffold1]['Scaffold_length'].iloc[0])  < 100000:
                    windows=100000
                elif min(subtab.loc[subtab['Scaffold_name']==EVE_scaffold2]['Scaffold_length'].iloc[0],subtab.loc[subtab['Scaffold_name']==EVE_scaffold1]['Scaffold_length'].iloc[0])  < 1000000:
                    windows=1000000
                elif min(subtab.loc[subtab['Scaffold_name']==EVE_scaffold2]['Scaffold_length'].iloc[0],subtab.loc[subtab['Scaffold_name']==EVE_scaffold1]['Scaffold_length'].iloc[0])  < 10000000:
                    windows=10000000
                else :
                    windows=100000000
                # Focus on the concerned scaffolds
                subSynteny_table=subSynteny_table.loc[~(subSynteny_table['qstart'].le(EVE_start1) & subSynteny_table['qend'].ge(EVE_start1) | subSynteny_table['qstart'].le(EVE_end1) & subSynteny_table['qend'].ge(EVE_end1))]
                # Focus only on HSP of the size of the defined windows around the EVE 
                subSynteny_table=subSynteny_table.loc[ subSynteny_table['qstart'].ge(EVE_start1-windows) & subSynteny_table['qend'].le(EVE_end1+windows) & subSynteny_table['tstart'].ge(EVE_start2-windows) & subSynteny_table['tend'].le(EVE_end2+windows) ][['query','target','qstart','qend','tstart','tend','qstrand','tstrand','alnlen','bits','evalue']]
                subSynteny_table
                #print(subSynteny_table)
                #print("\n")
                if len(subSynteny_table) >0 :
                                    len_homologous.append(len(subSynteny_table))
                else: 
                                    len_homologous.append(0)
                if len(subSynteny_table) >0 :
                                    count_success+=1
                if len(subSynteny_table) >10 :
                                    count_success_high+=1 
                if windows == 10000:
                    threshold=180
                elif windows == 100000:
                    threshold = 1545 # 1332
                elif windows == 1000000:
                    threshold = 49159 # 21000
                elif windows == 10000000:
                    threshold = 176179 # 74262
                elif windows == 100000000:
                    threshold = 191504 #172737
                if sum(subSynteny_table['alnlen']) >= threshold:
                    PASS="YES"
                else:
                    PASS="NO"
                #print( "Length homologous :" , sum(subSynteny_table['alnlen']))
                subResult_table=subResult_table.append({"Windows":windows,"Cluster_hmmer":Cluster_hmmer,"Prot_name":Prot_name,"Scaffold_name":EVE_scaffold2, "Species": EVE_species2,"Nb_HSPs": subSynteny_table.shape[0],"Tot_alnlen":sum(subSynteny_table['alnlen']),"PASS":PASS}, ignore_index=True)
                subResult_table
        try:
            Mean_homologous=mean(len_homologous)
        except:
            Mean_homologous=0
        try:
            Median_homologous=median(len_homologous)
        except:
            Median_homologous=0
        #print(len_homologous)
        Result_table=Result_table.append({"Windows":windows,"Cluster_hmmer":Cluster_hmmer,"Prot_name":Prot_name,"Mean_homologous":Mean_homologous,"Median_homologous":Median_homologous,"Tot_alnlen":sum(subSynteny_table['alnlen'])}, ignore_index=True)
        #print(Result_table.loc[[i]])
        i+=1
        count_success=0
        count_success_high=0 
        #print(group_name, " : ", sum(subSynteny_table['alnlen']) , " : ", windows)
        

subResult_table['PASS2'] = (subResult_table['PASS']=="YES").groupby([subResult_table['Cluster_hmmer'], subResult_table['Species']]).transform('any').map({True:"YES", False:"NO"})

#Print number of Cluster,Event without synteny around 

subResult_table=subResult_table.sort_values(by='Tot_alnlen', ascending=False)
subResult_table.drop_duplicates(subset = ['Cluster_hmmer','Species'], keep = 'first',inplace=True)
subResult_table.drop_duplicates(subset = ['Cluster_hmmer'], keep = 'first')['PASS2'].value_counts()

subResult_table.columns=['Synteny_Windows', 'Cluster_hmmer', 'Scaffold_name','Species_name', 'Synteny_Nb_HSPs', 'Synteny_Tot_alnlen', 'Synteny_PASS', 'Prot_name','Synteny_PASS2']
#to get a nice table 
subResult_table.pivot_table(columns='Prot_name',index='Species_name',values="Synteny_PASS2", aggfunc=','.join)


Main_table=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event_dNdS.tab",sep=";")

to_save=pd.merge(Main_table,subResult_table[['Cluster_hmmer','Species_name','Scaffold_name','Synteny_PASS2']], on=['Cluster_hmmer','Species_name','Scaffold_name'],how='left')


# Add start-end of filamentous viruses 

to_save.loc[to_save['Names'].str.contains("Lef5",na=False),"Cluster_hmmer"]='Cluster559-Cluster482-Cluster457-Cluster353-Cluster2746-Cluster2234-Cluster1338-Cluster1713-Cluster1176-Cluster1383-Cluster2222-Cluster2543-Cluster735-Cluster167'

to_save=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event_dNdS_Syntenic.tab",sep=";")

grouped = to_save.groupby(['Names2'])

LbFV_names=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Filamentous_core_genes.tab",sep="\t")
for group_name, group in grouped:
  if group_name in ['LbFV','LhFV','EfFV','CcFV1','CcFV2','PoFV','PcFV']:
    if group_name=="LhFV":
      GFF_tab=pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Genomes/"+group_name+"/Predicted_orfs/Final_ORF_prediction_"+group_name+".gff",sep="\t",comment="#",header=None)
      GFF_tab.columns=['Scaffold_name','soft','type','start','end','point1','strand','point2','Description']
      GFF_tab['Scaffold_name']=GFF_tab['Scaffold_name'].str.replace("%.*","")
      for index, row in group.iterrows():
        if "orf" in row['Names']:
          scaffold_name=row['Names'].split("_",3)[1]+"_"+row['Names'].split("_",3)[2]
          protein_name=row['Names']
          #print(protein_name)
          if group_name =="LbFV":
            if row['Names'] =="Lef5_LbFV":
              scaffold_name="NC_033778.1"
              continue
            else:
              scaffold_name="NC_033778.1"
              protein_name= LbFV_names.loc[LbFV_names['Names'].eq(re.sub("_LbFV","",row['Names']))]['YP_format'].iloc[0]
              New_prot_name=LbFV_names.loc[LbFV_names['Names'].eq(re.sub("_LbFV","",row['Names']))]['Prot_name'].iloc[0]
          subGFF_tab=GFF_tab.loc[GFF_tab['Scaffold_name'].eq(scaffold_name)]
          ORF_number=row['Names'].split("_",3)[3]
          start= GFF_tab.loc[GFF_tab['Description'].str.contains(row['Names'].split("_",3)[1]+"_"+row['Names'].split("_",3)[2]+"_"+row['Names'].split("_",3)[3])]['start'].iloc[0]
          end=GFF_tab.loc[GFF_tab['Description'].str.contains(row['Names'].split("_",3)[1]+"_"+row['Names'].split("_",3)[2]+"_"+row['Names'].split("_",3)[3])]['end'].iloc[0]
          strand=GFF_tab.loc[GFF_tab['Description'].str.contains(row['Names'].split("_",3)[1]+"_"+row['Names'].split("_",3)[2]+"_"+row['Names'].split("_",3)[3])]['strand'].iloc[0]
        else:
          protein_name=row['Names']
          scaffold_name=row['Names'].split("_", 3)[0]+"_"+ row['Names'].split("_", 3)[1]
          ORF_number="NA"
          start=row['Names'].split("_", 3)[2]
          end =row['Names'].split("_", 3)[3]
          strand =row['Names'].split("_", 5)[4]
        record=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Genomes/"+group_name+"/"+group_name+"_free.fna","fasta"))
        Scaff_length=len(record[scaffold_name].seq)
        to_save.loc[to_save['Names'].eq(row['Names']),"Scaffold_length"]=Scaff_length
        #Correct coordinates within filamentous free-living genomes
        to_save.loc[to_save['Names'].eq(row['Names']),"start"]= start
        to_save.loc[to_save['Names'].eq(row['Names']),"end"]= end
        to_save.loc[to_save['Names'].eq(row['Names']),"strand"]=strand
        to_save.loc[to_save['Names'].eq(row['Names']),"Scaffold_name"]=scaffold_name
        to_save.loc[to_save['Names'].eq(row['Names']),"Species_name"]=group_name
    else:
      print(group_name)
      #for index, row in group.iterrows():
      GFF_tab=pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Genomes/"+group_name+"/Predicted_orfs/Final_ORF_prediction_"+group_name+".gff",sep="\t",comment="#",header=None)
      GFF_tab.columns=['Scaffold_name','soft','type','start','end','point1','strand','point2','Description']
      GFF_tab['Scaffold_name']=GFF_tab['Scaffold_name'].str.replace("%.*","")
      #print(GFF_tab)
      for index, row in group.iterrows():
        #print(row['Names'])
        scaffold_name=row['Names']
        scaffold_name=re.sub("_orf.*","",scaffold_name)
        tag= ".*"+group_name+"_"
        scaffold_name=re.sub(tag,"",scaffold_name)
        protein_name=row['Names']
        #print(protein_name)
        if group_name =="LbFV":
          if row['Names'] =="Lef5_LbFV":
            print("ok")
          else:
            scaffold_name="NC_033778.1"
            protein_name= LbFV_names.loc[LbFV_names['Names'].eq(re.sub("_LbFV","",row['Names']))]['YP_format'].iloc[0]
            New_prot_name=LbFV_names.loc[LbFV_names['Names'].eq(re.sub("_LbFV","",row['Names']))]['Prot_name'].iloc[0]
        else:
          subGFF_tab=GFF_tab.loc[GFF_tab['Scaffold_name'].eq(scaffold_name)]
          ORF_number=re.sub(".*_orf","_orf",protein_name)
          ORF_number=re.sub("_putative.*","",ORF_number)
          start=subGFF_tab.loc[subGFF_tab['Description'].str.contains(ORF_number)]['start'].iloc[0]
          end=subGFF_tab.loc[subGFF_tab['Description'].str.contains(ORF_number)]['end'].iloc[0]
          strand=subGFF_tab.loc[subGFF_tab['Description'].str.contains(ORF_number)]['strand'].iloc[0]
          record=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Genomes/"+group_name+"/"+group_name+"_free.fna","fasta"))
          Scaff_length=len(record[scaffold_name].seq)
        if row['Names']=="Lef5_LbFV":
          print(row['Names'])
          start=36651
          end=39301
          strand="-"
          Scaff_length="111453"
        to_save.loc[to_save['Names'].eq(row['Names']),"Scaffold_length"]=Scaff_length
        #Correct coordinates within filamentous free-living genomes
        to_save.loc[to_save['Names'].eq(row['Names']),"start"]= start
        to_save.loc[to_save['Names'].eq(row['Names']),"end"]= end
        to_save.loc[to_save['Names'].eq(row['Names']),"strand"]=strand
        to_save.loc[to_save['Names'].eq(row['Names']),"Scaffold_name"]=scaffold_name
        to_save.loc[to_save['Names'].eq(row['Names']),"Species_name"]=group_name
        if group_name =="LbFV":
          to_save.loc[to_save['Names'].eq(row['Names']),"New_prot_name"]=New_prot_name

f = lambda x: '-'.join(dict.fromkeys(x.dropna()))
to_save['New_prot_name'] = to_save.groupby('Cluster_hmmer')['Prot_name'].transform(f).replace('', np.nan)

to_save.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event_dNdS_Syntenic.tab",sep=";",index=False)


##### Generate dataframe for synteny plot

Output_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event_dNdS_Syntenic.tab"
Main_table=pd.read_csv(Output_file,sep=";")

#Filter
Main_table=Main_table.loc[(Main_table['Cynip_count'].gt(0) & Main_table['Filamentous_count'].gt(0) )| (Main_table['Names2'].isin['LhFV','LbFV']) ]

Main_table=Main_table.loc[Main_table['Names2'].isin(['LbEFV','LhEFV','LcEFV','RhEFV','ThEFV','TrEFV','LbFV','LhFV'])]
Main_table['Type']="EVE"

Main_table.loc[Main_table['Names'].eq("contig_701_7350_6886_-_LhFV"),"start"]="7350"
Main_table.loc[Main_table['Names'].eq("contig_701_7350_6886_-_LhFV"),"end"]="6886"

Main_table['len']=Main_table['end'].astype(float).astype(int)-Main_table['start'].astype(float).astype(int)
Main_table=Main_table[['Names','Names2','start','end','strand','len','Type','Scaffold_name','Cluster_hmmer','Scaffold_length','Prot_name']]
Main_table.columns=['Names','name','start','end','strand','length','code','Scaffold_name','Clusters','Scaffold_length','Protein_name']
Main_table['Species_name'] = Main_table['Names'].str.replace(".*:","")
Main_table['Scaffold_species_name']=Main_table['Species_name']+":"+Main_table['Scaffold_name']

Length_table = pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds_length.tab",sep=";")

Main_table=Main_table.merge(Length_table,on="Scaffold_species_name",how="left")

Main_table['Scaffold_length_y'].fillna(Main_table['Scaffold_length_x'], inplace=True)

Main_table.loc[Main_table['Scaffold_species_name'].eq("Rhoptromeris:scaffold32196|size3216",'Scaffold_length_y')]="size3216"

del Main_table['Scaffold_length_x']

Main_table.rename(columns={'Scaffold_length_y': 'Scaffold_length'}, inplace=True)

Main_table=Main_table.loc[~Main_table['Scaffold_length'].isna()]


#Add all the ORFs for LbFV 
#GFF_tab=pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Genomes/LbFV/Predicted_orfs/Final_ORF_prediction_LbFV.gff",sep="\t",comment="#",header=None)
#GFF_tab.columns=['Scaffold_name','soft','type','start','end','point1','strand','point2','Description']
#GFF_tab=GFF_tab.loc[GFF_tab['soft'].eq("RefSeq") & GFF_tab['type'].eq("CDS")]
#Main_table['start']=Main_table['start'].astype(str)
#GFF_tab['start']=GFF_tab['start'].astype(str)
#GFF_tab=GFF_tab.loc[~GFF_tab['start'].isin(Main_table.loc[Main_table['Scaffold_name'].eq("NC_033778.1")]['start'])]

#Main_table['start']=Main_table['start'].astype(int)
#GFF_tab['start']=GFF_tab['start'].astype(int)
#GFF_tab=GFF_tab[['start','end','strand','Description']]
#GFF_tab['name']="LbFV"
##GFF_tab['length']=GFF_tab['end']-GFF_tab['start']
#GFF_tab['code']="EVE"
#GFF_tab['Scaffold_name']="NC_033778.1"
#GFF_tab['Clusters']="Unknown"
#GFF_tab['Scaffold_length']=111453.0
#GFF_tab['Protein_name']="NA"


#LbFV_names=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Filamentous_core_genes.tab",sep="\t")

#for index, row in LbFV_names.iterrows():
#	if len(GFF_tab.loc[GFF_tab['Description'].str.contains(row['YP_format'])])>0:
#		GFF_tab.loc[GFF_tab['Description'].str.contains(row['YP_format']),"Protein_name"] = row['Prot_name']

#GFF_tab=GFF_tab[['name','start','end','strand','length','code','Scaffold_name','Clusters','Scaffold_length','Protein_name']]
#Main_table=Main_table.append(GFF_tab)

#GFF_tab=pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Genomes/LhFV/Predicted_orfs/Final_ORF_prediction_LhFV.gff",sep="\t",comment="#",header=None)
#GFF_tab.columns=['Scaffold_name','soft','type','start','end','point1','strand','point2','Description']
#GFF_tab=GFF_tab.loc[GFF_tab['soft'].eq("Geneious") & GFF_tab['type'].eq("CDS")]
#Main_table['start']=Main_table['start'].astype(str)
#GFF_tab['start']=GFF_tab['start'].astype(str)
#GFF_tab=GFF_tab.loc[~GFF_tab['start'].isin(Main_table.loc[Main_table['name'].eq("LhFV")]['start'])]

#Main_table['start']=Main_table['start'].astype(int)
#GFF_tab['start']=GFF_tab['start'].astype(int)
#GFF_tab=GFF_tab[['Scaffold_name','start','end','strand','Description']]
#GFF_tab['name']="LhFV"
#GFF_tab['length']=GFF_tab['end']-GFF_tab['start']
#GFF_tab['code']="EVE"
#GFF_tab['Clusters']="Unknown"
#GFF_tab['Scaffold_length']="NA"
#GFF_tab['Protein_name']="NA"

#GFF_tab=GFF_tab[['name','start','end','strand','length','code','Scaffold_name','Clusters','Scaffold_length','Protein_name']]
#Main_table=Main_table.append(GFF_tab)

#Remove transcript scaffolds from LhEFV

Main_table=Main_table.loc[~ (Main_table['Scaffold_name'].str.contains("XM_"))]


# Add genes predictions 

Metaeuk_gff=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/Metaeuk_preds_results.gff",sep="\t",header=None)
Metaeuk_tsv=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/Taxa_Metaeuk_preds_results_tax_per_pred.tsv",sep="\t",header=None)

Metaeuk_gff['Gene_name']= Metaeuk_gff[8].str.replace(";.*","")
Metaeuk_gff['Gene_name']= Metaeuk_gff['Gene_name'].str.replace("Target_ID=","")

Metaeuk_tsv['Gene_name'] = Metaeuk_tsv[0].str.replace("\\|.*","")

Metaeuk_table=Metaeuk_gff.merge(Metaeuk_tsv[[3,4,'Gene_name']],on="Gene_name")
Metaeuk_table.columns=['Species_and_Scaffold','soft','type','start','end','bits','strand','point','description','Gene_name','order','taxlineage']
Metaeuk_table=Metaeuk_table.loc[Metaeuk_table['taxlineage'].str.contains("Euk",na=False) & Metaeuk_table['type'].eq("gene")]
Metaeuk_table=Metaeuk_table.sort_values(by='bits', ascending=False)
Metaeuk_table=Metaeuk_table.loc[Metaeuk_table['bits'].ge(50)]
Metaeuk_table = Metaeuk_table.drop_duplicates(subset=['Species_and_Scaffold', 'start','end','strand'], keep='first')

# Create a fasta file with all the eukaryotic genes in order to run a blastp all vs all to build connexions in the plot 

list_cynipoids=['Leptopilina_boulardi_GCA_019393585.1','Leptopilina_clavipes','Leptopilina_heterotoma_GCA_015476425.1','Rhoptromeris','Trybliographa','Thrichoplasta']

from Bio import SeqIO
import re

with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/MEtaeuk_sequences.dna","w") as output :
	for sp in list_cynipoids:
		record_dict=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+sp+"/assembly/"+sp+"_final_assembly.fna","fasta"))
		for index, row in Metaeuk_table.loc[Metaeuk_table['Species_and_Scaffold'].str.contains(sp)].iterrows():
			if row['strand'] == "+":
				print('>'+row['Species_and_Scaffold'],":",row['start'],"-",row['end'],"("+row['strand'],")",sep="",file=output)
				#print('>',record_dict[re.sub(".*:","",row['Species_and_Scaffold'])].id,sep="")
				print(str(record_dict[re.sub(".*:","",row['Species_and_Scaffold'])].seq[row['start']:row['end']]),file=output)
			if row['strand'] == "-":
				print('>'+row['Species_and_Scaffold'],":",row['start'],"-",row['end'],"("+row['strand'],")",sep="",file=output)
				#print('>',record_dict[re.sub(".*:","",row['Species_and_Scaffold'])].id,sep="")
				print(str(record_dict[re.sub(".*:","",row['Species_and_Scaffold'])].seq[row['start']:row['end']].reverse_complement()),file=output)


# Run blastn all vs all 
import subprocess

subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/MEtaeuk_sequences.dna /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/MEtaeuk_sequences_db", shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  search /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/MEtaeuk_sequences_db /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/MEtaeuk_sequences_db /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/MEtaeuk_sequences_all_vs_all_result /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/MEtaeuk_sequences_all_vs_all_tpm --threads 15 -e 0.000001 --search-type 3", shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs  convertalis  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tlen,qlen,qcov,tcov' --search-type 3 /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/MEtaeuk_sequences_db /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/MEtaeuk_sequences_db /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/MEtaeuk_sequences_all_vs_all_result /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/MEtaeuk_sequences_all_vs_all_result.m8 --threads 15", shell=True)

Metaeuk_blast_tab=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/MEtaeuk_sequences_all_vs_all_result.m8",sep="\t",header=None)
Metaeuk_blast_tab.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']

Metaeuk_blast_tab[['target_species','Scaffold_target','target_coordinates']] = Metaeuk_blast_tab['target'].str.split(":", expand=True)
Metaeuk_blast_tab[['query_species','Scaffold_query','query_coordinates']] = Metaeuk_blast_tab['query'].str.split(":", expand=True)

#remove sefl matchs  and filter 
Metaeuk_blast_tab=Metaeuk_blast_tab.loc[~ (Metaeuk_blast_tab['target_species'] == Metaeuk_blast_tab['query_species'])]
Metaeuk_blast_tab=Metaeuk_blast_tab.loc[Metaeuk_blast_tab['bits'].ge(50)]
Metaeuk_blast_tab=Metaeuk_blast_tab.loc[Metaeuk_blast_tab['alnlen'].ge(100)]
Metaeuk_blast_tab=Metaeuk_blast_tab.loc[Metaeuk_blast_tab['tcov'].ge(0.25)  & Metaeuk_blast_tab['qcov'].ge(0.25) ]

#Only keep results for scaffold containing EVEs 
Metaeuk_blast_tab=Metaeuk_blast_tab.loc[Metaeuk_blast_tab['Scaffold_target'].isin(Main_table['Scaffold_name']) &  Metaeuk_blast_tab['Scaffold_query'].isin(Main_table['Scaffold_name'])]

Metaeuk_blast_tab.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/Metaeuk_blast_table.tab",sep=";",index=False)
# add EUK coordinates to the Main_table 

Metaeuk_table['name']="NA"
Metaeuk_table['name'].mask(Metaeuk_table["Species_and_Scaffold"].str.contains("boulardi"), 'LbEFV', inplace=True)
Metaeuk_table['name'].mask(Metaeuk_table["Species_and_Scaffold"].str.contains("clavipes"), 'LcEFV', inplace=True)
Metaeuk_table['name'].mask(Metaeuk_table["Species_and_Scaffold"].str.contains("heterotoma"), 'LhEFV', inplace=True)
Metaeuk_table['name'].mask(Metaeuk_table["Species_and_Scaffold"].str.contains("Rhoptromeris"), 'RhEFV', inplace=True)
Metaeuk_table['name'].mask(Metaeuk_table["Species_and_Scaffold"].str.contains("Trybliographa"), 'TrEFV', inplace=True)
Metaeuk_table['name'].mask(Metaeuk_table["Species_and_Scaffold"].str.contains("Thrichoplasta"), 'ThEFV', inplace=True)
Metaeuk_table['code']="metaeuk"
Metaeuk_table['Clusters']="no"
Metaeuk_table['Scaffold_length']="NA"
Metaeuk_table['length']=Metaeuk_table['end'] - Metaeuk_table['start']
Metaeuk_table['Scaffold_name']=Metaeuk_table['Species_and_Scaffold'].str.replace(".*:","")

Metaeuk_table=Metaeuk_table.loc[Metaeuk_table['Scaffold_name'].isin(Main_table['Scaffold_name'])]

for index, row in Metaeuk_table.iterrows():
	#print(row['Scaffold_name'])
	Scaffold_length=Main_table.loc[Main_table['Scaffold_name'].eq(row['Scaffold_name'])]['Scaffold_length'].iloc[0]
	Metaeuk_table.loc[Metaeuk_table['Scaffold_name'].eq(row['Scaffold_name']),"Scaffold_length"]=Scaffold_length


Metaeuk_table=Metaeuk_table[['name','start','end','strand','length','code','Scaffold_name','Clusters','Scaffold_length','Gene_name']]
Metaeuk_table.columns=['name', 'start', 'end', 'strand', 'length', 'code', 'Scaffold_name', 'Clusters', 'Scaffold_length', 'Protein_name']


Main_table=Main_table.append(Metaeuk_table)


Main_table.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/Table_for_synteny_plot.txt",sep=";",index=False)


##### Generate the comparison table 
# Open the all vs all huge file 

All_vs_all=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses_result.m8",sep="\t",header=None)
All_vs_all.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']


All_vs_all['query_sp']= All_vs_all['query'].str.replace(".*:","")
All_vs_all['target_sp']= All_vs_all['target'].str.replace(".*:","")


All_vs_all['query_SP']="NA"
All_vs_all['query_SP'].mask(All_vs_all["query"].str.contains("LhFV"), 'LhFV', inplace=True)
All_vs_all['query_SP'].mask(All_vs_all["query"].str.contains("LbFV"), 'LbFV', inplace=True)
All_vs_all['query_SP'].mask(All_vs_all["query"].str.contains("Thrichoplasta"), 'ThEFV', inplace=True)
All_vs_all['query_SP'].mask(All_vs_all["query"].str.contains("Trybliographa"), 'TrEFV', inplace=True)
All_vs_all['query_SP'].mask(All_vs_all["query"].str.contains("Rhoptromeris"), 'RhEFV', inplace=True)
All_vs_all['query_SP'].mask(All_vs_all["query"].str.contains("clavipes"), 'LcEFV', inplace=True)
All_vs_all['query_SP'].mask(All_vs_all["query"].str.contains("boulardi"), 'LbEFV', inplace=True)
All_vs_all['query_SP'].mask(All_vs_all["query"].str.contains("heterotoma"), 'LhEFV', inplace=True)
All_vs_all['target_SP']="NA"
All_vs_all['target_SP'].mask(All_vs_all["target"].str.contains("LhFV"), 'LhFV', inplace=True)
All_vs_all['target_SP'].mask(All_vs_all["target"].str.contains("LbFV"), 'LbFV', inplace=True)
All_vs_all['target_SP'].mask(All_vs_all["target"].str.contains("Thrichoplasta"), 'ThEFV', inplace=True)
All_vs_all['target_SP'].mask(All_vs_all["target"].str.contains("Trybliographa"), 'TrEFV', inplace=True)
All_vs_all['target_SP'].mask(All_vs_all["target"].str.contains("Rhoptromeris"), 'RhEFV', inplace=True)
All_vs_all['target_SP'].mask(All_vs_all["target"].str.contains("clavipes"), 'LcEFV', inplace=True)
All_vs_all['target_SP'].mask(All_vs_all["target"].str.contains("boulardi"), 'LbEFV', inplace=True)
All_vs_all['target_SP'].mask(All_vs_all["target"].str.contains("heterotoma"), 'LhEFV', inplace=True)

All_vs_all=All_vs_all.loc[All_vs_all['target_SP'].isin(['LhFV','LbFV','LbEFV','LcEFV','LhEFV','RhEFV','TrEFV','ThEFV']) & All_vs_all['query_SP'].isin(['LhFV','LbFV','LbEFV','LcEFV','LhEFV','RhEFV','TrEFV','ThEFV'])] 

All_vs_all=All_vs_all.loc[~(All_vs_all['query'].str.contains('PoFV|PcFV|EfFV|CcFV') | All_vs_all['target'].str.contains('PoFV|PcFV|EfFV|CcFV'))]
#Manuall correction stuff 
#All_vs_all['query'].mask(All_vs_all["query"].str.contains("scaffold97118:4995-6326"), 'scaffold97118:4995-7124(+):Trybliographa', inplace=True)
#All_vs_all['target'].mask(All_vs_all["target"].str.contains("scaffold97118:4995-6326"), 'scaffold97118:4995-7124(+):Trybliographa', inplace=True)

#All_vs_all['query'].mask(All_vs_all["query"].str.contains("scaffold7282\\|size12268:1200-1859"), 'sscaffold7282|size12268:1200-2222(-):Rhoptromeris', inplace=True)
#All_vs_all['target'].mask(All_vs_all["target"].str.contains("scaffold7282\\|size12268:1200-1859"), 'scaffold7282|size12268:1200-2222(-):Rhoptromeris', inplace=True)

# Add cluster informations 
Main_table=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event_dNdS_Syntenic.tab",sep=";")
#Main_table=Main_table.loc[Main_table['Cynip_count'].gt(0) & Main_table['Filamentous_count'].gt(0)]

sub1_All_clusters=All_vs_all.merge(Main_table[['Cluster_hmmer','Prot_name','Names']],right_on="Names",left_on="query", how='left')

sub2_All_clusters=sub1_All_clusters.loc[sub1_All_clusters['Cluster_hmmer'].isna()]
del sub2_All_clusters['Cluster_hmmer']
del sub2_All_clusters['Prot_name']
del sub2_All_clusters['Names']
sub2_All_clusters=sub2_All_clusters.merge(Main_table[['Cluster_hmmer','Prot_name','Names']],right_on="Names",left_on="target", how='left')

sub_All_clusters=sub1_All_clusters.append(sub2_All_clusters)
sub_All_clusters = sub_All_clusters.drop_duplicates()

All_clusters=sub_All_clusters.loc[~sub_All_clusters['Cluster_hmmer'].isna()]

#generate start and end bases on the name of the sequences 

import numpy as np 
All_clusters['qstart']=All_clusters['query'].str.replace("\\(.*","")
All_clusters['qstart']=All_clusters['qstart'].str.replace(".*:","")
All_clusters['qend']=All_clusters['qstart'].str.replace(".*-","")
All_clusters['qstart']=All_clusters['qstart'].str.replace("-.*","")
All_clusters['qstrand']=np.where(All_clusters["query"].str.contains("\\+"),'+','-')

All_clusters['tstart']=All_clusters['target'].str.replace("\\(.*","")
All_clusters['tstart']=All_clusters['tstart'].str.replace(".*:","")
All_clusters['tend']=All_clusters['tstart'].str.replace(".*-","")
All_clusters['tstart']=All_clusters['tstart'].str.replace("-.*","")
All_clusters['tstrand']=np.where(All_clusters["target"].str.contains("\\+"),'+','-')

Filamentous_tab=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Filamentous_core_genes.tab",sep="\t")
# Deal with coordinates of LhFV and LbFV 
for sp in ['LhFV','LbFV']:
  print(sp)
  for query in All_clusters.loc[All_clusters['query'].str.contains(sp)]['query'].unique():
    try:
      start=int(Main_table.loc[Main_table['Names'].str.contains(query)]['start'].iloc[0])
      end=int(Main_table.loc[Main_table['Names'].str.contains(query)]['end'].iloc[0])
      strand=Main_table.loc[Main_table['Names'].str.contains(query)]['strand'].iloc[0]
    except:
      GFF_tab=pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Genomes/"+sp+"/Predicted_orfs/Final_ORF_prediction_"+sp+".gff",sep="\t",comment="#",header=None)
      GFF_tab.columns=['Scaffold_name','soft','type','start','end','point1','strand','point2','Description']
      try:
        start=int(GFF_tab.loc[GFF_tab['Description'].str.contains(re.sub(sp+"_","",query))]['start'].iloc[0])
        end=int(GFF_tab.loc[GFF_tab['Description'].str.contains(re.sub(sp+"_","",query))]['end'].iloc[0])
        strand=GFF_tab.loc[GFF_tab['Description'].str.contains(re.sub(sp+"_","",query))]['strand'].iloc[0]
      except:
        try:
          start=query.split("_", 3)[2]
          end =query.split("_", 4)[3]
          strand =query.split("_", 5)[4]
        except:
          if query=="Lef5_LbFV":
            start=36651
            end=39301
            strand="-"
          else:
            query2=Filamentous_tab.loc[Filamentous_tab['Names'].str.contains(re.sub("_.*","",query))]['YP_format'].iloc[0]
            start=int(GFF_tab.loc[GFF_tab['Description'].str.contains(re.sub(sp+"_","",query2))]['start'].iloc[0])
            end=int(GFF_tab.loc[GFF_tab['Description'].str.contains(re.sub(sp+"_","",query2))]['end'].iloc[0])
            strand=GFF_tab.loc[GFF_tab['Description'].str.contains(re.sub(sp+"_","",query2))]['strand'].iloc[0]
    All_clusters.loc[All_clusters['query'].eq(query),"qstart"]=start
    All_clusters.loc[All_clusters['query'].eq(query),"qend"]=end
    All_clusters.loc[All_clusters['query'].eq(query),"qstrand"]=strand
  for target in All_clusters.loc[All_clusters['target'].str.contains(sp)]['target'].unique():
    print(target)
    try:
      start=int(Main_table.loc[Main_table['Names'].str.contains(target)]['start'].iloc[0])
      end=int(Main_table.loc[Main_table['Names'].str.contains(target)]['end'].iloc[0])
      strand=Main_table.loc[Main_table['Names'].str.contains(target)]['strand'].iloc[0]
    except:
      GFF_tab=pd.read_csv("/beegfs/data/bguinet/LbFV_family_project/Genomes/"+sp+"/Predicted_orfs/Final_ORF_prediction_"+sp+".gff",sep="\t",comment="#",header=None)
      GFF_tab.columns=['Scaffold_name','soft','type','start','end','point1','strand','point2','Description']
      try:
        start=int(GFF_tab.loc[GFF_tab['Description'].str.contains(re.sub(sp+"_","",target))]['start'].iloc[0])
        end=int(GFF_tab.loc[GFF_tab['Description'].str.contains(re.sub(sp+"_","",target))]['end'].iloc[0])
        strand=GFF_tab.loc[GFF_tab['Description'].str.contains(re.sub(sp+"_","",target))]['strand'].iloc[0]
      except:
        try:
          start=target.split("_", 3)[2]
          end =target.split("_", 4)[3]
          strand =target.split("_", 5)[4]
        except:
          if target=="Lef5_LbFV":
            start=36651
            end=39301
            strand="-"
          else:
            target2=Filamentous_tab.loc[Filamentous_tab['Names'].str.contains(re.sub("_.*","",target))]['YP_format'].iloc[0]
            start=int(GFF_tab.loc[GFF_tab['Description'].str.contains(re.sub(sp+"_","",target2))]['start'].iloc[0])
            end=int(GFF_tab.loc[GFF_tab['Description'].str.contains(re.sub(sp+"_","",target2))]['end'].iloc[0])
            strand=GFF_tab.loc[GFF_tab['Description'].str.contains(re.sub(sp+"_","",target2))]['strand'].iloc[0]
    #print("start : ", start)
    #print("end : ", end)
    All_clusters.loc[All_clusters['target'].eq(target),"tstart"]=start
    All_clusters.loc[All_clusters['target'].eq(target),"tend"]=end
    All_clusters.loc[All_clusters['target'].eq(target),"tstrand"]=strand

All_clusters.loc[All_clusters['query'].str.contains("AQQ80017")]

All_clusters.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/Synteny_table_for_comparaison_table.txt",sep=";",index=False)
#Generate the comparaison format table from GenoplotR 

All_clusters = pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/Synteny_table_for_comparaison_table.txt",sep=";")
import itertools

x = ["LhFV","LbFV","LbEFV","LcEFV","LhEFV","RhEFV","ThEFV","TrEFV"]

pairwise_list=list(itertools.combinations(x, 2))


Comparaison_table=pd.DataFrame(columns=('species1','start1','end1','species2','start2','end2','name1','name2','scaff1','scaff2','per_id','aln_len','mism','gaps','e_value','bit_score','direction','pair_group','code'))


pairwise_list=[['LbFV','LhFV'],['LbFV','LbEFV'],['LhFV','LbEFV'],['LbEFV','LhEFV'],['LhEFV','LcEFV'],['LcEFV','RhEFV'],['RhEFV','ThEFV'],['ThEFV','TrEFV']]

pairwise_list=[["LbFV","LhFV"],["LhFV","LbEFV"],["LbEFV","LhEFV"],["LhEFV","LcEFV"],["LcEFV","RhEFV"],["RhEFV","ThEFV"],["ThEFV","TrEFV"]]

All_clusters=All_clusters.loc[~All_clusters['Cluster_hmmer'].isin(['Cluster29_redefined'])]

All_clusters['Scaffold_target']=All_clusters['target'].str.replace(":.*","")
All_clusters['Scaffold_target']=All_clusters['Scaffold_target'].str.replace("_orf.*","")
All_clusters['Scaffold_target']=All_clusters['Scaffold_target'].str.replace("LhFV_","")
All_clusters['Scaffold_target'].mask(All_clusters['Scaffold_target'].str.contains("LbFV"),"NC_033778.1",inplace=True)

All_clusters['Scaffold_query']=All_clusters['query'].str.replace(":.*","")
All_clusters['Scaffold_query']=All_clusters['Scaffold_query'].str.replace("_orf.*","")
All_clusters['Scaffold_query']=All_clusters['Scaffold_query'].str.replace("LhFV_","")
All_clusters['Scaffold_query'].mask(All_clusters['Scaffold_query'].str.contains("LbFV"),"NC_033778.1",inplace=True)

Metaeuk_blast_table = pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/Metaeuk_blast_table.tab",sep=";")
Metaeuk_blast_table['query_SP']="NA"
Metaeuk_blast_table['query_SP'].mask(Metaeuk_blast_table["query"].str.contains("Thrichoplasta"), 'ThEFV', inplace=True)
Metaeuk_blast_table['query_SP'].mask(Metaeuk_blast_table["query"].str.contains("Trybliographa"), 'TrEFV', inplace=True)
Metaeuk_blast_table['query_SP'].mask(Metaeuk_blast_table["query"].str.contains("Rhoptromeris"), 'RhEFV', inplace=True)
Metaeuk_blast_table['query_SP'].mask(Metaeuk_blast_table["query"].str.contains("clavipes"), 'LcEFV', inplace=True)
Metaeuk_blast_table['query_SP'].mask(Metaeuk_blast_table["query"].str.contains("boulardi"), 'LbEFV', inplace=True)
Metaeuk_blast_table['query_SP'].mask(Metaeuk_blast_table["query"].str.contains("heterotoma"), 'LhEFV', inplace=True)
Metaeuk_blast_table['target_SP']="NA"
Metaeuk_blast_table['target_SP'].mask(Metaeuk_blast_table["target"].str.contains("Thrichoplasta"), 'ThEFV', inplace=True)
Metaeuk_blast_table['target_SP'].mask(Metaeuk_blast_table["target"].str.contains("Trybliographa"), 'TrEFV', inplace=True)
Metaeuk_blast_table['target_SP'].mask(Metaeuk_blast_table["target"].str.contains("Rhoptromeris"), 'RhEFV', inplace=True)
Metaeuk_blast_table['target_SP'].mask(Metaeuk_blast_table["target"].str.contains("clavipes"), 'LcEFV', inplace=True)
Metaeuk_blast_table['target_SP'].mask(Metaeuk_blast_table["target"].str.contains("boulardi"), 'LbEFV', inplace=True)
Metaeuk_blast_table['target_SP'].mask(Metaeuk_blast_table["target"].str.contains("heterotoma"), 'LhEFV', inplace=True)

Metaeuk_blast_table['qend']= Metaeuk_blast_table['query'].str.replace(".*:","")
Metaeuk_blast_table['qstart']= Metaeuk_blast_table['qend'].str.replace("-.*","")
Metaeuk_blast_table['qend']= Metaeuk_blast_table['qend'].str.replace("\\(.*","")
Metaeuk_blast_table['qend']= Metaeuk_blast_table['qend'].str.replace(".*-","")

Metaeuk_blast_table['tend']= Metaeuk_blast_table['target'].str.replace(".*:","")
Metaeuk_blast_table['tstart']= Metaeuk_blast_table['tend'].str.replace("-.*","")
Metaeuk_blast_table['tend']= Metaeuk_blast_table['tend'].str.replace("\\(.*","")
Metaeuk_blast_table['tend']= Metaeuk_blast_table['tend'].str.replace(".*-","")

for sp in pairwise_list:
  print(sp)
  sub_tab1=All_clusters.loc[All_clusters['query_SP'].str.contains(sp[0]) & All_clusters['target_SP'].str.contains(sp[1])]
  for index, row in sub_tab1.iterrows():
    species1=sp[0]
    start1=row['qstart']
    end1=row['qend']
    name1=row['query']
    scaff1=row['Scaffold_query']
    species2=sp[1]
    start2=row['tstart']
    end2=row['tend']
    name2=row['target']
    scaff2=row['Scaffold_target']
    per_id=row['pident']
    aln_len=row['alnlen']
    mism=row['mismatch']
    gaps=row['gapopen']
    e_value=row['evalue']
    bit_score=row['bits']
    direction=1
    code="EVE"
    Comparaison_table=Comparaison_table.append({"species1":species1,"start1":start1,"end1":end1,'species2':species2,'start2':start2,'end2':end2,'name1':name1,'name2':name2,'scaff1':scaff1,'scaff2':scaff2,'per_id':per_id,'aln_len':aln_len,'mism':mism,'gaps':gaps,'e_value':e_value,'bit_score':bit_score,'direction':direction,'pair_group':sp[0]+"-"+sp[1],"code":code}, ignore_index=True)
  #Comparaison_table = Comparaison_table.reset_index(level=0)
  sub_tab2=All_clusters.loc[All_clusters['query_SP'].str.contains(sp[1]) & All_clusters['target_SP'].str.contains(sp[0])]
  for index2, row2 in sub_tab2.iterrows():
    species1=sp[0]
    start1=row2['tstart']
    end1=row2['tend']
    name1=row2['target']
    scaff1=row2['Scaffold_target']
    start2=row2['qstart']
    species2=sp[1]
    end2=row2['qend']
    name2=row2['query']
    scaff2=row2['Scaffold_query']
    per_id=row2['pident']
    aln_len=row2['alnlen']
    mism=row2['mismatch']
    gaps=row2['gapopen']
    e_value=row2['evalue']
    bit_score=row2['bits']
    direction=1
    code="EVE"
    Comparaison_table=Comparaison_table.append({"species1":species1,"start1":start1,"end1":end1,'species2':species2,'start2':start2,'end2':end2,'name1':name1,'name2':name2,'scaff1':scaff1,'scaff2':scaff2,'per_id':per_id,'aln_len':aln_len,'mism':mism,'gaps':gaps,'e_value':e_value,'bit_score':bit_score,'direction':direction,'pair_group':sp[0]+"-"+sp[1],"code":code}, ignore_index=True)
  sub_tab3=Metaeuk_blast_table.loc[Metaeuk_blast_table['query_SP'].str.contains(sp[0]) & Metaeuk_blast_table['target_SP'].str.contains(sp[1])]
  for index, row in sub_tab3.iterrows():
    species1=sp[0]
    start1=row['qstart']
    end1=row['qend']
    name1=row['query']
    scaff1=row['Scaffold_query']
    species2=sp[1]
    start2=row['tstart']
    end2=row['tend']
    name2=row['target']
    scaff2=row['Scaffold_target']
    per_id=row['pident']
    aln_len=row['alnlen']
    mism=row['mismatch']
    gaps=row['gapopen']
    e_value=row['evalue']
    bit_score=row['bits']
    direction=1
    code="metaeuk"
    Comparaison_table=Comparaison_table.append({"species1":species1,"start1":start1,"end1":end1,'species2':species2,'start2':start2,'end2':end2,'name1':name1,'name2':name2,'scaff1':scaff1,'scaff2':scaff2,'per_id':per_id,'aln_len':aln_len,'mism':mism,'gaps':gaps,'e_value':e_value,'bit_score':bit_score,'direction':direction,'pair_group':sp[0]+"-"+sp[1],"code":code}, ignore_index=True)
    #Comparaison_table = Comparaison_table.reset_index(level=0)
  sub_tab4=Metaeuk_blast_table.loc[Metaeuk_blast_table['query_SP'].str.contains(sp[1]) & Metaeuk_blast_table['target_SP'].str.contains(sp[0])]
  for index2, row2 in sub_tab4.iterrows():
    species1=sp[0]
    start1=row2['tstart']
    end1=row2['tend']
    name1=row2['target']
    scaff1=row2['Scaffold_target']
    start2=row2['qstart']
    species2=sp[1]
    end2=row2['qend']
    name2=row2['query']
    scaff2=row2['Scaffold_query']
    per_id=row2['pident']
    aln_len=row2['alnlen']
    mism=row2['mismatch']
    gaps=row2['gapopen']
    e_value=row2['evalue']
    bit_score=row2['bits']
    direction=1
    code="metaeuk"
    Comparaison_table=Comparaison_table.append({"species1":species1,"start1":start1,"end1":end1,'species2':species2,'start2':start2,'end2':end2,'name1':name1,'name2':name2,'scaff1':scaff1,'scaff2':scaff2,'per_id':per_id,'aln_len':aln_len,'mism':mism,'gaps':gaps,'e_value':e_value,'bit_score':bit_score,'direction':direction,'pair_group':sp[0]+"-"+sp[1],"code":code}, ignore_index=True)

Comparaison_table=Comparaison_table[['start1', 'end1', 'start2', 'end2', 'name1', 'name2','scaff1','scaff2','species1','species2', 'per_id', 'aln_len', 'mism', 'gaps', 'e_value', 'bit_score', 'direction', 'pair_group','code']]
Comparaison_table=Comparaison_table.drop_duplicates()
Comparaison_table.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/Comparison_table.txt",sep="\t",index=False)


#Comparaison_table=Comparaison_table[['start1', 'end1', 'start2', 'end2', 'name1', 'name2','species1','species2', 'per_id', 'aln_len', 'mism', 'gaps', 'e_value', 'bit_score', 'direction', 'pair_group']]
#Comparaison_table.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/Comparison_table.txt",sep="\t",index=False)


