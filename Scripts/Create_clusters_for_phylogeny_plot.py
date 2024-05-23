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

# example usage
"""
python3 Create_clusters_for_phylogeny_plot.py -b /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event_dNdS.tab
"""

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

Blast_table.loc[Blast_table['Names'].str.contains("MdBV",na=False),"Names2"]="MdENV"
Blast_table.loc[Blast_table['Names'].str.contains("CcBV",na=False),"Names2"]="CcENV"

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
Blast_table=Blast_table.loc[~Blast_table['Names'].str.contains("orseo")]

Blast_table['strand'] = np.where(Blast_table['Names'].str.contains("\\+"), "+", "-")
Blast_table['start-end']=Blast_table['Names'].str.replace("\\(.*","")
Blast_table['start-end']=Blast_table['start-end'].str.replace(".*:","")
Blast_table['start']=Blast_table['start-end'].str.replace("-.*","")                  
Blast_table['end']=Blast_table['start-end'].str.replace(".*-","")

#Blast_table=Blast_table.loc[Blast_table['Cluster_hmmer'].str.contains("279")]


#Load NR results 
Nr_tab=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/NR_analyse/NR_homology_result.m8",sep="\t",header=None)

Nr_tab.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qcov','tcov','taxid','taxname','taxlineage']

#Sort evalues 
Nr_tab=Nr_tab.sort_values(by='evalue', ascending=True)

# We keep only one species for each query 
Nr_tab = Nr_tab.drop_duplicates(subset=['query', 'taxname'], keep='first')

Nr_tab=Nr_tab.loc[Nr_tab['bits'].ge(50)]

sub_Nr_tab=Nr_tab.loc[Nr_tab['query'].isin(Blast_table['Names'].unique())]
#Remove seq already present in analysis 

toremove=['Autographa californica nucleopolyhedrovirus','Autographa californica multiple nucleopolyhedrovirus','Lymantria xylina nucleopolyhedrovirus','Cydia pomonella granulovirus',
'Culex nigripalpus nucleopolyhedrovirus','Heliothis zea nudivirus','Helicoverpa zea nudivirus 2','Gryllus bimaculatus nudivirus','Oryctes rhinoceros nudivirus',
'Tipula oleracea nudivirus','Drosophila innubila nudivirus','Drosophila-associated filamentous virus','Kallithea virus','Esparto virus','Tomelloso virus','Mauternbach virus',
'Penaeus monodon nudivirus','Homarus gammarus nudivirus','Dikerogammarus haemobaphes nudivirus','Glossina pallidipes salivary gland hypertrophy virus','Musca domestica salivary gland hypertrophy virus',
'Leptopilina boulardi filamentous virus','Rachiplusia ou MNPV']

sub_Nr_tab=sub_Nr_tab.loc[~sub_Nr_tab['taxname'].isin(toremove)]
sub_Nr_tab['taxname2']=sub_Nr_tab['taxname'].str.replace(" ","_")
sub_Nr_tab['taxname2']=sub_Nr_tab['taxname2'].str.split("_", expand=True)[0]+ '_'  + sub_Nr_tab['taxname2'].str.split("_", expand=True)[1]
sub_Nr_tab['genus']=sub_Nr_tab['taxname2'].str.split("_", expand=True)[0]
sub_Nr_tab['family'] = sub_Nr_tab['taxlineage'].str.replace(".*;f_","") 
sub_Nr_tab['family'] = sub_Nr_tab['family'].str.replace(";.*","") 
sub_Nr_tab['order'] = sub_Nr_tab['taxlineage'].str.replace(".*;o_","") 
sub_Nr_tab['order'] = sub_Nr_tab['order'].str.replace(";.*","") 

sub_Nr_tab.loc[~sub_Nr_tab['order'].str.contains("-_cellular organisms"),"order")="NA"
sub_Nr_tab.loc[~sub_Nr_tab['family'].str.contains("ae"),"family"]="NA"

sub_Nr_tab = sub_Nr_tab.drop_duplicates(subset=['query', 'taxname2'], keep='first')

from Bio import Entrez
Entrez.email = "benjamin.guinet95@gmail.com"
Entrez.api_key = "30bf99cff0e43d6827934fa6ab127f3b5f09"


Added_sequences1=[]
Added_sequences2=[]
Blast_table_grouped = Blast_table.groupby('Cluster_hmmer')
for group_name, df_group in Blast_table_grouped:
  # if there are viruses and more than two sequence within a cluster, create it !
  genus=[]
  family=[]
  if df_group['Virus_count'].iloc[0]>0:
          if len(df_group['Names'].unique())>=2:
           with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_plot.aa","w") as output_aa:
            print('Processing ... ',group_name)
            print(family)
            for row_index, row in sub_Nr_tab.loc[sub_Nr_tab['query'].isin(df_group['Names'])].iterrows():
             if row['genus'] in genus:
              continue
             else:
              if row['family'] in family:
               continue
              else:
               #print(row)
               #Add NR results 
               try:
                handle=Entrez.efetch(db="protein", id=row['target'], rettype="fasta_cds_na", retmode="text")
                record = SeqIO.read(handle, "fasta")
                if "Virus" in row['taxlineage']:
                 add=" [VIRUS]"
                elif "Euk" in row['taxlineage']:
                 add=" [EUKARYOTA]"
                elif "Bacteria" in row['taxlineage']:
                 add=" [BACTERIA]"
                else:
                 add=" [Unknown]"
                print('>',re.sub(" ","_",row['taxname'])+" [" +row['order']+']' +add,sep="",file=output_aa)
                print(str(record.seq.translate()),file=output_aa)
                genus.append(row['genus'])
                if "Hymeno" in row['order']:
                 continue
                else:
                 family.append(row['family'])
               except:
                print(row['target'], " did not work...")
            print("\n")
            #print("writting cluster : ", group_name)
            # EVE protein loci 
            # ALL sequences
            #with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_plot.aa","w") as output_aa:
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
                        print('>',row['Names']+ " [VIRUS]",sep="",file=output_aa)
                        print(AA_records_origin[row['Names']].seq,file=output_aa)
                        Added_sequences1.append(row['Names'])

print ("ALL AA files with all sequences written to : /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/ with extension plot.aa")
