
import os
import argparse
import networkx as nx
import pandas as pd 
import numpy as np

# Print out a message when the program is initiated. traitrise
print('----------------------------------------------------------------------------------------------\n')
print('                        Filter Blast file with Nr                        .\n')
print('----------------------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Create Binary and cluster file ')
parser.add_argument("-Nr", "--Nr_file", help="The Nr file")
parser.add_argument("-o", "--output", help="The updated blast output file")
args = parser.parse_args()


#Usage example python3 /beegfs/data/bguinet/Cynipoidea_project/Scripts/Filter_NR_results.py -Nr /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/NR_analyse/NR_homology_result.m8 -o /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/NR_analyse/Blast_file_HSPs_NR.tab

Nr_file=args.Nr_file
output_cluster_file=args.output


Species_list=['Chelonus_insularis','Trybliographa','Rhoptromeris','Thrichoplasta','Leptopilina_clavipes','Leptopilina_boulardi_GCA_003121605','Leptopilina_heterotoma_GCA_010016045','Leptopilina_heterotoma_GCA_009602685',
'Leptopilina_heterotoma_GCA_009026005','Leptopilina_heterotoma_GCA_009025955','Leptopilina_boulardi_GCA_015476485',
'Leptopilina_boulardi_GCA_011634795','Leptolamina','Phanerotoma','Dolichomitus','Melanaphis_sacchari','Eurytoma_brunniventris','Leptopilina_heterotoma_GCA_015476425.1','Leptopilina_boulardi_GCA_019393585.1','Platygaster_orseoliae','Platygaster_equestris']

output_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/NR_analyse/Blast_file_HSPs_NR.tab"
Nr_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/NR_analyse/NR_homology_result.m8"


#Open and merge all blast files : 

Blast_tab=pd.DataFrame(columns=['query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tlen', 'qlen', 'qcov', 'tcov', 'new_start', 'new_end', 'strand', 'count_duplicate', 'diff_length_target', 'diff_length_query', 'HSP_group', 'HSP_min_target', 'HSP_max_target', 'HSP_min_query', 'HSP_max_query', 'Overlapp_HSP_query_group', 'HSP_close_query_groups', 'Overlapp_HSP_target_group', 'count_Overlapp_HSP_target_group', 'Overlapp_group', 'Tot_length', 'full_name']) 
for species in Species_list:
  subtab=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/run_mmseqs2/Mmseqs_dsDNA_EVEs_results_filtred.m8",sep=";")
  Blast_tab=pd.concat([Blast_tab,subtab])

## NR blast file 
Nr_blast="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/NR_analyse/NR_homology_result.m8"
NR_blast=pd.read_csv(Nr_blast,sep="\t",header=None)

NR_blast.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qcov','tcov','taxid','taxname','taxlineage']
Number_loci=len(NR_blast['query'].unique())
#Remove ichnovirus and wrongly annotated nudiviruses 
NR_blast=NR_blast.loc[~NR_blast['taxname'].str.contains("ichno")]

NR_blast=NR_blast.loc[~NR_blast['taxname'].str.contains("Cotesia|Microplitis|Chelonus")]


#Sort evalues 
NR_blast=NR_blast.sort_values(by='evalue', ascending=True)

# We keep only one species for each query 
NR_blast = NR_blast.drop_duplicates(subset=['query', 'taxname'], keep='first')

NR_blast=NR_blast.loc[NR_blast['bits'].ge(50)]

#Take into account only non-hymenoptera matches 
NR_blast['count_Bacteria'] = NR_blast['taxlineage'].str.contains('Bacter').groupby(NR_blast['query']).transform('sum')
NR_blast['count_Archaea'] = NR_blast['taxlineage'].str.contains('Archaea').groupby(NR_blast['query']).transform('sum')
NR_blast['count_Eukaryota'] = NR_blast['taxlineage'].str.contains('Eukaryota').groupby(NR_blast['query']).transform('sum')
NR_blast['count_Virus'] = NR_blast['taxlineage'].str.contains('Virus').groupby(NR_blast['query']).transform('sum')
NR_blast['Count_Bacteria_Eukaryota_Archeae'] = NR_blast['count_Eukaryota'] + NR_blast['count_Bacteria']+ NR_blast['count_Archaea']

NR_blast['Perc_viral_vs_other'] = NR_blast['count_Virus'] / NR_blast['Count_Bacteria_Eukaryota_Archeae']

NR_blast['Viral_confidence'] = "NA"

NR_blast.loc[NR_blast['Count_Bacteria_Eukaryota_Archeae'].eq(0),"Viral_confidence"] = "Viral"
NR_blast.loc[NR_blast['Count_Bacteria_Eukaryota_Archeae'].gt(10) & NR_blast['Perc_viral_vs_other'].gt(1),"Viral_confidence"] = "Uncertain"
NR_blast.loc[NR_blast['Perc_viral_vs_other'].gt(1) & NR_blast['Viral_confidence'].str.contains("NA") ,"Viral_confidence"] = "Viral"
NR_blast.loc[NR_blast['Perc_viral_vs_other'].lt(1) & NR_blast['count_Virus'].eq(0) ,"Viral_confidence"] = "Not_viral"
NR_blast.loc[NR_blast['count_Virus'].eq(0)  & NR_blast['Count_Bacteria_Eukaryota_Archeae'].eq(0) ,"Viral_confidence"] = "Viral"
NR_blast.loc[NR_blast['Viral_confidence'].str.contains("NA") ,"Viral_confidence"] = "Uncertain"
NR_blast.loc[NR_blast['count_Virus'].lt(10)  & NR_blast['Count_Bacteria_Eukaryota_Archeae'].gt(50) ,"Viral_confidence"] = "Not_viral"

NR_blast=NR_blast.loc[~NR_blast['Viral_confidence'].eq("Viral")]

m = NR_blast['taxlineage'].str.contains('Virus')

NR_blast = NR_blast.reset_index()

NR_blast['evalue']=NR_blast['evalue'].replace([0.0,], 2.747000e-310)

NR_blast['ln_evalue'] = NR_blast['evalue'].apply(np.log10)

# subm=subNR_blast['taxlineage'].str.contains('Virus')
#m = NR_blast['taxlineage'].str.contains('Virus')

subNR_blast1=NR_blast.loc[NR_blast['count_Virus'].gt(0) & NR_blast['Count_Bacteria_Eukaryota_Archeae'].gt(0)]
m = subNR_blast1['taxlineage'].str.contains('Virus')
subNR_blast1['diff_evalue'] = (subNR_blast1.groupby(['query'])['ln_evalue'].transform(lambda d: d.mask(m).min() -  d.where(m).min() ))

subNR_blast2=NR_blast.loc[~(NR_blast['count_Virus'].gt(0) & NR_blast['Count_Bacteria_Eukaryota_Archeae'].gt(0))]
#NR_blast['diff_evalue'] = (NR_blast.groupby(['query'])['ln_evalue'].transform(lambda d: d.mask(m).min() -  d.where(m).min() ))

subNR_blast1=subNR_blast1.loc[~subNR_blast1['diff_evalue'].gt(0)]

NR_blast=pd.concat([subNR_blast1,subNR_blast2])

# Save NR filtred results :
NR_blast.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/NR_analyse/NR_homology_result_filtred.m8",sep=";",index=False)

print("Number loci removed : ", len(NR_blast['query'].unique()) , " / ",  Number_loci, " (",len(NR_blast['query'].unique())/Number_loci*100,"%)",sep="")

# Remove all remaning non_viral loci 
Blast_tab2=Blast_tab.loc[~Blast_tab['full_name'].isin(NR_blast['query'].unique())]
print(len(Blast_tab2))
Blast_tab2=Blast_tab2.merge(NR_blast[['query','diff_evalue','Viral_confidence']], right_on="query",left_on="full_name",how="left")
del Blast_tab2['query_y']
Blast_tab2.rename(columns={'query_x': 'query'}, inplace=True)
print(len(Blast_tab2))
Blast_tab2.to_csv(output_file,sep=";",index=False)

print("Output Nr filtred file written to ", output_file)


