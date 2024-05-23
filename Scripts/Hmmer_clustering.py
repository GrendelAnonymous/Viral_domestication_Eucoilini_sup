
import pandas as pd 
import sys 
import argparse
import os
from collections import defaultdict
import pandas as pd
from Bio import SearchIO
import re 

# Print out a message when the program is initiated.
print('------------------------------------------------------------------------------\n')
print('                        Gather HMMEr files                                   .\n')
print('------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Gather HMMEr files')
parser.add_argument("-p", "--hmmer_path", help="The path with all hmmer files")
parser.add_argument("-b", "--blast_cluster_file", help="The blast cluster file")
parser.add_argument("-o", "--output_file", help="The hmmer cluster output  file")

#Example usage : 

"""
#
python3 Hmmer_clustering.py -p /beegfs/data/bguinet/Cynipoidea_project/Clustering/Cluster_hmmer -b /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-cluster.tab \
    -o /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster.tab
"""

# Variable that stores fasta sequences
args = parser.parse_args()
path=args.hmmer_path
blast=args.blast_cluster_file
output=args.output_file

path="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_hmmer"
blast="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-cluster.tab"
output="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster.tab"


attribs = ['accession', 'bias', 'bitscore', 'description', 'cluster_num', 'domain_exp_num',  'domain_included_num', 'domain_obs_num', 'domain_reported_num', 'env_num', 'evalue', 'id', 'overlap_num', 'region_num']


Main_Hmmer_tab = pd.DataFrame(columns=attribs)
#open the hmmer result tables
for filename in os.listdir(path):
    if filename.endswith(".tab"):
      change="no"
      hits = defaultdict(list)
      with open(path+"/"+filename) as handle:
          for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
            #print(queryresult.id)
            #print(queryresult.accession)
            #print(queryresult.description)
            for hit in queryresult.hits:
              for attrib in attribs:
                hits[attrib].append(getattr(hit, attrib))
      hmmer_tab=pd.DataFrame.from_dict(hits)
      if "redefined" in filename:
        change="yes"
      else:
        change="no"
      Cluster_name=re.sub("Hmmer_","",filename)
      Cluster_name=re.sub("_.*","",Cluster_name)
      if change =="yes":
        Cluster_name=Cluster_name+"_redefined"
      hmmer_tab['Cluster']=Cluster_name
      Main_Hmmer_tab=pd.concat([Main_Hmmer_tab,hmmer_tab])
      #print(hmmer_tab)

#


#Re-assign the clusters depending on hmmer results 
Main_Hmmer_tab_save=Main_Hmmer_tab.copy()
Main_Hmmer_tab_save2=Main_Hmmer_tab.copy()
Blast_table=pd.read_csv(blast,sep=";")
# Remove self cluster matchs
Blast_table=Blast_table[['Cluster', 'Names', 'Length', 'Names2', 'Number_unique_Names2', 'Nb_duplicated', 'Ratio_duplicated_unique']]

Blast_table.columns=['Cluster_blast', 'Names', 'Length', 'Names2', 'Number_unique_Names2', 'Nb_duplicated', 'Ratio_duplicated_unique']



Main_Hmmer_tab=Main_Hmmer_tab.merge(Blast_table[['Cluster_blast', 'Names']],right_on="Names",left_on="id",how="left")
Main_Hmmer_tab['Cluster']=Main_Hmmer_tab['Cluster'].str.replace(".tab","")
Main_Hmmer_tab=Main_Hmmer_tab.loc[ ~ (Main_Hmmer_tab['Cluster'] == Main_Hmmer_tab['Cluster_blast'])]


Main_Hmmer_tab=Main_Hmmer_tab.loc[Main_Hmmer_tab['evalue'].lt(8.0300000e-05)] #since in control the lowest = 7.300000e-05

# If to much duplicated free-living viruses, then do not merge using Hmmer ! 


# Manually add the ORF72 LbFV EVE within the long read Lboulardi assembly : JADEYJ010000214.1:5429830-5430153(-):Leptopilina_boulardi_GCA_019393585.1 and put the same information as in PHTE01000008.1:237233-237979(+):Leptopilina_boulardi_GCA_011634795 


# Look for non-assigned sequence within cluster and add them within cluster 
New_hmmer_Blast_table=Main_Hmmer_tab.loc[Main_Hmmer_tab['Cluster'].str.contains("Cluster") & Main_Hmmer_tab['Cluster_blast'].isna()]
New_hmmer_Blast_table['Length']="NA"
New_hmmer_Blast_table['Number_unique_Names2']="NA"
New_hmmer_Blast_table['Nb_duplicated']="NA"
New_hmmer_Blast_table['Ratio_duplicated_unique']="NA"

import numpy as np 
New_hmmer_Blast_table['Names2'] = np.where(New_hmmer_Blast_table.id.str.contains(":") , New_hmmer_Blast_table['id'].str.replace(".*:",""), New_hmmer_Blast_table['id'].str.replace(".*_",""))

New_hmmer_Blast_table=New_hmmer_Blast_table[['Cluster','id','Length','Names2','Number_unique_Names2', 'Nb_duplicated','Ratio_duplicated_unique']]
New_hmmer_Blast_table.columns=['Cluster_blast', 'Names', 'Length', 'Names2', 'Number_unique_Names2', 'Nb_duplicated', 'Ratio_duplicated_unique']

Blast_table=Blast_table.append(New_hmmer_Blast_table)

Main_Hmmer_tab=Main_Hmmer_tab.loc[~(Main_Hmmer_tab['Cluster'].str.contains("Cluster") & Main_Hmmer_tab['Cluster_blast'].isna())]
import networkx as nx

G = nx.from_pandas_edgelist(Main_Hmmer_tab, 'Cluster', 'Cluster_blast')
d = {k: '-'.join(c) for c in nx.connected_components(G) for k in c}

# Remove self cluster matchs
Blast_table['Cluster_hmmer'] = Blast_table['Cluster_blast'].replace(d)


# Find wrongly assembled clusters 
Blast_table2=Blast_table.copy()
g = Blast_table2.groupby('Cluster_hmmer')
# get unique values
Blast_table2['Number_unique_Names2'] = g['Names2'].transform('nunique')
# get non-duplicates
non_dup = g['Names2'].transform(lambda x: (~x.duplicated(False)).sum())
# duplicates = unique - non-duplicates
Blast_table2['Nb_duplicated'] = Blast_table2['Number_unique_Names2'] - non_dup


Blast_table2['Ratio_duplicated_unique'] = Blast_table2['Nb_duplicated']/Blast_table2['Number_unique_Names2']
Blast_table_duplicated_clusters=Blast_table2.loc[Blast_table2['Ratio_duplicated_unique'].gt(0.40) & Blast_table2['Nb_duplicated'].ge(4)]

# Re-define more stringent clusters for those sequences 

Main_Hmmer_tab_save2=Main_Hmmer_tab_save2.merge(Blast_table_duplicated_clusters[['Cluster_blast', 'Names']],right_on="Names",left_on="id",how="right")
Main_Hmmer_tab_save2['Cluster']=Main_Hmmer_tab_save2['Cluster'].str.replace(".tab","")
Main_Hmmer_tab_save2=Main_Hmmer_tab_save2.loc[ ~ (Main_Hmmer_tab_save2['Cluster'] == Main_Hmmer_tab_save2['Cluster_blast'])]

Main_Hmmer_tab_save2=Main_Hmmer_tab_save2.loc[Main_Hmmer_tab_save2['evalue'].lt(1e-40)] 
Main_Hmmer_tab_save2=Main_Hmmer_tab_save2.loc[~(Main_Hmmer_tab_save2['Cluster_blast'].str.contains("Cluster") & Main_Hmmer_tab_save2['Cluster'].isna())]

if len(Main_Hmmer_tab_save2)==0:
  Blast_table_duplicated_clusters['Cluster_hmmer']  = Blast_table_duplicated_clusters['Cluster_blast']
else:
  G2 = nx.from_pandas_edgelist(Main_Hmmer_tab_save2, 'Cluster', 'Cluster_blast')
  d2 = {k: '-'.join(c) for c in nx.connected_components(G2) for k in c}
  Blast_table_duplicated_clusters['Cluster_hmmer'] = Blast_table_duplicated_clusters['Cluster_blast'].replace(d2)


#
# Add new result to previous Blast file 

Blast_table=Blast_table[~Blast_table['Cluster_blast'].isin(Blast_table_duplicated_clusters['Cluster_blast'])]
Blast_table=Blast_table.append(Blast_table_duplicated_clusters)



blast_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses_result.m8"


df=pd.read_csv(blast_file,sep="\t",header=None)
df.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qcov','tcov','qlen']
df2=df[['target','query','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','qlen','qcov','tcov','tlen']]
df2.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qcov','tcov','qlen']
df3=df.append(df2)
df3=df3.sort_values(by=['evalue','bits'], ascending=[True,False])
df3=df3.loc[~(df3['query'] == df3['target'])]
df3 = df3.drop_duplicates(subset=['query'], keep='first')
df3=df3[['query','qlen']]
df3.columns=['Names','Length']
del Blast_table['Length']
Blast_table=Blast_table.merge(df3,on="Names",how="left")



# Check duplicates 
Blast_table2=Blast_table.copy()
g = Blast_table2.groupby('Cluster_hmmer')
# get unique values
Blast_table2['Number_unique_Names2'] = g['Names2'].transform('nunique')
# get non-duplicates
non_dup = g['Names2'].transform(lambda x: (~x.duplicated(False)).sum())
# duplicates = unique - non-duplicates
Blast_table2['Nb_duplicated'] = Blast_table2['Number_unique_Names2'] - non_dup

Blast_table2['Ratio_duplicated_unique'] = Blast_table2['Nb_duplicated']/Blast_table2['Number_unique_Names2']

Blast_table_duplicated_clusters=Blast_table2.loc[Blast_table2['Ratio_duplicated_unique'].gt(0.40) & Blast_table2['Nb_duplicated'].ge(4)]


# Re-define more stringent clusters for those sequences 

Main_Hmmer_tab_save=Main_Hmmer_tab_save.merge(Blast_table_duplicated_clusters[['Cluster_blast', 'Names']],right_on="Names",left_on="id",how="right")
Main_Hmmer_tab_save['Cluster']=Main_Hmmer_tab_save['Cluster'].str.replace(".tab","")
Main_Hmmer_tab_save=Main_Hmmer_tab_save.loc[ ~ (Main_Hmmer_tab_save['Cluster'] == Main_Hmmer_tab_save['Cluster_blast'])]

Main_Hmmer_tab_save=Main_Hmmer_tab_save.loc[Main_Hmmer_tab_save['evalue'].lt(1e-50)] 
Main_Hmmer_tab_save=Main_Hmmer_tab_save.loc[~(Main_Hmmer_tab_save['Cluster_blast'].str.contains("Cluster") & Main_Hmmer_tab_save['Cluster'].isna())]

if len(Main_Hmmer_tab_save)==0:
  Blast_table_duplicated_clusters['Cluster_hmmer']  = Blast_table_duplicated_clusters['Cluster_blast']
else:
  G2 = nx.from_pandas_edgelist(Main_Hmmer_tab_save, 'Cluster', 'Cluster_blast')
  d2 = {k: '-'.join(c) for c in nx.connected_components(G2) for k in c}
  Blast_table_duplicated_clusters['Cluster_hmmer'] = Blast_table_duplicated_clusters['Cluster_blast'].replace(d2)

# Add new result to previous Blast file 

Blast_table=Blast_table[~Blast_table['Cluster_blast'].isin(Blast_table_duplicated_clusters['Cluster_blast'])]
Blast_table=Blast_table.append(Blast_table_duplicated_clusters)

g = Blast_table.groupby('Cluster_hmmer')
# get unique values
Blast_table['Number_unique_Names2'] = g['Names2'].transform('nunique')
# get non-duplicates
non_dup = g['Names2'].transform(lambda x: (~x.duplicated(False)).sum())
# duplicates = unique - non-duplicates
Blast_table['Nb_duplicated'] = Blast_table['Number_unique_Names2'] - non_dup

Blast_table['Ratio_duplicated_unique'] = Blast_table['Nb_duplicated']/Blast_table['Number_unique_Names2']

# Manually add the ORF72 LbFV EVE within the long read Lboulardi assembly : >JADEYJ010000123.1:877336-878082(-):Leptopilina_boulardi_GCA_019393585.1 and put the same information as in JABAIF010000692.1:288490-289236(+):Leptopilina_boulardi_GCA_015476485
#Cluster_Hmmer_to_asign=Blast_table.loc[Blast_table['Names'].str.contains("JABAIF010000692.1:288490-289236")]['Cluster_hmmer'].iloc[0]
#Cluster_Blast_to_asign=Blast_table.loc[Blast_table['Names'].str.contains("JABAIF010000692.1:288490-289236")]['Cluster_blast'].iloc[0]
#Length_to_assign=249

#Blast_table=Blast_table.append({'Cluster_blast':Cluster_Blast_to_asign, 'Names':"JADEYJ010000123.1:877336-878082(-):Leptopilina_boulardi_GCA_019393585.1", 'Names2':"LbEFV", 'Number_unique_Names2':12, 'Nb_duplicated':0, 'Ratio_duplicated_unique':0.0, 'Cluster_hmmer':Cluster_Hmmer_to_asign, 'Length':Length_to_assign}, ignore_index=True)

# Manually add the ORF72 LbFV EVE within the long read  Lheterotoma assembly : NW_025111203.1:2150605-2150865(+):Leptopilina_heterotoma_GCA_015476425.1 and put the same information as in   QYUB01233969.1:605-1352(-):Leptopilina_heterotoma_GCA_009025955
#Cluster_Hmmer_to_asign=Blast_table.loc[Blast_table['Names'].str.contains("QYUB01233969.1:605-1352")]['Cluster_hmmer'].iloc[0]
#Cluster_Blast_to_asign=Blast_table.loc[Blast_table['Names'].str.contains("QYUB01233969.1:605-1352")]['Cluster_blast'].iloc[0]
#Length_to_assign=87


#Bast_table=Blast_table.append({'Cluster_blast':Cluster_Blast_to_asign, 'Names':"NW_025111203.1:2150605-2150865(+):Leptopilina_heterotoma_GCA_015476425.1", 'Names2':"LhEFV", 'Number_unique_Names2':12, 'Nb_duplicated':0, 'Ratio_duplicated_unique':0.0, 'Cluster_hmmer':Cluster_Hmmer_to_asign, 'Length':Length_to_assign}, ignore_index=True)

# Manually correct two ORFs caused by  sequencing mistakes 
#Trybliographa

#Blast_table.loc[Blast_table['Names'].str.contains("scaffold97118:4995-6326"),"Length"]=710
#Blast_table.loc[Blast_table['Names'].str.contains("scaffold97118:4995-6326"),"Names"]="scaffold97118:4995-7124(+):Trybliographa"
#Blast_table=Blast_table.loc[~Blast_table['Names'].str.contains("scaffold97118:6330-7124")]

#Rhoptromeris 
#Blast_table.loc[Blast_table['Names'].str.contains("scaffold7282\\|size12268:1200-1859"),"Length"]=341
#Blast_table.loc[Blast_table['Names'].str.contains("scaffold7282\\|size12268:1200-1859"),"Names"]="scaffold7282|size12268:1200-2222(-):Rhoptromeris"
#Blast_table=Blast_table.loc[~Blast_table['Names'].str.contains("scaffold7282\\|size12268:1878-2222")]



#####




# Add best hits informations for each queries 
#Revious blast all vs all file 

All_vs_All=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses_result.m8",sep="\t",header=None)
All_vs_All.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qcov','tcov','qlen']


All_vs_All['query2'] = np.where(All_vs_All['query'].str.contains("_orf") , All_vs_All['query'].str.replace("_.*",""),All_vs_All['query'])
All_vs_All['query2'] = np.where(All_vs_All['query'].str.contains(":") , All_vs_All['query2'],All_vs_All['query2'].str.replace(".*_",""))
All_vs_All['query2'] = np.where(All_vs_All['query'].str.contains(":") ,All_vs_All['query2'].str.replace(".*:",""), All_vs_All['query2'])

All_vs_All['target2'] = np.where(All_vs_All['target'].str.contains("_orf") , All_vs_All['target'].str.replace("_.*",""),All_vs_All['target'])
All_vs_All['target2'] = np.where(All_vs_All['target'].str.contains(":") , All_vs_All['target2'],All_vs_All['target2'].str.replace(".*_",""))
All_vs_All['target2'] = np.where(All_vs_All['target'].str.contains(":") ,All_vs_All['target2'].str.replace(".*:",""), All_vs_All['target2'])



#Remove self matchs  and keep in priority  match with free-living viruses 
All_vs_All=All_vs_All.loc[~(All_vs_All['query'] == All_vs_All['target'])]

All_EVEs=['MdBV', 'Chelonus_insularis',
       'Leptopilina_heterotoma_GCA_009602685',
       'Leptopilina_heterotoma_GCA_010016045', 'Melanaphis_sacchari',
       'Leptopilina_boulardi_GCA_019393585.1',
       'Leptopilina_heterotoma_GCA_015476425.1', 'Dolichomitus',
       'Leptopilina_boulardi_GCA_003121605',
       'Leptopilina_boulardi_GCA_011634795', 'Leptopilina_clavipes',
       'Leptopilina_heterotoma_GCA_009025955',
       'Leptopilina_heterotoma_GCA_009026005', 'Thrichoplasta',
       'Trybliographa', 'Eurytoma_brunniventris',
       'Leptopilina_boulardi_GCA_015476485',
       'VcBV', 'Phanerotoma',
       'Rhoptromeris', 'Leptolamina','FaBV', 
       'BtBV', 'CcBV',
       'PoEFV', 'NlBV', 'CiBV',
       'ApBV']

#

from heapq import merge 
list_loci1=All_vs_All['query'].unique()
list_loci=list(merge(list_loci1,All_vs_All['target'].unique()))

All_vs_All_bis=All_vs_All.loc[~(All_vs_All['query2'].isin(All_EVEs) & All_vs_All['target2'].isin(All_EVEs))]
All_vs_All_bis1=All_vs_All_bis.copy()

All_vs_All_bis1.columns=['target', 'query', 'pident', 'alnlen', 'mismatch', 'gapopen', 'tstart', 'tend', 'qstart', 'qend', 'evalue', 'bits', 'qlen', 'tlen', 'tcov', 'qcov', 'target2', 'query2']
All_vs_All_bis1=All_vs_All_bis1[['query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tlen', 'qlen', 'qcov', 'tcov', 'query2', 'target2']]
All_vs_All_bis=All_vs_All_bis.append(All_vs_All_bis1)
All_vs_All_bis=All_vs_All_bis.sort_values(by='evalue', ascending=True)
All_vs_All_bis = All_vs_All_bis.drop_duplicates(subset = "query",keep="first")


#Deal for loci that have only match with EVEs
filtred_loci1=All_vs_All_bis['query'].unique()
filtred_loci=list(merge(filtred_loci1,All_vs_All_bis['target'].unique()))

List_loci_without_FL_matches=set(list_loci) - set(filtred_loci)
All_vs_All_without_FL_matches=All_vs_All.loc[All_vs_All['query'].isin(List_loci_without_FL_matches) | All_vs_All['target'].isin(List_loci_without_FL_matches)]
All_vs_All_without_FL_matches=All_vs_All_without_FL_matches.sort_values(by='evalue', ascending=True)
All_vs_All_without_FL_matches1=All_vs_All_without_FL_matches.copy()

All_vs_All_without_FL_matches1.columns=['target', 'query', 'pident', 'alnlen', 'mismatch', 'gapopen', 'tstart', 'tend', 'qstart', 'qend', 'evalue', 'bits', 'qlen', 'tlen', 'tcov', 'qcov', 'target2', 'query2']
All_vs_All_without_FL_matches1=All_vs_All_without_FL_matches1[['query', 'target', 'pident', 'alnlen', 'mismatch', 'gapopen', 'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'tlen', 'qlen', 'qcov', 'tcov', 'query2', 'target2']]
All_vs_All_without_FL_matches=All_vs_All_without_FL_matches.append(All_vs_All_without_FL_matches1)
All_vs_All_without_FL_matches=All_vs_All_without_FL_matches.sort_values(by='evalue', ascending=True)
All_vs_All_without_FL_matches = All_vs_All_without_FL_matches.drop_duplicates(subset = "query",keep="first")
All_vs_All_without_FL_matches=All_vs_All_without_FL_matches.loc[~All_vs_All_without_FL_matches['query'].isin(All_vs_All_bis['query'])]



#Merge all

All_vs_All_best_hits=All_vs_All_bis.append(All_vs_All_without_FL_matches)
All_vs_All_best_hits.loc[All_vs_All_best_hits['query'].str.contains("scaffold97118:4995-6326"),"query"]="scaffold97118:4995-7124(+):Trybliographa"
All_vs_All_best_hits.loc[All_vs_All_best_hits['query'].str.contains("scaffold97118:4995-712"),"pident"]=29.9
All_vs_All_best_hits.loc[All_vs_All_best_hits['query'].str.contains("scaffold97118:4995-712"),"qlen"]=710
All_vs_All_best_hits.loc[All_vs_All_best_hits['query'].str.contains("scaffold97118:4995-712"),"alnlen"]=688

All_vs_All_best_hits.loc[All_vs_All_best_hits['query'].str.contains("scaffold7282\\|size12268:1200-1859"),"query"]="scaffold7282|size12268:1200-2222(-):Rhoptromeris"
All_vs_All_best_hits.loc[All_vs_All_best_hits['query'].str.contains("scaffold7282\\|size12268:1200-2222"),"pident"]=42.8
All_vs_All_best_hits.loc[All_vs_All_best_hits['query'].str.contains("scaffold7282\\|size12268:1200-2222"),"qlen"]=341
All_vs_All_best_hits.loc[All_vs_All_best_hits['query'].str.contains("scaffold7282\\|size12268:1200-2222"),"alnlen"]=331

Blast_table= Blast_table.merge(All_vs_All_best_hits[['query','target','pident', 'alnlen', 'mismatch', 'gapopen', 'tstart', 'tend', 'qstart', 'qend', 'evalue', 'bits', 'qlen', 'tlen', 'tcov', 'qcov']],left_on="Names",right_on="query",how="left")



# Check 

# Count number filamentous in each clusters 
list_cynipoidea =['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV']
list_filamentous = ['LbFV','EfFV','PoFV','PcFV','CcFV1','CcFV2','LhFV']

count_cluster= Blast_table[['Cluster_hmmer','Names2']]
count_cluster = count_cluster.drop_duplicates(subset=['Cluster_hmmer', 'Names2'], keep='first')

count_cluster['Cynip_count'] = count_cluster['Names2'].isin(list_cynipoidea).groupby(count_cluster['Cluster_hmmer']).transform('sum')
count_cluster['Filamentous_count'] = count_cluster['Names2'].isin(list_filamentous).groupby(count_cluster['Cluster_hmmer']).transform('sum')
count_cluster = count_cluster.drop_duplicates(subset=['Cluster_hmmer'], keep='first')

Blast_table=Blast_table.merge(count_cluster[['Cluster_hmmer','Cynip_count','Filamentous_count']],on="Cluster_hmmer",how="left")


# Add gene name when available 
Gene_names=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Filamentous_gene_names.tab",sep="\t")
Gene_names['Names']=Gene_names['Names']+"_LbFV"
Gene_names_cluster=Blast_table.merge(Gene_names,on="Names",how="left")
Gene_names_cluster=Gene_names_cluster.loc[~Gene_names_cluster['Prot_name'].isna()]
Blast_table= Blast_table.drop_duplicates()
Blast_table=Blast_table.merge(Gene_names_cluster[['Cluster_hmmer','Prot_name']], on="Cluster_hmmer",how="left")


# Manually add three loci 


Manual_df={'Cluster_blast': {920: 'Cluster258'}, 'Names': {920: 'JADEYJ010000123.1:877336-878082(-):Leptopilina_boulardi_GCA_019393585.1'}, 'Names2': {920: 'LbEFV'}, 'Number_unique_Names2': {920: 12}, 'Nb_duplicated': {920: 0}, 'Ratio_duplicated_unique': {920: 0.0}, 'Cluster_hmmer': {920: 'Cluster258'}, 'Length': {920: 249.0}, 'query': {920: 'nan'}, 'target': {920: 'PHTE01000008.1:237233-237979(+):Leptopilina_boulardi_GCA_011634795.1'}, 'pident': {920: 96.8}, 'alnlen': {920: 249.0}, 'mismatch': {920: 8.0}, 'gapopen': {920: 0.0}, 'tstart': {920: 1.0}, 'tend': {920: 249.0}, 'qstart': {920: 1.0}, 'qend': {920: 249.0}, 'evalue': {920: 2.213e-165}, 'bits': {920: 510.0}, 'qlen': {920: 249.0}, 'tlen': {920: 249.0}, 'tcov': {920: 1.0}, 'qcov': {920: 1.0}, 'Cynip_count': {920: 6}, 'Filamentous_count': {920: 1}, 'Prot_name': {920: 'LbFVorf72'}}
Blast_table=pd.concat([Blast_table,pd.DataFrame.from_dict(Manual_df)])
#Manual_df={'Cluster_blast': {920: 'Cluster258'}, 'Names': {920: 'JABAIF010000692.1:288490-289236(+):Leptopilina_boulardi_GCA_015476485.1'}, 'Names2': {920: 'LbEFV'}, 'Number_unique_Names2': {920: 12}, 'Nb_duplicated': {920: 0}, 'Ratio_duplicated_unique': {920: 0.0}, 'Cluster_hmmer': {920: 'Cluster258'}, 'Length': {920: 249.0}, 'query': {920: 'nan'}, 'target': {920: 'PHTE01000008.1:237233-237979(+):Leptopilina_boulardi_GCA_011634795.1'}, 'pident': {920: 96.8}, 'alnlen': {920: 249.0}, 'mismatch': {920: 8.0}, 'gapopen': {920: 0.0}, 'tstart': {920: 1.0}, 'tend': {920: 249.0}, 'qstart': {920: 1.0}, 'qend': {920: 249.0}, 'evalue': {920: 2.213e-165}, 'bits': {920: 510.0}, 'qlen': {920: 249.0}, 'tlen': {920: 249.0}, 'tcov': {920: 1.0}, 'qcov': {920: 1.0}, 'Cynip_count': {920: 6}, 'Filamentous_count': {920: 1}, 'Prot_name': {920: 'LbFVorf72'}}
#Blast_table=pd.concat([Blast_table,pd.DataFrame.from_dict(Manual_df)])
Manual_df={'Cluster_blast': {923: 'Cluster258'}, 'Names': {923: 'NW_025111203.1:2150605-2150865(+):Leptopilina_heterotoma_GCA_015476425.1'}, 'Names2': {923: 'LhEFV'}, 'Number_unique_Names2': {923: 12}, 'Nb_duplicated': {923: 0}, 'Ratio_duplicated_unique': {923: 0.0}, 'Cluster_hmmer': {923: 'Cluster258'}, 'Length': {923: 249.0}, 'query': {923: "nan"}, 'target': {923: 'QYUC01227125.1:2219-2966(+):Leptopilina_heterotoma_GCA_009026005.1'}, 'pident': {923: 93.3}, 'alnlen': {923: 249.0}, 'mismatch': {923: 17.0}, 'gapopen': {923: 0.0}, 'tstart': {923: 1.0}, 'tend': {923: 249.0}, 'qstart': {923: 1.0}, 'qend': {923: 249.0}, 'evalue': {923: 4.001e-158}, 'bits': {923: 489.0}, 'qlen': {923: 249.0}, 'tlen': {923: 249.0}, 'tcov': {923: 1.0}, 'qcov': {923: 1.0}, 'Cynip_count': {923: 6}, 'Filamentous_count': {923: 1}, 'Prot_name': {923: 'LbFVorf72'}}
Blast_table=pd.concat([Blast_table,pd.DataFrame.from_dict(Manual_df)])
#Manual_df={'Cluster_blast': {323: 'Cluster126_redefined'}, 'Names': {323: 'JABAIF010007716.1:423933-424160(-):Leptopilina_boulardi_GCA_015476485.1'}, 'Names2': {323: 'LbEFV'}, 'Number_unique_Names2': {323: 47}, 'Nb_duplicated': {323: 22}, 'Ratio_duplicated_unique': {323: 0.4680851063829787}, 'Cluster_hmmer': {323: 'Cluster126_redefined'}, 'Length': {323: 76.0}, 'query': {323: 'JABAIF010007716.1:423933-424160(-):Leptopilina_boulardi_GCA_015476485.1'}, 'target': {323: 'B0YLN3_GpSGHV'}, 'pident': {323: 37.6}, 'alnlen': {323: 76.0}, 'mismatch': {323: 45.0}, 'gapopen': {323: 0.0}, 'tstart': {323: 559.0}, 'tend': {323: 631.0}, 'qstart': {323: 1.0}, 'qend': {323: 76.0}, 'evalue': {323: 1.205e-09}, 'bits': {323: 51.0}, 'qlen': {323: 76.0}, 'tlen': {323: 954.0}, 'tcov': {323: 0.077}, 'qcov': {323: 1.0}, 'Cynip_count': {323: 6}, 'Filamentous_count': {323: 7}, 'Prot_name': {323: 'LbFVorf58 (DNApol)'}}
#Blast_table=pd.concat([Blast_table,pd.DataFrame.from_dict(Manual_df)])

# Remove overlapping queries 
sub_Blast_table=Blast_table.loc[Blast_table['Names'].str.contains("\\(",na=False) & ~( Blast_table['Names'].str.contains("VcBV|PoFV"))]
#sub_Blast_table = sub_Blast_table.drop_duplicates()
sub_Blast_table['Scaffold_name']=sub_Blast_table['Names'].str.replace(":.*","")
sub_Blast_table['start-end']=sub_Blast_table['Names'].str.replace("\\(.*","")
sub_Blast_table['start-end']=sub_Blast_table['start-end'].str.replace(".*:","")
sub_Blast_table['start']=sub_Blast_table['start-end'].str.replace("-.*","")
sub_Blast_table['end']=sub_Blast_table['start-end'].str.replace(".*-","")
sub_Blast_table['end']=sub_Blast_table['end'].astype(int)
sub_Blast_table['start']=sub_Blast_table['start'].astype(int)

sub_Blast_table = sub_Blast_table.reset_index()

sub_Blast_table=sub_Blast_table.sort_values(['Cluster_hmmer','Scaffold_name', 'start'], ascending=[True, True,True])
is_overlapped = lambda x: x['start'] >= x['end'].shift(fill_value=-1)
sub_Blast_table['Overlapp_group'] = sub_Blast_table.sort_values(['Cluster_hmmer','Scaffold_name','start','end']) \
                .groupby(['Cluster_hmmer','Scaffold_name'],as_index=False).apply(is_overlapped).droplevel(0).cumsum()


sub_Blast_table['Tot_length'] = abs(sub_Blast_table['end'] - sub_Blast_table['start'])
sub_Blast_table['count_duplicate_overlapp'] = sub_Blast_table.groupby(['Overlapp_group'])['Overlapp_group'].transform('size')

#sub_Blast_table = sub_Blast_table.sort_values(['Cluster_hmmer','Scaffold_name', 'Tot_length'], ascending=[True, False,False]) \
#        .groupby(sub_Blast_table['Overlapp_group']).head(1)

# Remove still overlpping loci 

g = sub_Blast_table.groupby(['Cluster_hmmer','Scaffold_name'], group_keys=False)
group = sub_Blast_table['start'].gt(g['end'].apply(lambda s: s.shift().cummax())).cumsum()

sub_Blast_table['New_overlapping_group'] = (sub_Blast_table.groupby(['Cluster_hmmer','Scaffold_name', group])
                     .ngroup().add(1).astype(str)
                     .radd('G')
                  )
sub_Blast_table['len']=sub_Blast_table['end']-sub_Blast_table['start']
sub_Blast_table = sub_Blast_table.sort_values(['Cluster_hmmer','Scaffold_name','len'], ascending=[True, False, False]) \
        .groupby(sub_Blast_table['New_overlapping_group']).head(1)

Blast_table=Blast_table.loc[~(Blast_table['Names'].str.contains("\\(",na=False)) | (Blast_table['Names'].str.contains("VcBV|PoFV"))]

Blast_table=pd.concat([Blast_table,sub_Blast_table])
# Check candidate 

list_ORF_to_test=['LbFVorf72','LbFVorf83','LbFVorf108','LbFVorf87','LbFVorf94','LbFVorf5','LbFVorf60 (lcat)','LbFVorf107 (lef-4)','LbFVorf96 (lef-8)','LbFVorf78 (lef-9)','LbFVorf58 (DNApol)','LbFVorf68 (helicase)','LbFVorf92','LbFVorf85 (Ac81)']

for orf in list_ORF_to_test:
  print(orf)
  print(len(Blast_table.loc[Blast_table['Prot_name'].eq(orf) & Blast_table['Names2'].isin(['LhEFV','LbEFV','LcEFV','ThEFV','TrEFV','RhEFV'])]['Names'].unique()))
  print("\n")


Blast_table.to_csv(output,sep=";",index=False)
print("File saved to : ",output)

"""
# Check pseudogenes 
AA_records=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses.faa","fasta"))

for index, row in Blast_table.iterrows():
	if '*' in AA_records[row['Names']].seq:
		if row['Prot_name'] in  list_ORF_to_test:
			if row['Names2'] in list_cynipoidea:
				print(row['Prot_name'])
				print('>', row['Names'],sep="")
				print(AA_records[row['Names']].seq)
				print("\n")
		print()
"""

