import pandas as pd 
import sys 
import argparse
import os


# Print out a message when the program is initiated.
print('------------------------------------------------------------------------------\n')
print('                        Gather HMMEr files and add multiple informations     .\n')
print('------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Gather HMMEr files')
parser.add_argument("-p", "--hmmer_path", help="The path with all hmmer files")
parser.add_argument("-b", "--blast_cluster_file", help="The blast cluster file")
parser.add_argument("-out1", "--out1", help="The first output containing clusters informations")
parser.add_argument("-out2", "--out2", help="The first output containing clusters  and coverages informations ")
parser.add_argument("-out3", "--out3", help="The first output containing clusters  and coverages and NR informations ")
parser.add_argument("-out4", "--out4", help="The first output containing clusters  and coverages and NR and ORFs informations ")
args = parser.parse_args()

#Example usage : 

"""
#
python3 Gather_hmmer_clusters_and_add_informations.py -p /beegfs/data/bguinet/Cynipoidea_project/Clustering/Cluster_hmmer -b /beegfs/data/bguinet/Cynipoidea_project/Clustering/ALL_Predicted_and_known_ORFs_cluster.tab \
-out1 dsDNA_hmmer_clusters.tab -out2 dsDNA_hmmer_clusters_cov.tab -out3 dsDNA_hmmer_clusters_cov_NR.tab -out4 dsDNA_hmmer_clusters_cov_NR_ORF.tab


"""
# Variable that stores fasta sequences

path=args.hmmer_path
blast=args.blast_cluster_file
out1=args
path="/beegfs/data/bguinet/Cynipoidea_project/Clustering/Cluster_hmmer"
blast="/beegfs/data/bguinet/Cynipoidea_project/Clustering/ALL_Predicted_and_known_ORFs_cluster.tab"
out1="/beegfs/data/bguinet/Cynipoidea_project/Clustering/dsDNA_hmmer_clusters.tab"
out2="/beegfs/data/bguinet/Cynipoidea_project/Clustering/dsDNA_hmmer_clusters_cov.tab"
out3="/beegfs/data/bguinet/Cynipoidea_project/Clustering/dsDNA_hmmer_clusters_cov_NR.tab"
out4="/beegfs/data/bguinet/Cynipoidea_project/Clustering/dsDNA_hmmer_clusters_cov_NR_ORF.tab"

from collections import defaultdict
import pandas as pd
from Bio import SearchIO
import re 
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
      Main_Hmmer_tab=Main_Hmmer_tab.append(hmmer_tab)
      #print(hmmer_tab)


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


Main_Hmmer_tab=Main_Hmmer_tab.loc[Main_Hmmer_tab['evalue'].lt(9.000000e-05)] #since in control the lowest = 7.300000e-05


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
Blast_table_duplicated_clusters=Blast_table2.loc[Blast_table2['Ratio_duplicated_unique'].gt(0.60) & Blast_table2['Nb_duplicated'].ge(4)]

# Re-define more stringent clusters for those sequences 

Main_Hmmer_tab_save2=Main_Hmmer_tab_save2.merge(Blast_table_duplicated_clusters[['Cluster_blast', 'Names']],right_on="Names",left_on="id",how="right")
Main_Hmmer_tab_save2['Cluster']=Main_Hmmer_tab_save2['Cluster'].str.replace(".tab","")
Main_Hmmer_tab_save2=Main_Hmmer_tab_save2.loc[ ~ (Main_Hmmer_tab_save2['Cluster'] == Main_Hmmer_tab_save2['Cluster_blast'])]

Main_Hmmer_tab_save2=Main_Hmmer_tab_save2.loc[Main_Hmmer_tab_save2['bitscore'].gt(100)] 
Main_Hmmer_tab_save2=Main_Hmmer_tab_save2.loc[~(Main_Hmmer_tab_save2['Cluster_blast'].str.contains("Cluster") & Main_Hmmer_tab_save2['Cluster'].isna())]

if len(Main_Hmmer_tab_save2)==0:
  Blast_table_duplicated_clusters['Cluster_hmmer']  = Blast_table_duplicated_clusters['Cluster_blast']
else:
  G2 = nx.from_pandas_edgelist(Main_Hmmer_tab_save2, 'Cluster', 'Cluster_blast')
  d2 = {k: '-'.join(c) for c in nx.connected_components(G2) for k in c}
  Blast_table_duplicated_clusters['Cluster_hmmer'] = Blast_table_duplicated_clusters['Cluster_blast'].replace(d2)

# Add new result to previous Blast file 

Blast_table=Blast_table[~Blast_table['Cluster_blast'].isin(Blast_table_duplicated_clusters['Cluster_blast'])]
Blast_table=Blast_table.append(Blast_table_duplicated_clusters)



blast_file="/beegfs/data/bguinet/Cynipoidea_project/Clustering/ALL_vs_ALL_Predicted_and_known_ORFs_and_EVEs_and_new_result.m8"

df=pd.read_csv(blast_file,sep="\t",header=None)
df.columns=['query','target','pident','alnlen','mismatch','gapopen','tstart','tend','qstart','qend','evalue','bits','tlen','qlen','qcov','tcov']
df2=df[['target','query','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','qlen','tlen','qcov','tcov']]
df2.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']
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

Blast_table_duplicated_clusters=Blast_table2.loc[Blast_table2['Ratio_duplicated_unique'].gt(0.60) & Blast_table2['Nb_duplicated'].ge(4)]

# Re-define more stringent clusters for those sequences 

Main_Hmmer_tab_save=Main_Hmmer_tab_save.merge(Blast_table_duplicated_clusters[['Cluster_blast', 'Names']],right_on="Names",left_on="id",how="right")
Main_Hmmer_tab_save['Cluster']=Main_Hmmer_tab_save['Cluster'].str.replace(".tab","")
Main_Hmmer_tab_save=Main_Hmmer_tab_save.loc[ ~ (Main_Hmmer_tab_save['Cluster'] == Main_Hmmer_tab_save['Cluster_blast'])]

Main_Hmmer_tab_save=Main_Hmmer_tab_save.loc[Main_Hmmer_tab_save['bitscore'].gt(300)] 
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



#Manuall stuff 

# Manually add the ORF72 LbFV EVE within the long read Lboulardi assembly : >JADEYJ010000123.1:877336-878082(-):Leptopilina_boulardi_GCA_019393585.1 and put the same information as in JABAIF010000692.1:288490-289236(+):Leptopilina_boulardi_GCA_015476485
Cluster_Hmmer_to_asign=Blast_table.loc[Blast_table['Names'].str.contains("JABAIF010000692.1:288490-289236")]['Cluster_hmmer'].iloc[0]
Cluster_Blast_to_asign=Blast_table.loc[Blast_table['Names'].str.contains("JABAIF010000692.1:288490-289236")]['Cluster_blast'].iloc[0]
Length_to_assign=249

Blast_table=Blast_table.append({'Cluster_blast':Cluster_Blast_to_asign, 'Names':"JADEYJ010000123.1:877336-878082(-):Leptopilina_boulardi_GCA_019393585.1", 'Names2':"LbEFV", 'Number_unique_Names2':12, 'Nb_duplicated':0, 'Ratio_duplicated_unique':0.0, 'Cluster_hmmer':Cluster_Hmmer_to_asign, 'Length':Length_to_assign}, ignore_index=True)

# Manually add the ORF72 LbFV EVE within the long read  Lheterotoma assembly : NW_025111203.1:2150605-2150865(+):Leptopilina_heterotoma_GCA_015476425.1 and put the same information as in   QYUB01233969.1:605-1352(-):Leptopilina_heterotoma_GCA_009025955
Cluster_Hmmer_to_asign=Blast_table.loc[Blast_table['Names'].str.contains("QYUB01233969.1:605-1352")]['Cluster_hmmer'].iloc[0]
Cluster_Blast_to_asign=Blast_table.loc[Blast_table['Names'].str.contains("QYUB01233969.1:605-1352")]['Cluster_blast'].iloc[0]
Length_to_assign=87

Blast_table=Blast_table.append({'Cluster_blast':Cluster_Blast_to_asign, 'Names':"NW_025111203.1:2150605-2150865(+):Leptopilina_heterotoma_GCA_015476425.1", 'Names2':"LhEFV", 'Number_unique_Names2':12, 'Nb_duplicated':0, 'Ratio_duplicated_unique':0.0, 'Cluster_hmmer':Cluster_Hmmer_to_asign, 'Length':Length_to_assign}, ignore_index=True)

# Manually correct two ORFs caused by  sequencing mistakes 
#Trybliographa

Blast_table.loc[Blast_table['Names'].str.contains("scaffold97118:4995-6326"),"Length"]=710
Blast_table.loc[Blast_table['Names'].str.contains("scaffold97118:4995-6326"),"Names"]="scaffold97118:4995-7124(+):Trybliographa"
Blast_table=Blast_table.loc[~Blast_table['Names'].str.contains("scaffold97118:6330-7124")]

#Rhoptromeris 
Blast_table.loc[Blast_table['Names'].str.contains("scaffold7282\\|size12268:1200-1859"),"Length"]=341
Blast_table.loc[Blast_table['Names'].str.contains("scaffold7282\\|size12268:1200-1859"),"Names"]="scaffold7282|size12268:1200-2222(-):Rhoptromeris"
Blast_table=Blast_table.loc[~Blast_table['Names'].str.contains("scaffold7282\\|size12268:1878-2222")]



# Add best hits informations for each queries 
#Revious blast all vs all file 

All_vs_All=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Clustering/ALL_vs_ALL_Predicted_and_known_ORFs_and_EVEs_and_new_result.m8",sep="\t",header=None)
All_vs_All.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']


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




Blast_table.to_csv(out1,sep=";",index=False)
print("File saved to : ",out1)

#for ids in ['AQQ79998','AQQ79975','AQQ80022','AQQ79957','AQQ79983','AQQ79939','AQQ80012','AQQ79978','AQQ79925','AQQ79952','AQQ80016','AQQ80016','AQQ80007','AQQ80027',
#'AQQ79972','AQQ79980','AQQ79988','AQQ80001','AQQ79922','AQQ79940','AQQ79963','AQQ79987','AQQ79994','AQQ80005','AQQ80019','AQQ80017','AQQ80003','AQQ80014','AQQ79992']:
#	Blast_table.loc[Blast_table['Cluster_hmmer'].str.contains(Blast_table.loc[Blast_table['Names'].str.contains(ids)]['Cluster_hmmer'].iloc[0])].iloc[0:2]
#                          EVE   FL                             core 
# AQQ79998 > lef-9        6/6      7/7                         yes
# AQQ79975 > pif-3        0/6      7/7 
# AQQ80022 > LbFVorf102   0/6      7/7 
# AQQ79957 > Ac38         0/6      7/7 
# AQQ79983 > p74          0/6      7/7
# AQQ79939 > 38k          1/6      7/7  (RhEFV)
# AQQ80012 > LbFVorf92    6/6      7/7                         yes
# AQQ79978 > DNAPOL       6/6      7/7                         yes
# AQQ79925 > LbFVorf5     6/6      7/7                         yes
# AQQ79952 > PIF1         0/6      7/7
# AQQ80016 > LEF8         6/6      7/7                         yes
# AQQ80007 > LbFVorf87    6/6      7/7                         yes
# AQQ80027 > LEF4         6/6      7/7                         yes
# AQQ79972 > PIF2         1/6      7/7 (RhEFV)
# AQQ79980 > LCAT         6/6      7/7                         yes
# AQQ79988 > Helicase     6/6      7/7                         yes
# AQQ80001 > ATPase       6/6      6/7  (PoFV)                 yes
# AQQ79922 > Integrase    6/6      7/7                         yes
# AQQ79940 > LbFVorf20    0/6      7/7
# AQQ79963 > PIF5         0/6      7/7
# AQQ79987 > PDDEXK       0/6      7/7
# AQQ79994 > p69          0/6      1/7         
# AQQ80005 > Ac81         6/6      7/7                          yes
# AQQ80019 > LbFVorf99    0/6      7/7
# AQQ80017 > P33          0/6      7/7
# AQQ79992 > LbFVorf72   
#Pas coeur 
# AQQ80003 > LbFVorf83    6/6      4/7  (PcFV,CcFV2, CcFV1)     yes
# AQQ80014 > LbFVorf94    6/6      7/7   (creuser)              yes
# AQQ79992 > LbFVorf72    6/6      1/7                          yes



## Add coverage information 

Blast_table = pd.read_csv(out1,sep=";")

import subprocess
from statistics import mean

Blast_table['Scaffold_name']=Blast_table['Names'].str.replace(":.*","")
Blast_table['Species_name']=Blast_table['Names'].str.replace(".*:","")

list_cov_scaff_viral=[]
for sp in list_cynipoidea:
	print("\n")
	print("Processing ", sp)
	sub_Blast_table=Blast_table.loc[Blast_table['Names2'].str.contains(sp)]
	species=sub_Blast_table['Names'].str.replace(".*:","").iloc[0]
	print(species)
	#In order to open the file correctly
	subprocess.call(" sed -i 's@# Busco@Busco@g' /beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/assembly/BUSCO_results_hymenoptera/run_hymenoptera_odb10/full_table.tsv", shell=True)
	if sp == "LhFV":
		BUSCO_table=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/assembly/BUSCO_results_hymenoptera/run_hymenoptera_odb10/full_table.tsv",comment="#",sep="\t")
		heterotoma_tab=pd.DataFrame(columns=['Scaffold_name','Info'])
		for scf in SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/Leptopilina_heterotoma_GCA_015476425.1/assembly/Leptopilina_heterotoma_GCA_015476425.1_final_assembly.fna", "fasta"):
			heterotoma_tab=heterotoma_tab.reset_index(drop=True)
			heterotoma_tab = heterotoma_tab.append({'Scaffold_name':scf.id, 'Info':scf.description},ignore_index=True)
		BUSCO_table_bis=BUSCO_table.loc[BUSCO_table['Sequence'].isin(heterotoma_tab.loc[heterotoma_tab['Info'].str.contains("whole genome shotgun sequence")]['Scaffold_name'])]
		BUSCO_table=BUSCO_table.loc[~BUSCO_table['Species_name'].str.contains("hetero")]
		BUSCO_table=BUSCO_table.append(BUSCO_table_bis)
		Blast_table_bis=Blast_table.loc[Blast_table['Scaffold_name'].isin(heterotoma_tab.loc[heterotoma_tab['Info'].str.contains("whole genome shotgun sequence")]['Scaffold_name']) & Blast_table['Species_name'].str.contains("hetero")]
		Blast_table=Blast_table.loc[~Blast_table['Species_name'].str.contains("hetero")]
		Blast_table=Blast_table.append(Blast_table_bis)
	else:
		BUSCO_table=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/assembly/BUSCO_results_hymenoptera/run_hymenoptera_odb10/full_table.tsv",comment="#",sep="\t")
	#Keep only non-missing busco scaffolds
	BUSCO_table=BUSCO_table.loc[~(BUSCO_table['Status'].eq('Missing'))]
	if os.path.exists("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/assembly/Mapping/"+species+"_mapping_sorted_reduced.coverage") :
		Coverage_tab_reduced=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/assembly/Mapping/"+species+"_mapping_sorted_reduced.coverage",sep="\t")
	else:
		Coverage_tab=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/assembly/Mapping/"+species+"_mapping_sorted.coverage",sep="\t",header=None)
		Coverage_tab.columns=['Scaffold_name','position','coverage']
		#Calculate median coverage and mean coverage for each scaffolds 
		Coverage_tab2=Coverage_tab.groupby("Scaffold_name")["coverage"].mean()
		Coverage_tab3=Coverage_tab.groupby("Scaffold_name")["coverage"].median()
		Coverage_tab2=Coverage_tab2.to_frame()
		Coverage_tab2.columns=['Mean_cov']
		Coverage_tab3=Coverage_tab3.to_frame()
		Coverage_tab3.columns=['Median_cov']
		Coverage_tab_reduced=pd.merge(Coverage_tab3,Coverage_tab2,on="Scaffold_name")
		# Add GC content
		from Bio.Seq import Seq
		from Bio.SeqRecord import SeqRecord
		from Bio.SeqUtils import GC
		Assembly=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/assembly/"+species+"_final_assembly.fna", "fasta"))
		Coverage_tab_reduced2=pd.json_normalize([{'Scaffold_name':k,
			'Scaffold_length': len(s.seq),
			'GC': GC(s.seq)}
			for k, s in Assembly.items()])
		# Save the coverage 
		Coverage_tab_reduced=Coverage_tab_reduced.merge(Coverage_tab_reduced2,on="Scaffold_name")
		Coverage_tab_reduced.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/assembly/Mapping/"+species+"_mapping_sorted_reduced.coverage",sep="\t",index=False)
	# Write BUSCO file with all informations for plot 
	if sp =="TrEFV":
		BUSCO_table['Sequence']=BUSCO_table['Sequence'].str.replace("\\|.*","")
	BUSCO_table=BUSCO_table.merge(Coverage_tab_reduced,left_on="Sequence",right_on="Scaffold_name",how="left")
	Mean_busco_GC=mean(BUSCO_table['GC'])
	Mean_busco_cov=mean(BUSCO_table['Mean_cov'])
	Median_busco_cov=mean(BUSCO_table['Median_cov'])
	BUSCO_table.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/assembly/BUSCO_results_hymenoptera/run_hymenoptera_odb10/full_table_cov_gc.tsv",sep=";",index=False)
	nb_scaff=len(BUSCO_table['Sequence'].unique())
	n=1
	BUSCO_cov_distribution_median=[]
	BUSCO_cov_distribution_mean=[]
	for busco_scaff in BUSCO_table['Sequence'].unique():
		try:
			if sp == "TrEFV":
				busco_scaff=re.sub("\\|.*","",busco_scaff)
			Median_cov=Coverage_tab_reduced.loc[Coverage_tab_reduced['Scaffold_name'].eq(busco_scaff)]['Median_cov'].iloc[0]
			BUSCO_cov_distribution_median.append(Median_cov)
			Mean_cov=Coverage_tab_reduced.loc[Coverage_tab_reduced['Scaffold_name'].eq(busco_scaff)]['Mean_cov'].iloc[0]
			BUSCO_cov_distribution_mean.append(Mean_cov)
			#print(n,"/",nb_scaff)
			n+=1
		except:
			print(busco_scaff, " did not have coverage")
	print("Running cov dev pvalue assessment .....")
	sub_Blast_table=sub_Blast_table.merge(Coverage_tab_reduced,on="Scaffold_name",how="left")
	for index, row in sub_Blast_table.iterrows():
		if row['Median_cov'] > Mean_busco_cov:
			pvalue_median=(sum(i > row['Median_cov'] for i in BUSCO_cov_distribution_median)/len(BUSCO_cov_distribution_median))
		else:
			pvalue_median=(sum(i < row['Median_cov'] for i in BUSCO_cov_distribution_median)/len(BUSCO_cov_distribution_median))
		if row['Mean_cov'] > Mean_busco_cov:
                        pvalue_mean=(sum(i > row['Mean_cov'] for i in BUSCO_cov_distribution_mean)/len(BUSCO_cov_distribution_mean))
		else:
			pvalue_mean=(sum(i < row['Mean_cov'] for i in BUSCO_cov_distribution_mean)/len(BUSCO_cov_distribution_mean))
		list_cov_scaff_viral.append({'Species_name':species,'Scaffold_name':row['Scaffold_name'],'Scaffold_length':row['Scaffold_length'],"Mean_GC_candidat":row['GC'], "Mean_GC_BUSCO":Mean_busco_GC,'Mean_cov_depth_candidat':row['Mean_cov'],'Median_cov_depth_candidat':row['Median_cov'],'Mean_cov_depth_BUSCO':Mean_busco_cov,'Median_cov_depth_BUSCO':Median_busco_cov,'pvalue_cov_mean':pvalue_mean, 'pvalue_cov_median':pvalue_median,})

candidates_cov_scaff=pd.DataFrame(list_cov_scaff_viral)


# FDR analysis

from multipy.fdr import tst

list_nb=[0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01]

for i in list_nb:
       candidates_cov_scaff['FDR_pvalue_cov_mean']=tst(candidates_cov_scaff['pvalue_cov_mean'], q=i)
       print(i,len(candidates_cov_scaff.loc[candidates_cov_scaff['FDR_pvalue_cov_mean'].astype(str)=="False"]))

#Here we take 0.04

del candidates_cov_scaff['FDR_pvalue_cov_mean']

FDR_with_cov_inf= candidates_cov_scaff.drop_duplicates(subset=['Scaffold_name'], keep='first')
FDR_with_cov_inf['FDR_pvalue_cov_mean']=tst(FDR_with_cov_inf['pvalue_cov_mean'], q=0.04)
FDR_with_cov_inf['FDR_pvalue_cov_median']=tst(FDR_with_cov_inf['pvalue_cov_median'], q=0.04)
FDR_with_cov_inf=FDR_with_cov_inf[['Scaffold_name','Species_name','FDR_pvalue_cov_mean','FDR_pvalue_cov_median']]

Blast_table=Blast_table.merge(FDR_with_cov_inf,on=['Scaffold_name','Species_name'],how="left")

Blast_table = Blast_table.merge(candidates_cov_scaff, on=["Scaffold_name","Species_name"],how="left")	

Blast_table = Blast_table.drop_duplicates()

Blast_table.to_csv(out2,sep=";",index=False)
print("File saved to : ",out2)

# Save the file for plotting heatmap gene content 


# Add NR results 


## NR blast file 

Nr_blast="/beegfs/data/bguinet/Cynipoidea_project/Clustering/NR_homology_ALL_Predicted_and_known_ORFs_and_EVEs_and_new_result.m8"

NR_blast=pd.read_csv(Nr_blast,sep="\t",header=None)
NR_blast.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qcov','tcov','taxid','taxname','taxlineage']

NR_blast.loc[NR_blast['query'].str.contains("scaffold97118:4995-6326"),"query"]="scaffold97118:4995-7124(+):Trybliographa"
NR_blast.loc[NR_blast['query'].str.contains("scaffold7282\\|size12268:1200-1859"),"Names"]="scaffold7282|size12268:1200-2222(-):Rhoptromeris"


#Remove ichnovirus and wrongly annotated nudiviruses 
NR_blast=NR_blast.loc[~NR_blast['taxname'].str.contains("ichno")]

NR_blast=NR_blast.loc[~NR_blast['taxname'].str.contains("Cotesia|Microplitis|Chelonus")]


#Sort evalues 
NR_blast=NR_blast.sort_values(by='evalue', ascending=True)

# We keep only one species for each query 
NR_blast = NR_blast.drop_duplicates(subset=['query', 'taxname'], keep='first')

NR_blast=NR_blast.loc[NR_blast['bits'].ge(50)]

###


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
NR_blast.loc[ NR_blast['count_Virus'].eq(0)  & NR_blast['Count_Bacteria_Eukaryota_Archeae'].eq(0) ,"Viral_confidence"] = "Viral"
NR_blast.loc[NR_blast['Viral_confidence'].str.contains("NA") ,"Viral_confidence"] = "Uncertain"


#First NR best hit 
NR_blast=NR_blast.sort_values(by='evalue', ascending=True)
Best_NR_blast= NR_blast.drop_duplicates(subset=['query'], keep='first')

Best_NR_blast['Best_NR']="Unknown"
for index, row in Best_NR_blast.iterrows():
	if "acter" in row['taxlineage']:
		Best_NR_blast.loc[Best_NR_blast['query'].eq(row['query']),'Best_NR']="Bacteria"
	elif "irus" in row['taxlineage']:
		Best_NR_blast.loc[Best_NR_blast['query'].eq(row['query']),'Best_NR']="Virus"
	elif "rchaea" in row['taxlineage']:
		Best_NR_blast.loc[Best_NR_blast['query'].eq(row['query']),'Best_NR']="Archaea"
	elif "ukaryota" in row['taxlineage']:
                Best_NR_blast.loc[Best_NR_blast['query'].eq(row['query']),'Best_NR']="Eukaryota"

Blast_table=pd.read_csv(out2,sep=";")
Blast_table=Blast_table.merge(Best_NR_blast[['query','Best_NR']],right_on="query",left_on="Names",how="left")

Blast_table=Blast_table.merge(NR_blast[['query','count_Archaea','count_Eukaryota','count_Bacteria','count_Virus','Count_Bacteria_Eukaryota_Archeae','Perc_viral_vs_other','Viral_confidence']],left_on="Names",right_on="query",how="left")
Blast_table=Blast_table.drop_duplicates(subset = "Names")
Blast_table.loc[Blast_table['Viral_confidence'].isna() ,"Viral_confidence"] = "Viral"


Blast_table['Cluster_count_Viral'] = Blast_table['Viral_confidence'].str.contains('Viral').groupby(Blast_table['Cluster_hmmer']).transform('sum')
Blast_table['Cluster_count_Uncertain'] = Blast_table['Viral_confidence'].str.contains('Uncertain').groupby(Blast_table['Cluster_hmmer']).transform('sum')
Blast_table['Cluster_count_Not_Viral'] = Blast_table['Viral_confidence'].str.contains('Not_viral').groupby(Blast_table['Cluster_hmmer']).transform('sum')

Blast_table['Cluster_Perc_viral_vs_other'] = Blast_table['Cluster_count_Viral'] / ( Blast_table['Cluster_count_Not_Viral'])


# Add protein names
Table_names=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Clustering/Filamentous_core_genes.tab",sep="\t")

#del Blast_table['Prot_name']

LbFV_Blast_table=Blast_table.loc[Blast_table['Names2'].str.contains("LbFV")]
LbFV_Blast_table['Names3']=LbFV_Blast_table['Names'].str.replace("_.*","")

LbFV_Blast_table=LbFV_Blast_table.merge(Table_names,right_on="Names",left_on="Names3", how="left")


Blast_table=Blast_table.merge(LbFV_Blast_table[['Cluster_hmmer','Prot_name']],on="Cluster_hmmer")

#del Blast_table['Names3']
Blast_table.to_csv(out3,sep=";",index=False)

print("File saved to : ",out3)


##############
## Run ORF analysis 

# Extract the scaffolds of interest 
from Bio import SeqIO
import pandas as pd 
Blast_table=pd.read_csv(out3,sep=";")

Blast_table['Species_name'].mask(Blast_table['Species_name'].str.contains("PoEFV") ,'PoEFV',inplace=True)
Blast_table['Scaffold_name'].mask(Blast_table['Species_name'].str.contains("PoEFV") , Blast_table['Scaffold_name'].str.replace("-.*",""),inplace=True)
Blast_table['Scaffold_name'].mask(Blast_table['Species_name'].str.contains("PoEFV") , Blast_table['Scaffold_name'].str.split('_').str[:-1].str.join('_'),inplace=True)


list_cynipoidea =['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV']

with open("/beegfs/data/bguinet/Cynipoidea_project/All_candidate_scaffold.fna","w") as output:
	for sp in Blast_table.loc[Blast_table['Names2'].isin(['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV','PoEFV'])]['Species_name'].unique():
		subTable=Blast_table.loc[Blast_table['Species_name'].eq(sp)]
		if sp =="PoEFV":
			record_dict=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/Platygaster_orseoliae_corrected2.fa","fasta"))
		else:
			record_dict=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+sp+"/assembly/"+sp+"_final_assembly.fna","fasta"))
		for scaffold in subTable['Scaffold_name'].unique():
			print('>',sp,":",record_dict[scaffold].id,sep="",file=output)
			print(record_dict[scaffold].seq,file=output)

# Run python ORF finder 
import subprocess
import pathlib  


subprocess.run("/beegfs/data/bguinet/Bguinet_conda/bin/orfipy /beegfs/data/bguinet/Cynipoidea_project/All_candidate_scaffold.fna --ignore-case  --procs 5 --outdir /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis --bed /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy.bed --min 150 --start ATG --dna /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy.fna", shell=True)
#/beegfs/data/bguinet/Bguinet_conda/bin/orfipy /beegfs/data/bguinet/Cynipoidea_project/All_candidate_scaffold.fna --procs 5 --outdir /beegfs/data/bguinet/Cynipoidea_project/ --bed /beegfs/data/bguinet/Cynipoidea_project/All_candidate_scaffold_orfipy.bed --min 150 --start ATG --dna /beegfs/data/bguinet/Cynipoidea_project/All_candidate_scaffold_orfipy.fna
# Translate 

ORF_bed=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy.bed",sep="\t",header=None)
ORF_bed.columns=['Scaffold_name','ORF_start','ORF_end','ORF_name','zero','ORF_strand']

ORF_bed['ORF_name2']=ORF_bed['ORF_name'].str.replace(";.*","")
ORF_bed['ORF_name2']=ORF_bed['ORF_name2'].str.replace("ID=","")

with open("/beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy.faa","w") as output:
	record_dict=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy.fna","fasta"))
	for species in ORF_bed['Scaffold_name'].unique():
		for index, row in ORF_bed.loc[ORF_bed['Scaffold_name'].str.contains(species)].iterrows():
			print(">",row['ORF_name2'],';',row['ORF_start'],"-",row['ORF_end'],"(",row['ORF_strand'],")",sep="",file=output)
			print(record_dict[row['ORF_name2']].seq.translate(),file=output)

# Run mmseqs between ORFs and previously selected candidates

record_dict=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Clustering/ALL_Predicted_and_known_ORFs_and_EVEs_and_new.faa","fasta"))

Blast_table['Names3']=Blast_table['Names'].str.replace("\\(.*","")
Blast_table['end']=Blast_table['Names3'].str.replace(".*-","")
Blast_table['start']=Blast_table['Names3'].str.replace(".*:","")
Blast_table['start']=Blast_table['start'].str.replace("-.*","")
Blast_table['strand']=Blast_table['Names'].str.replace("\\).*","")
Blast_table['strand']=Blast_table['strand'].str.replace(".*\\(","")
#Deal with PoEFV 


Blast_table.loc[Blast_table['Names'].str.contains("_\\+_"),"strand"]="+"
Blast_table.loc[Blast_table['Names'].str.contains("_-_"),"strand"]="-"
Blast_table['start-end']=Blast_table['Names'].str.replace("_\\+.*","")
Blast_table['start-end']=Blast_table['start-end'].str.replace("_-.*","")
Blast_table['start-end']=Blast_table['start-end'].str.replace(".*_","")
Blast_table['end'].mask(Blast_table['Species_name'].str.contains("PoEFV"),Blast_table['start-end'].str.replace(".*-",""),inplace=True)
Blast_table['start'].mask(Blast_table['Species_name'].str.contains("PoEFV"),Blast_table['start-end'].str.replace("-.*",""),inplace=True)
del Blast_table['start-end']




list_cynipoidea =['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV','PoEFV']

with open ("/beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_EVEs.faa","w") as output:
	for index, row in Blast_table.loc[Blast_table['Names2'].isin(list_cynipoidea)].iterrows():
		if row['Species_name']=="PoEFV":
			seq_name=row['Scaffold_name']+'_'+str(int(row['start']))+'-'+str(int(row['end']))+'_'+row['strand']+'_'+row['Species_name']
		else:
			seq_name=row['Scaffold_name']+':'+str(int(row['start']))+'-'+str(int(row['end']))+'('+row['strand']+'):'+row['Species_name']
		print('>',record_dict[seq_name].id,sep="",file=output)
		print(record_dict[seq_name].seq,file=output)



subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_EVEs.faa /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_EVEs_db", shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy.faa /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy_db", shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs search /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_EVEs_db /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy_db /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy_EVEs_results /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy_EVEs_tpm --threads 10", shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs convertalis --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tlen,qcov,tcov' /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_EVEs_db /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy_db /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy_EVEs_results /beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy_EVEs_result.m8 --threads 10", shell=True)

ORF_vs_EVE_table=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/ORF_analysis/All_candidate_scaffold_orfipy_EVEs_result.m8",sep="\t",header=None)
ORF_vs_EVE_table.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qcov','tcov']
# Keep only self matching hits 
ORF_vs_EVE_table=ORF_vs_EVE_table.loc[ORF_vs_EVE_table['bits'].ge(50)]

ORF_vs_EVE_table['query_EVE_species']=ORF_vs_EVE_table['query'].str.replace(".*:","")
ORF_vs_EVE_table['query_EVE_species'].mask(ORF_vs_EVE_table['query'].str.contains("PoEFV"),ORF_vs_EVE_table['query'].str.replace(".*_",""),inplace=True)
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['target'].str.replace(".*;","")
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['ORF_start'].str.replace("\\(\\+\\)","")
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['ORF_start'].str.replace("\\(-\\)","")
ORF_vs_EVE_table['ORF_strand']="NA"
ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query'].str.contains("\\+"),'ORF_strand'] ="+"
ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query'].str.contains("-"),'ORF_strand'] ="-"
#ORF_vs_EVE_table['target'].str.replace(".*;","")
ORF_vs_EVE_table['ORF_end']=ORF_vs_EVE_table['ORF_start'].str.replace(".*-","").astype(int)
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['ORF_start'].str.replace("-.*","").astype(int)
ORF_vs_EVE_table['target_ORF_species']=ORF_vs_EVE_table['target'].str.replace(":.*","")
ORF_vs_EVE_table['query_EVE_scaffold']=ORF_vs_EVE_table['query'].str.replace(":.*","")
ORF_vs_EVE_table['query_EVE_scaffold'].mask(ORF_vs_EVE_table['query'].str.contains("PoEFV") , ORF_vs_EVE_table['query'].str.replace("-.*",""),inplace=True)
ORF_vs_EVE_table['query_EVE_scaffold'].mask(ORF_vs_EVE_table['query'].str.contains("PoEFV") , ORF_vs_EVE_table['query_EVE_scaffold'].str.split('_').str[:-1].str.join('_'),inplace=True)
ORF_vs_EVE_table['target_ORF_scaffold']=ORF_vs_EVE_table['target'].str.replace(".*:","")
ORF_vs_EVE_table['target_ORF_scaffold']=ORF_vs_EVE_table['target_ORF_scaffold'].str.replace("_ORF.*","")

ORF_vs_EVE_table=ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query_EVE_species'] == ORF_vs_EVE_table['target_ORF_species']]
ORF_vs_EVE_table=ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query_EVE_scaffold'] == ORF_vs_EVE_table['target_ORF_scaffold']]

# Remove duplicate 
ORF_vs_EVE_table = ORF_vs_EVE_table.drop_duplicates()

# find overlapping ORFs within the same scaffold
ORF_vs_EVE_table.rename(columns={'query': 'ORF_query',
                   'target': 'ORF_target'},
          inplace=True, errors='raise')

Table_with_ORFs=Blast_table.merge(ORF_vs_EVE_table[['ORF_query','ORF_target','ORF_start','ORF_end','ORF_strand']],left_on="Names",right_on="ORF_query",how="left")


#Table_with_ORFs['start']= Table_with_ORFs['start'].astype(int)
#Table_with_ORFs['end']= Table_with_ORFs['end'].astype(int)

import numpy as np

Table_with_ORFs['Overlapp_ORF_EVEs']= np.nan
Table_with_ORFs['ORF_perc']= np.nan
for index, row in Table_with_ORFs.loc[Table_with_ORFs['Names2'].isin(list_cynipoidea) & Table_with_ORFs['ORF_start'].ge(0)].iterrows():
	overlapp="no"
	if (int(row['start']) <= int(row['ORF_start'])) & (int(row['end']) <= int(row['ORF_end'])) & (int(row['end']) > int(row['ORF_start'])):
		#print("left")
		#print(row['start']," : ", row['end'])
		#print(row['ORF_start']," : ", row['ORF_end'])
		overlapp = "yes-left"
	elif (int(row['start']) >= int(row['ORF_start'])) &  (int(row['end']) >= int(row['ORF_end'])) & (int(row['start']) < int(row['ORF_end'])):
		#print("right")
                #print(row['start']," : ", row['end'])
                #print(row['ORF_start']," : ", row['ORF_end'])
		overlapp = "yes-right"
	elif (int(row['start']) >= int(row['ORF_start'])) &  (int(row['end']) <= int(row['ORF_end'])):
		#if "JABAIF010000738.1:2152" in row['Names']:
	#		print("inside")
#			print(row['start']," : ", row['end'])#
#			print(row['ORF_start']," : ", row['ORF_end'])
		overlapp = "yes-inside"
	elif (int(row['start']) <= int(row['ORF_start'])) & (int(row['end']) >= int(row['ORF_end'])) :
		#print("outside")
                #print(row['start']," : ", row['end'])
                #print(row['ORF_start']," : ", row['ORF_end'])
		overlapp = "yes-outside"
	else:
		overlapp = "no"
	Table_with_ORFs.loc[Table_with_ORFs['ORF_query'].eq(row['ORF_query']) & Table_with_ORFs['ORF_target'].eq(row['target']) & Table_with_ORFs['ORF_end'].eq(row['ORF_end']) & Table_with_ORFs['end'].eq(row['end']),"Overlapp_ORF_EVEs"]=overlapp
	if overlapp =="no":
		continue
	else:
		Table_with_ORFs.loc[Table_with_ORFs['ORF_query'].eq(row['ORF_query']) & Table_with_ORFs['ORF_target'].eq(row['ORF_target']) & Table_with_ORFs['ORF_end'].eq(row['ORF_end']) & Table_with_ORFs['end'].eq(row['end']),"ORF_perc"]= (int(row['ORF_end'])- int(row['ORF_start']))/ (int(row['end'])- int(row['start']))
	#print("\n")


# remove non-overlapping ORFS with EVEs 
import numpy as np
Table_with_ORFs.loc[Table_with_ORFs['Overlapp_ORF_EVEs'].eq("no"),"ORF_start"] = np.nan
Table_with_ORFs.loc[Table_with_ORFs['Overlapp_ORF_EVEs'].eq("no"),"ORF_end"] = np.nan

Table_with_ORFs.loc[Table_with_ORFs['Overlapp_ORF_EVEs'].eq("yes-outside") & Table_with_ORFs['ORF_perc'].lt(0.5),"ORF_start"] = np.nan
Table_with_ORFs.loc[Table_with_ORFs['Overlapp_ORF_EVEs'].eq("yes-outside") & Table_with_ORFs['ORF_perc'].lt(0.5),"ORF_end"] = np.nan

# find overlapping ORFs within the same scaffold

is_overlapped = lambda x: x['ORF_start'] >= x['ORF_end'].shift(fill_value=-1)
Table_with_ORFs['Overlapp_group'] = Table_with_ORFs.sort_values(['Scaffold_name', 'ORF_start', 'ORF_end']) \
                .groupby(['Species_name','Scaffold_name'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()

Table_with_ORFs.loc[Table_with_ORFs['ORF_start'].isna(),"Overlapp_group"] = np.nan

Table_with_ORFs=Table_with_ORFs.merge(Table_with_ORFs.groupby(["Overlapp_group"]).size().reset_index(name="N_Overlapp"),on="Overlapp_group",how="left")
Table_with_ORFs=Table_with_ORFs.merge(Table_with_ORFs.groupby(["Species_name","Scaffold_name","Cluster_hmmer"]).size().reset_index(name="N_dup"),on=["Species_name","Scaffold_name","Cluster_hmmer"],how="left")


Table_with_ORFs.loc[Table_with_ORFs['Prot_name'].eq("LbFVorf13"),"Prot_name"] = "LbFVorf13 (JmJC)"
Table_with_ORFs.loc[Table_with_ORFs['Prot_name'].eq("LbFVorf11"),"Prot_name"] = "LbFVorf11 (JmJC)"
#Table_with_ORFs.loc[Table_with_ORFs['Prot_name'].eq("LbFVorf44"),"Prot_name"] ="LbFVorf44 ()"

Cynipid_pseudo_table=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Clustering/HSP_cynipids_table.txt",sep="\t")
Table_with_ORFs['Pseudogenized']="no"
Table_with_ORFs.loc[Table_with_ORFs['ORF_start'].ge(1) & Table_with_ORFs['ORF_perc'].lt(0.5),"Pseudogenized"] = "yes"
Table_with_ORFs.loc[Table_with_ORFs['Names'].isin(Cynipid_pseudo_table.loc[Cynipid_pseudo_table['New_start'].eq("pseudo")]['Names']),"Pseudogenized"] = "yes"


#Ad Platygaster orseoliae EVEs informations 

#Porseoliae_tab=pd.read_csv("/beegfs/data/bguinet/these/Genomes/Platygaster_orseoliae/run_mmseqs2_with_PoFV/Porseoliae_to_PoFV.tab",sep=";")

#Porseoliae_tab=Porseoliae_tab[['scaf_name','count_repeat','count_eucaryote','cov_depth_candidat','cov_depth_BUSCO','pvalue_cov','pseudogenized']]

#Porseoliae_Table_with_ORFs=Table_with_ORFs.loc[Table_with_ORFs['Names2'].str.contains("PoEFV")]
#Porseoliae_Table_with_ORFs['Scaffold_name2']= Porseoliae_Table_with_ORFs['Scaffold_name'].str.replace("-.*","")
#Porseoliae_Table_with_ORFs['Scaffold_name']=Porseoliae_Table_with_ORFs['Scaffold_name2'].str.split('_').str[:-1].str.join('_')

#Porseoliae_Table_with_ORFs.drop(['count_repeat','count_eucaryote','cov_depth_candidat','cov_depth_BUSCO','pvalue_cov','pseudogenized'], axis=1, inplace=True)

#Porseoliae_Table_with_ORFs.merge(Porseoliae_tab,left_on="Scaffold_name",right_on="scaf_name",how="left")
#Porseoliae_Table_with_ORFs

#Remove the wrong HSP


Table_with_ORFs.to_csv(out4,sep=";",index=False)
print("File saved to : ",out4)


