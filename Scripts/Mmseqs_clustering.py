
import os
import argparse
import networkx as nx
import pandas as pd 
import numpy as np
import subprocess
from Bio import SeqIO 

# Print out a message when the program is initiated. traitrise
print('----------------------------------------------------------------------------------------------\n')
print('                        Extract filtred loci and run Clystering analysis                      .\n')
print('----------------------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Create Binary and cluster file ')
parser.add_argument("-b", "--blast_file", help="The Nr file")
parser.add_argument("-o", "--output", help="The updated blast output file")
args = parser.parse_args()


#Usage example python3 /beegfs/data/bguinet/Cynipoidea_project/Scripts/Mmseqs_clustering.py -b /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/NR_analyse/Blast_file_HSPs_NR.tab -o /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-cluster.tab

output_cluster_file=args.output

output_cluster_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-cluster.tab"


blast_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/NR_analyse/Blast_file_HSPs_NR.tab"
blast_tab=pd.read_csv(blast_file,sep=";")

All_seq_aa="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses.faa"
All_seq_db="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses_db"
All_seq_result="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses_result"
All_seq_tpm="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses_tpm"
All_seq_tab="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses_result.m8"

# Save all virus protein and candidate loci within a fasta file to run ALL vs ALL blastp

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

record_loci = to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/All_seqs/Fasta_viral_loci_ALL.faa", "fasta"))
with open(All_seq_aa,"w") as output :
  for virus in SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_Predicted_and_known_ORFs_new.fa", "fasta"):
    print('>',virus.id,sep="",file=output)
    print(virus.seq,file=output)
  for loci in blast_tab['full_name'].unique():
    print('>',record_loci[loci].id,sep="",file=output)
    print(record_loci[loci].seq,file=output)


# Run mmseqs all vs all analysis 

subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb " + All_seq_aa + " " + All_seq_db, shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs search " + All_seq_db + " " + All_seq_db + " " + All_seq_result + " " + All_seq_tpm +" --threads 10 -e 0.01 -s 7.5 ", shell=True)
subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs convertalis --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tlen,qcov,tcov,qlen' " + All_seq_db + " " + All_seq_db + " " + All_seq_result + " " + All_seq_tab, shell=True)

# Open the file 

df=pd.read_csv(All_seq_tab,sep="\t")
df.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qcov','tcov','qlen']

# Filter 
df=df.loc[df['qcov'].ge(0.25) & df['tcov'].ge(0.25) ]
df=df.loc[df['bits'].gt(50)]

# Create the graph from the dataframe
g = nx.Graph()

g.add_edges_from(df[['query','target']].itertuples(index=False))
new = list(nx.connected_components(g))

mapped =  {node: f'Cluster{cid + 1}' for cid, component in enumerate(new) for node in component}

cluster = pd.DataFrame({'Cluster': mapped.values(), 'Names':mapped.keys()})


#Add alignment informations 
cluster_with_ali_info=cluster.copy()
df2=df[['target','query','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']]
df2.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']
df3=df.append(df2)
df3=df3.sort_values(by=['evalue','bits'], ascending=[True,False])
df3=df3.loc[~(df3['query'] == df3['target'])]
df3 = df3.drop_duplicates(subset=['query'], keep='first')
cluster_with_ali_info=cluster_with_ali_info.merge(df3,right_on="query",left_on="Names")
#Get number of loci within clusters 
cluster_with_ali_info['Species_name'] =  cluster_with_ali_info['Names'].str.replace(".*_","")
cluster_with_ali_info['Nb_species'] = cluster_with_ali_info.groupby(['Cluster'])['Species_name'].transform('nunique')
cluster_with_ali_info_filtred = cluster_with_ali_info.loc[cluster_with_ali_info['Nb_species'].ge(4)]
print("Max evalue :", cluster_with_ali_info_filtred['evalue'].max())
print("Min pident :", cluster_with_ali_info_filtred['pident'].min())
print("Min bits :", cluster_with_ali_info_filtred['bits'].min())

#cluster.to_csv(output_cluster_file,sep=";",index=False)

#Create cluster files only for cluster with at least 4 species loci 
import re 
from Bio import SeqIO
#add length and Shorten the names :

cluster['Names2'] = np.where(cluster.Names.str.contains("_orf") , cluster['Names'].str.replace("_.*",""),cluster.Names)
cluster['Names2'] = np.where(cluster.Names.str.contains(":") , cluster.Names2,cluster['Names2'].str.replace(".*_",""))
cluster['Names2'] = np.where(cluster.Names.str.contains(":") ,cluster['Names2'].str.replace(".*:",""), cluster.Names2)

#Remove self match sequences 
#df=df[~(df['query']==df['target'])]

df_length1= df[['query','qlen','evalue']]
df_length1=df_length1.sort_values(by='evalue', ascending=True)
df_length1 = df_length1.drop_duplicates(subset = "query")
df_length1.columns=['Names','Length','evalue']

df_length2= df[['query','tlen','evalue']]
df_length2=df_length2.sort_values(by='evalue', ascending=True)
df_length2 = df_length2.drop_duplicates(subset = "query")
df_length2.columns=['Names','Length','evalue']

df_length=df_length1.append(df_length2)
df_length=df_length.sort_values(by='evalue', ascending=True)
df_length = df_length.drop_duplicates(subset = "Names")

cluster=cluster.merge(df_length,on="Names",how="left")

#Remove cluster without free-living viruses

Freeliving_viruses=['AcMNPV', 'LhFV', 'CcFV2', 'LbFV', 'DFV', 'CcFV1', 'PoFV', 'PcFV', 'esp', 'OrNV', 'tom', 'MdSGHV', 'DiNV', 'HgNV', 'GbNV', 'CuniNPV', 'EfFV', 'GpSGHV', 'DhNV', 'WSSV', 'HzNV-2', 'LdMNPV', 'PmNV', 'mau', 'ToNV', 'AmFV', 'kal', 'HzNV-1', 'CpV', 'NeseNPV']

out=cluster.groupby('Cluster')['Names2'].apply(lambda x: len(set(x) & set(Freeliving_viruses))).reset_index()
out=out.sort_values(by='Names2', ascending=False)
out['Cluster_contains_FL'] = np.where(out.Names2.ge(1), 'yes', 'no')
cluster = cluster.merge(out[['Cluster','Cluster_contains_FL']])

cluster=cluster.loc[cluster['Cluster_contains_FL'].eq("yes")]

#####

#Count number of duplicates species loci within each clusters in order to find possibly wrongly associated clusters 
# set up group
g = cluster.groupby('Cluster')
# get unique values
cluster['Number_unique_Names2'] = g['Names2'].transform('nunique')
# get non-duplicates
non_dup = g['Names2'].transform(lambda x: (~x.duplicated(False)).sum())
# duplicates = unique - non-duplicates
cluster['Nb_duplicated'] = cluster['Number_unique_Names2'] - non_dup

cluster['Ratio_duplicated_unique'] = cluster['Nb_duplicated']/cluster['Number_unique_Names2']
duplicated_clusters=cluster.loc[cluster['Ratio_duplicated_unique'].gt(0.30) & cluster['Nb_duplicated'].ge(4)]

# step1
# Redefine more stringent clusters for those sequences 

df_bis=pd.read_csv(All_seq_tab,sep="\t",header=None)
df_bis.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qcov','tcov','qlen']
df_bis=df_bis.loc[df_bis['qcov'].ge(0.25) | df_bis['tcov'].ge(0.25) ]
df_bis=df_bis.loc[df_bis['bits'].gt(47)]


# Create the graph from the dataframe
g_bis = nx.Graph()
g_bis.add_edges_from(df_bis[['query','target']].itertuples(index=False))
new_bis = list(nx.connected_components(g_bis))
mapped_bis =  {node: f'Cluster{cid + 1}' for cid, component in enumerate(new_bis) for node in component}
cluster_bis = pd.DataFrame({'Cluster': mapped_bis.values(), 'Names':mapped_bis.keys()})

#Add alignment informations 
cluster_with_ali_info_bis=cluster_bis.copy()
df2_bis=df_bis[['target','query','pident','alnlen','mismatch','gapopen','tstart','tend','qstart','qend','evalue','bits','qlen','tlen','qcov','tcov']]
df2_bis.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']
df3_bis=df_bis.append(df2_bis)
df3_bis=df3_bis.sort_values(by=['evalue','bits'], ascending=[True,False])
df3_bis=df3_bis.loc[~(df3_bis['query'] == df3_bis['target'])]
df3_bis = df3_bis.drop_duplicates(subset=['query'], keep='first')
cluster_with_ali_info_bis=cluster_with_ali_info_bis.merge(df3_bis,right_on="query",left_on="Names")
#Get number of loci within clusters 
cluster_with_ali_info_bis['Species_name'] =  cluster_with_ali_info_bis['Names'].str.replace(".*_","")
cluster_with_ali_info_bis['Nb_species'] = cluster_with_ali_info_bis.groupby(['Cluster'])['Species_name'].transform('nunique')
cluster_with_ali_info_filtred_bis = cluster_with_ali_info_bis.loc[cluster_with_ali_info_bis['Nb_species'].ge(4)]
print("Max evalue :", cluster_with_ali_info_filtred_bis['evalue'].max())
print("Min pident :", cluster_with_ali_info_filtred_bis['pident'].min())
print("Min bits :", cluster_with_ali_info_filtred_bis['bits'].min())


#Keep only cluster with too many duplicated sequences 
cluster_bis=cluster_bis.loc[cluster_bis['Names'].isin(duplicated_clusters['Names'].unique())]

#Remove from the original clusters the redifined clusters 
cluster=cluster.loc[~cluster['Names'].isin(cluster_bis['Names'].unique())]
cluster_bis['Cluster']=cluster_bis['Cluster']+'_redefined'
cluster=cluster.append(cluster_bis)

#cluster['Names2'] = np.where(cluster.Names.str.contains(":") , cluster['Names'].str.replace(".*:",""), cluster['Names'].str.replace(".*_",""))
cluster['Names2'] = np.where(cluster.Names.str.contains("_orf") , cluster['Names'].str.replace("_.*",""),cluster.Names)
cluster['Names2'] = np.where(cluster.Names.str.contains(":") , cluster.Names2,cluster['Names2'].str.replace(".*_",""))
cluster['Names2'] = np.where(cluster.Names.str.contains(":") ,cluster['Names2'].str.replace(".*:",""), cluster.Names2)


List_viral_sp= list(cluster['Names2'].unique())


# step 2 
# find still duplicated clusters 

#Count number of duplicates species loci within each clusters in order to find possibly wrongly associated clusters 
# set up group
g = cluster.groupby('Cluster')
# get unique values
cluster['Number_unique_Names2'] = g['Names2'].transform('nunique')
# get non-duplicates
non_dup = g['Names2'].transform(lambda x: (~x.duplicated(False)).sum())
# duplicates = unique - non-duplicates
cluster['Nb_duplicated'] = cluster['Number_unique_Names2'] - non_dup

cluster['Ratio_duplicated_unique'] = cluster['Nb_duplicated']/cluster['Number_unique_Names2']
duplicated_clusters=cluster.loc[cluster['Ratio_duplicated_unique'].gt(0.50) & cluster['Nb_duplicated'].ge(4)]

# Redefine more stringent clusters for those sequences 

df_bis=pd.read_csv(All_seq_tab,sep="\t",header=None)
df_bis.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qcov','tcov','qlen']
df_bis=df_bis.loc[df_bis['qcov'].ge(0.50) & df_bis['tcov'].ge(0.50) ]
df_bis=df_bis.loc[df_bis['bits'].gt(60)]


# Create the graph from the dataframe
g_bis = nx.Graph()
g_bis.add_edges_from(df_bis[['query','target']].itertuples(index=False))
new_bis = list(nx.connected_components(g_bis))
mapped_bis =  {node: f'Cluster{cid + 1}' for cid, component in enumerate(new_bis) for node in component}
cluster_bis = pd.DataFrame({'Cluster': mapped_bis.values(), 'Names':mapped_bis.keys()})

#Add alignment informations 
cluster_with_ali_info_bis=cluster_bis.copy()
df2_bis=df_bis[['target','query','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']]
df2_bis.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']
df3_bis=df_bis.append(df2_bis)
df3_bis=df3_bis.sort_values(by=['evalue','bits'], ascending=[True,False])
df3_bis=df3_bis.loc[~(df3_bis['query'] == df3_bis['target'])]
df3_bis = df3_bis.drop_duplicates(subset=['query'], keep='first')
cluster_with_ali_info_bis=cluster_with_ali_info_bis.merge(df3_bis,right_on="query",left_on="Names")




#Get number of loci within clusters 
cluster_with_ali_info_bis['Species_name'] =  cluster_with_ali_info_bis['Names'].str.replace(".*_","")
cluster_with_ali_info_bis['Nb_species'] = cluster_with_ali_info_bis.groupby(['Cluster'])['Species_name'].transform('nunique')
cluster_with_ali_info_filtred_bis = cluster_with_ali_info_bis.loc[cluster_with_ali_info_bis['Nb_species'].ge(4)]
print("Max evalue :", cluster_with_ali_info_filtred_bis['evalue'].max())
print("Min pident :", cluster_with_ali_info_filtred_bis['pident'].min())
print("Min bits :", cluster_with_ali_info_filtred_bis['bits'].min())

# Keep only cluster with too many duplicated sequences 
cluster_bis=cluster_bis.loc[cluster_bis['Names'].isin(duplicated_clusters['Names'].unique())]

#Remove from the original clusters the redifined clusters 
cluster=cluster.loc[~cluster['Names'].isin(cluster_bis['Names'].unique())]
cluster_bis['Cluster']=cluster_bis['Cluster']+'_redefined'
cluster=cluster.append(cluster_bis)

#cluster['Names2'] = np.where(cluster.Names.str.contains(":") , cluster['Names'].str.replace(".*:",""), cluster['Names'].str.replace(".*_",""))

cluster['Names2'] = np.where(cluster.Names.str.contains("_orf") , cluster['Names'].str.replace("_.*",""),cluster.Names)
cluster['Names2'] = np.where(cluster.Names.str.contains(":") , cluster.Names2,cluster['Names2'].str.replace(".*_",""))
cluster['Names2'] = np.where(cluster.Names.str.contains(":") ,cluster['Names2'].str.replace(".*:",""), cluster.Names2)


List_viral_sp= list(cluster['Names2'].unique())

g = cluster.groupby('Cluster')
# get unique values
cluster['Number_unique_Names2'] = g['Names2'].transform('nunique')
# get non-duplicates
non_dup = g['Names2'].transform(lambda x: (~x.duplicated(False)).sum())
# duplicates = unique - non-duplicates
cluster['Nb_duplicated'] = cluster['Number_unique_Names2'] - non_dup

cluster['Ratio_duplicated_unique'] = cluster['Nb_duplicated']/cluster['Number_unique_Names2']

# Other round 

out=cluster.groupby('Cluster')['Names2'].apply(lambda x: len(set(x) & set(Freeliving_viruses))).reset_index()
out=out.sort_values(by='Names2', ascending=False)
out['Cluster_contains_FL'] = np.where(out.Names2.ge(1), 'yes', 'no')
cluster = pd.merge(left=cluster, right=out[['Cluster','Cluster_contains_FL']],on="Cluster",how="left")

#cluster=cluster.loc[cluster['Cluster_contains_FL'].eq("yes")]

#Number cluster with at least 27 sequences : 
#len(cluster.loc[cluster_bis['Core_genes'].str.contains("yes")]['Cluster'].unique())

# Count cluster with N shared viral sequences 
out=cluster.groupby('Cluster')['Names2'].apply(lambda x: len(set(x) & set(List_viral_sp))).reset_index()
out=out.sort_values(by='Names2', ascending=False)
out['Core_genes'] = np.where(out.Names2.ge(27), 'yes', 'no')

cluster = cluster.merge(out[['Cluster','Core_genes']])
#Does the core gene cluster capture all viruses ? 
print( len(cluster.loc[cluster['Core_genes'].eq("yes")]['Names2'].unique()), " / ", len(List_viral_sp))
# Species not present in the remeaning clusters 
list( set(List_viral_sp) - set(list(cluster.loc[cluster['Core_genes'].eq("yes")]['Names2'].unique())))



# Add gene name when available 
Gene_names=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Filamentous_gene_names.tab",sep="\t")
Gene_names['Names']=Gene_names['Names']+"_LbFV"
Gene_names_cluster=cluster.merge(Gene_names,on="Names",how="left")
Gene_names_cluster=Gene_names_cluster.loc[~Gene_names_cluster['Prot_name'].isna()]
cluster=cluster.merge(Gene_names_cluster[['Cluster','Prot_name']], on="Cluster",how="left" )


#Check clusters 

#list_ORF_to_test=['LbFVorf72','LbFVorf83','LbFVorf108','LbFVorf87','LbFVorf94','LbFVorf5','LbFVorf60 (lcat)','LbFVorf107 (lef-4)','LbFVorf96 (lef-8)','LbFVorf78 (lef-9)','LbFVorf58 (DNApol)','LbFVorf68 (helicase)','LbFVorf92','LbFVorf85 (Ac81)']

#for orf in list_ORF_to_test:
#  print(orf)
#  print(len(cluster.loc[cluster['Prot_name'].eq(orf) & cluster['Names2'].isin('LhEFV','LbEFV','L





import re 
from Bio import SeqIO
#add length and Shorten the names :

del cluster['Length']
del cluster['evalue']

cluster=cluster.merge(df_length,on="Names",how="left")

def to_dict_remove_dups(sequences):
      return {record.id: record for record in sequences}

#


record_dict = to_dict_remove_dups(SeqIO.parse(All_seq_aa,"fasta"))



# Correct pif-1 for CcFV1  according to Annie
#cluster.loc[cluster['Names'].eq("CcFV1_32661_33227_+_CcFV1"),"Cluster"]=cluster.loc[cluster['Names'].str.contains("CcFV2_45457_4612")]['Cluster'].iloc[0]

#COrect pif-3 according to annie 
cluster.loc[cluster['Names'].str.contains("CcFV1_orf028"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79975_LbFV")]['Cluster'].iloc[0]

#Mannually add WSSV  into right clusters

cluster.loc[cluster['Names'].str.contains("YP_009220649"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79978")]['Cluster'].iloc[0] #DNA_pol
cluster.loc[cluster['Names'].str.contains("YP_009220510"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79983")]['Cluster'].iloc[0] #PIF0
cluster.loc[cluster['Names'].str.contains("YP_009220545"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79952")]['Cluster'].iloc[0] #PIF1
cluster.loc[cluster['Names'].str.contains("YP_009220486"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79972")]['Cluster'].iloc[0] #PIF2
cluster.loc[cluster['Names'].str.contains("YP_009220581"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79975")]['Cluster'].iloc[0] #PIF3
cluster.loc[cluster['Names'].str.contains("YP_009220481"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79963")]['Cluster'].iloc[0] #PIF5
cluster.loc[cluster['Names'].str.contains("YP_009220590"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ80017")]['Cluster'].iloc[0] #P33

#Manually add AmFV into right clusters

cluster.loc[cluster['Names'].str.contains("AKY03143"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79978")]['Cluster'].iloc[0]
cluster.loc[cluster['Names'].str.contains("AKY03146"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79983")]['Cluster'].iloc[0]
cluster.loc[cluster['Names'].str.contains("AKY03129"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79952")]['Cluster'].iloc[0]
cluster.loc[cluster['Names'].str.contains("AKY03169"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79972")]['Cluster'].iloc[0]
cluster.loc[cluster['Names'].str.contains("AKY03157"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79975")]['Cluster'].iloc[0]
cluster.loc[cluster['Names'].str.contains("AKY03126"),'Cluster'] = cluster.loc[cluster['Names'].str.contains("AQQ79963")]['Cluster'].iloc[0]


#Replace some names 
cluster['Names2']=cluster['Names2'].str.replace("ApBV","ApENV")
cluster['Names2']=cluster['Names2'].str.replace("BtBV","BtENV")
cluster['Names2']=cluster['Names2'].str.replace("CcBV","CcENV")
cluster['Names2']=cluster['Names2'].str.replace("CiBV","CiENV") 
cluster['Names2']=cluster['Names2'].str.replace("DnBV","DnENV")
cluster['Names2']=cluster['Names2'].str.replace("Platygaster_orseoliae","PoEFV")
cluster['Names2']=cluster['Names2'].str.replace("Dolichomitus","DoEFV")
cluster['Names2']=cluster['Names2'].str.replace("Eurytoma_brunniventris","EbENV")
cluster['Names2']=cluster['Names2'].str.replace("FaBV","FaENV")
cluster['Names2']=cluster['Names2'].str.replace("Leptopilina_boulardi_GCA_019393585.1","LbEFV")
cluster['Names2']=cluster['Names2'].str.replace("Leptopilina_clavipes","LcEFV")
cluster['Names2']=cluster['Names2'].str.replace("Leptopilina_heterotoma_GCA_015476425.1","LhEFV") 
cluster['Names2']=cluster['Names2'].str.replace("MdBV","MdENV")
cluster['Names2']=cluster['Names2'].str.replace("Melanaphis_sacchari","MsENV")
cluster['Names2']=cluster['Names2'].str.replace("NlBV","NlENV")
cluster['Names2']=cluster['Names2'].str.replace("Phanerotoma","PhENV")
cluster['Names2']=cluster['Names2'].str.replace("Rhoptromeris","RhEFV")
cluster['Names2']=cluster['Names2'].str.replace("Thrichoplasta","ThEFV")
cluster['Names2']=cluster['Names2'].str.replace("Trybliographa","TrEFV")
cluster['Names2']=cluster['Names2'].str.replace("VcBV","VcENV")
cluster['Names2']=cluster['Names2'].str.replace("TnBV","TnENV")


#Check clusters 

list_ORF_to_test=['LbFVorf72','LbFVorf83','LbFVorf108','LbFVorf87','LbFVorf94','LbFVorf5','LbFVorf60 (lcat)','LbFVorf107 (lef-4)','LbFVorf96 (lef-8)','LbFVorf78 (lef-9)','LbFVorf58 (DNApol)','LbFVorf68 (helicase)','LbFVorf92','LbFVorf85 (Ac81)']

for orf in list_ORF_to_test:
  print(orf)
  print(len(cluster.loc[cluster['Prot_name'].eq(orf) & cluster['Names2'].isin(['LhEFV','LbEFV','LcEFV','ThEFV','TrEFV','RhEFV'])]['Names'].unique()))
  print("\n")




list_cynipoidea =['LbEFV','LcEFV','LhEFV','RhEFV','PoEFV','DoEFV','ThEFV','TrEFV']

count_cluster= cluster[['Cluster','Names2']]
count_cluster = count_cluster.drop_duplicates(subset=['Cluster', 'Names2'], keep='first')

count_cluster['Cynip_count'] = count_cluster['Names2'].isin(list_cynipoidea).groupby(count_cluster['Cluster']).transform('sum')
count_cluster = count_cluster.drop_duplicates(subset=['Cluster'], keep='first')

cluster=cluster.merge(count_cluster[['Cluster','Cynip_count']],on="Cluster",how="left")


cluster.to_csv(output_cluster_file,sep=";",index=False)
print("Cluster tab written to : ", output_cluster_file)

#Write the cluster files and keep the longest sequence when paralogs  
grouped = cluster.groupby('Cluster')
for group_name, group in grouped:
        subcluster=group.sort_values(by='Length',ascending=False)
        subcluster=subcluster.drop_duplicates(subset='Names2', keep="first")
        if len(subcluster)>=4:
                with open ("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_mmseqs/"+group_name+".aa","w") as output: 
                        for index, row in subcluster.iterrows():
                                if row['Names2'] == "CoBV" :
                                        continue
                                else:
                                        print('>',row['Names2'],sep="",file=output)
                                        print(re.sub("\\*","",str(record_dict[row['Names']].seq)),file=output)

#
print( " All cluster fasta files written to : /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_mmseqs/" )
