# Create clusters 
import pandas as pd 
import re 
from Bio import SeqIO
import argparse
import os 
import subprocess


Output_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Syntenie.tab"
Main_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score.tab"
Main_table=pd.read_csv(Main_file,sep=";")

print("Blast table oppened")
# keep only scaffolds of interest

Main_table=Main_table.loc[Main_table['Scaffold_score'].isin(["A","B","C","D"])]
Main_table=Main_table.loc[Main_table['Names2'].isin(['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV'])]


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
  for loci in SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Genomes/LbFV/LbFV_free.fna","fasta"):
    print(">LbFV:",loci.id,sep="",file=output)
    print(loci.seq,file=output)
  for loci in SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Genomes/LhFV/LhFV_free.fna","fasta"):
    print(">LhFV:",loci.id,sep="",file=output)
    print(loci.seq,file=output)
  #Add also EVEs in dna format 
  for loci in SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses.fna","fasta"):
    for sp in ['Leptopilina_heterotoma_GCA_015476425.1','Leptopilina_boulardi_GCA_019393585.1', 'Trybliographa','Rhoptromeris', 'Leptopilina_clavipes', 'Thrichoplasta']:
      if sp in loci.id:
        #print(">",record_dict[loci].id,sep="",file=output)
        #print(record_dict[loci].seq,file=output)
        print(">",loci.id,sep="",file=output)
        print(loci.seq,file=output)
  
print("All scaffold extracted")

if os.path.exists("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result.m8"):
  print("Blast file already exist")
else:
  subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred.dna /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_db", shell=True)
  subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs search /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_db /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_db /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_tpm --threads 20 --search-type 4", shell=True)
  subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs convertalis --search-type 4  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tlen,qlen,qcov,tcov' /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_db /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_db /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result.m8" ,shell=True)

# Open blast synteny 

print("openning the huge synteny blast all vs all file...")
Synteny_table=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result.m8",sep="\t",header=None)
Synteny_table.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']
print("Huge synteny blast all vs all file oppened ...")

print("processing...1")
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
print("processing...2")
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

print("processing...3")

is_overlapped = lambda x: x['qstart'] >= x['qend'].shift(fill_value=-1)
Synteny_table['group'] = Synteny_table.sort_values(['query', 'qstart', 'qend']) \
    .groupby(['query','target'],as_index=False).apply(is_overlapped).droplevel(0).cumsum()

Synteny_table = Synteny_table.sort_values(['group', 'evalue', 'bits'], ascending=[True, True, False]) \
    .groupby(Synteny_table['group']).head(1)

print("processing...4")
is_overlapped = lambda x: x['tstart'] >= x['tend'].shift(fill_value=-1)
Synteny_table['group2'] = Synteny_table.sort_values(['target', 'tstart', 'tend']) \
    .groupby(['target','query'],as_index=False).apply(is_overlapped).droplevel(0).cumsum()

Synteny_table = Synteny_table.sort_values(['group2', 'evalue', 'bits'], ascending=[True, True, False]) \
    .groupby(Synteny_table['group2']).head(1)

Synteny_table.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Synteny_analysis/All_candidate_scaffold_filtred_all_vs_all_result_reduced.m8",sep=";",index=False)
#
