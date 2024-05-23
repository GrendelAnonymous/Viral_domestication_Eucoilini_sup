
import sys
import pandas as pd 
import argparse
import re 
import ast

print('-----------------------------------------------------------------------------------\n')
print('                        Add Monophyletic table.\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add taxonomy informationsa blast file')
parser.add_argument("-M","--Monophyletic_file", help="The Hsp file produced by the rule Merge_Hsp_analysis and the Merge_HSP_sequences_within_clusters.py code")
parser.add_argument("-b","--blast_file", help="The fasta file with all informations")
parser.add_argument("-o","--out", help="The oufile with HSp and scaffold score informations")
args = parser.parse_args()


Blast_tab=args.blast_file
Outfile=args.out
Monop_tab=args.Monophyletic_file

Blast_tab="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score.tab"
Monop_tab="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Monophyletic_assessment/Monophyletic_table.tab"
Outfile="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event.tab"
#Blast_tab="/beegfs/data/bguinet/these/Clustering3/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_HSP_Score.m8"

#Monop_tab="/beegfs/data/bguinet/these/Monophyletic_assessment3/Monophyletic_tab.tab"
#Ex usage : python3 /beegfs/home/bguinet/these_scripts_2/Add_Monophyletic_table.py -b /beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred.m8 -M /beegfs/data/bguinet/these/Monophyletic_assessment_filtred/Monophyletic_tab.tab -o  /beegfs/data/bguinet/these/Clustering4/Viralprot_vs_Viral_loci_result_all_match_and_cluster_taxid_ontology_2LCA_HMMER_Augustus_Repeat_Env_pseudo_HSP_Score_Filtred_Mono.m8

#Load Monophyletic_file
Monophyletic_table =pd.read_csv(Monop_tab,sep=";")

Monophyletic_table = Monophyletic_table.rename({'cluster': 'Clustername'}, axis='columns')

#Manually change for orf72 

#Adding of other columns
Monophyletic_table[["Nsp_MRCA", "Nsp"]] = Monophyletic_table[["Nsp_MRCA", "Nsp"]].apply(pd.to_numeric)

Monophyletic_table['Nsp_losses']=Monophyletic_table['Nsp_MRCA']-Monophyletic_table['Nsp']

Monophyletic_table['Clustername']=Monophyletic_table['Clustername'].str.replace("_AA.fna.treefile","")

Monophyletic_table['Tips4dNdS.array']=Monophyletic_table['Tips4dNdS.array'].str.replace("_\\+__PoEFV","_+_PoEFV")
Monophyletic_table['Tips4dNdS.array']=Monophyletic_table['Tips4dNdS.array'].str.replace("_-__PoEFV","_-_PoEFV")

#Monophyletic_table['Clustername']=Monophyletic_table['Clustername'].str.replace("_AA*.","")

#Remove duplicates 
Monophyletic_table=Monophyletic_table.drop_duplicates(['Clustername','Event'],keep= 'last')

#Load Blast file 
Blast_table=pd.read_csv(Blast_tab,sep=";")

Blast_table['Names_bis']=Blast_table['Names3'].str.replace(":|\(|\)", '_')
Blast_table['Names_bis']=Blast_table['Names_bis'].str.replace(" \\[ORF\\]", '')
#Add Monophyletic informations:

Monophyletic_table['Tips4dNdS.array']=Monophyletic_table['Tips4dNdS.array'].str.replace("'","")

Monophyletic_table['Tips4dNdS.array'] = Monophyletic_table['Tips4dNdS.array'].astype(str).str.strip('[|]').str.split(',')

Monophyletic_table['Names_bis'] = Monophyletic_table['Tips4dNdS.array'] 

Monophyletic_table=Monophyletic_table.explode('Names_bis')


del Monophyletic_table['Clustername']
Blast_mono =  pd.merge(Blast_table,Monophyletic_table, how="outer", on=["Names_bis"])


#fill event if not event 
"""
list_cluster=Monophyletic_table['Clustername'].unique()


list_cluster=  Blast_mono.loc[~Blast_mono['Clustername'].isin(list_cluster)]['Clustername'].unique()

def getgroup(List_Clustername):
    lst=[]
    for x in List_Clustername:
        m=Blast_mono['Clustername'].eq(x)
        if m.any():
            lst.append(Blast_mono[m].groupby(['Clustername','Species_name'],sort=False).ngroup()+1)
    return pd.concat(lst)

print(Blast_mono['Event'])
Blast_mono['Event']=pd.Series(Blast_mono.index.map(getgroup(list_cluster))).fillna(Blast_mono['Event']).astype(int)
"""
#Blast_mono = Blast_mono.loc[~(Blast_mono['Event'].isna() & Blast_mono['Scaffold_score'].isin(['A','B','C','D']) & Blast_mono['NSPtot'].gt(0))]

Blast_mono.to_csv(Outfile,index=False,sep=";")

print("Output file written to :", Outfile)
