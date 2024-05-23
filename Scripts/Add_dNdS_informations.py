
from Bio import SeqIO
import numpy as np
import pandas as pd 
import sys 
import argparse
import os,re
from multipy.fdr import tst


# Print out a message when the program is initiated.
print('------------------------------------------------------------------------------\n')
print('                        Add dN/dS informations                       .\n')
print('------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Gather HMMEr files')
parser.add_argument("-c", "--Cluster_file", help="The path of the Cluster file to add the dN/dS informations")
parser.add_argument("-o", "--Output_file", help="The desired output file name")
parser.add_argument("-fel", "--FEL_file", help="The FEL analyse table")

args = parser.parse_args()

#Example usage : 

"""
#
python3 Add_dNdS_informations.py -c /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event.tab \
  -o /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event_dNdS.tab \
  -fel /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/dNdS_analysis/FEL_analyse.tab

"""
# Variable that stores fasta sequences

Cluster_file=args.Cluster_file
Output_file=args.Output_file
Positive_FEL_tab=args.FEL_file
Cluster_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event.tab"
Output_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score_Event_dNdS.tab"
Positive_FEL_tab="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/dNdS_analysis/FEL_analyse.tab"

#Paths 
EVE_dNdS_path="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/dNdS_analysis/"
BUSCO_dNdS_path="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/BUSCO_dNdS/Paml_results/"

attribs = ['Clustername','Event','Mean_dNdS','Pvalue_dNdS','SE_dNdS']

print("Adding dN/dS informations for EVEs ...")
# Add all dN/dS out tab results 
EVE_dNdS_tab = pd.DataFrame(columns=attribs)
#open the hmmer result tables
for filename in os.listdir(EVE_dNdS_path):
    if filename.endswith(".out"):
      if 'dNds_Cluster' in filename:
        subtab=pd.read_csv(EVE_dNdS_path+filename,sep="\t")
        EVE_dNdS_tab=EVE_dNdS_tab.append(subtab)

EVE_dNdS_tab.columns=['Cluster_hmmer','Event','Mean_dNdS','Pvalue_dNdS','SE_dNdS']

Blast_table=pd.read_csv(Cluster_file,sep=";")
Blast_table=Blast_table.merge(EVE_dNdS_tab, on=['Cluster_hmmer','Event'],how='left')


#Blast_table.loc[Blast_table['Event'].ge(1)][['Names','Cluster_hmmer','Event','New_prot_name','Mean_dNdS','Pvalue_dNdS']].head(60)

# Add positive FEL analyse

Fel_table=pd.read_csv(Positive_FEL_tab,sep="\t")
Fel_table=Fel_table.loc[~Fel_table['codon'].isna()]

Blast_table['N_codon_Purifying']=np.nan
Blast_table['N_codon_Diversifying']=np.nan

Fel_table['p-value']=Fel_table['p-value'].astype(float)
Fel_table_grouped = Fel_table.groupby('Clustername_protname')


# FDR analysis 
list_nb=[0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01]

for i in list_nb:
         Fel_table['FDR_pvalue']=tst(Fel_table['p-value'], q=i)
         print(i,len(Fel_table.loc[Fel_table['FDR_pvalue'].astype(str)=="False"]))

# O.04 chosen 
Fel_table['FDR_pvalue']=tst(Fel_table['p-value'], q=0.04)

Fel_table_grouped = Fel_table.groupby('Clustername_protname')
# iterate over each group
for group_name, df_group in Fel_table_grouped:
    n_diversi=len(df_group.loc[(df_group['beta'] > df_group['alpha']) & (df_group['p-value'].le(0.05))])
    n_purif=len(df_group.loc[(df_group['beta'] < df_group['alpha']) & (df_group['p-value'].le(0.05))])
    name=re.sub("_LbFV.*","",df_group['Clustername_protname'].iloc[0])
    #Blast_table.loc[Blast_table['Cluster_hmmer'].eq(name)]
    Blast_table.loc[Blast_table['Cluster_hmmer'].eq(name),'N_codon_Purifying']=n_purif
    Blast_table.loc[Blast_table['Cluster_hmmer'].eq(name),'N_codon_Diversifying']=n_diversi


# Add Contrast FEL result 

Contrast_FEL_df = pd.DataFrame(columns=['Cluster_hmmer','Site','Partition','alpha','beta (background)','beta (Foreground)','alpha/beta foreground','alpha/beta background','subs (Foreground)','P-value (overall)','Q-value (overall)','Permutation p-value','Total branch length','Nb_pos_contrast_FEL_codon','Nb_neg_contrast_FEL_codon','alpha/beta ratio'])

for filename in os.listdir("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/dNdS_analysis"):
	if filename.endswith("FEL.csv"):
		subContrast_FEL_df=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/dNdS_analysis/"+filename,sep=",")
		filename=re.sub("_dNdS_.*","",filename)
		subContrast_FEL_df['Cluster_hmmer']=filename
		subContrast_FEL_df['alpha/beta foreground']= subContrast_FEL_df['beta (Foreground)']/subContrast_FEL_df['alpha']
		subContrast_FEL_df['alpha/beta background']= subContrast_FEL_df['beta (background)']/subContrast_FEL_df['alpha']
		Nb_codon_dNdS_sup1=len(subContrast_FEL_df.loc[subContrast_FEL_df['alpha/beta foreground'].gt(1.2) & subContrast_FEL_df['Q-value (overall)'].lt(0.2)])
		Nb_codon_dNdS_inf1=len(subContrast_FEL_df.loc[subContrast_FEL_df['alpha/beta foreground'].lt(0.8) & subContrast_FEL_df['Q-value (overall)'].lt(0.2)])
		subContrast_FEL_df['Nb_pos_contrast_FEL_codon']=Nb_codon_dNdS_sup1
		subContrast_FEL_df['Nb_neg_contrast_FEL_codon']=Nb_codon_dNdS_inf1
		subContrast_FEL_df['alpha/beta ratio']=subContrast_FEL_df['alpha/beta foreground']/subContrast_FEL_df['alpha/beta background']
		Contrast_FEL_df=pd.concat([subContrast_FEL_df,Contrast_FEL_df])

# Mean dN/dS of codon with negative dN/dS significantly higher than free-living viruses
len(Contrast_FEL_df.loc[Contrast_FEL_df['alpha/beta foreground'].lt(0.8) & Contrast_FEL_df['Q-value (overall)'].lt(0.2) & Contrast_FEL_df['alpha/beta ratio'].gt(0.1)])
Contrast_FEL_df.loc[Contrast_FEL_df['alpha/beta foreground'].lt(0.8) & Contrast_FEL_df['Q-value (overall)'].lt(0.2) & Contrast_FEL_df['alpha/beta ratio'].gt(0.1)]['alpha/beta foreground'].mean()
Contrast_FEL_df.loc[Contrast_FEL_df['alpha/beta foreground'].lt(0.8) & Contrast_FEL_df['Q-value (overall)'].lt(0.2) & Contrast_FEL_df['alpha/beta ratio'].gt(0.1)]['alpha/beta background'].mean()

# Mean dN/dS of codon with positive dN/dS significantly higher than free-living viruses
len(Contrast_FEL_df.loc[Contrast_FEL_df['alpha/beta foreground'].gt(1.2) & Contrast_FEL_df['Q-value (overall)'].lt(0.2) & Contrast_FEL_df['alpha/beta ratio'].gt(0.1)])
Contrast_FEL_df.loc[Contrast_FEL_df['alpha/beta foreground'].gt(1.2) & Contrast_FEL_df['Q-value (overall)'].lt(0.2) & Contrast_FEL_df['alpha/beta ratio'].gt(0.1)]['alpha/beta foreground'].mean()
Contrast_FEL_df.loc[Contrast_FEL_df['alpha/beta foreground'].gt(1.2) & Contrast_FEL_df['Q-value (overall)'].lt(0.2) & Contrast_FEL_df['alpha/beta ratio'].gt(0.1)]['alpha/beta background'].mean()

Contrast_FEL_df.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/dNdS_analysis/Contrast_FEL_table.txt",sep=";",index=False)

#Check stop codon within sequences
def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

AA_records_origin=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses.faa","fasta"))
AA_records_ORFs=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.aa","fasta"))

Blast_table['Pseudogenized']='NA'
for index, row in Blast_table.loc[~Blast_table['Names'].isna()].iterrows():
	if '(' in row['Names']:
		loci=row['Names']
		try:
			sequence=str(AA_records_origin[loci].seq)
		except:
			loci=re.sub("009025955.1","009025955",row['Names'])
			loci=re.sub("GCA_015476485.1","GCA_015476485",loci)
			sequence=str(AA_records_origin[loci].seq)
		if '*' in sequence:
			#print(">",loci,sep="")
			#print(sequence)
			pseudo="yes"
		else:
			pseudo="no"
		Blast_table.at[index,'Pseudogenized'] = pseudo


Blast_table.to_csv(Output_file,sep=";",index=False)
#
print("All dN/dS informations added to the file : ", Output_file)

# Create a reduced format for plotting 
Blast_table=Blast_table.sort_values(['Mean_dNdS', 'Cluster_hmmer'], ascending=[True, False])
Blast_table=Blast_table.loc[~Blast_table['Mean_dNdS'].isna()]
Blast_table = Blast_table.drop_duplicates(subset = ["Cluster_hmmer", "Event"])

Blast_table.loc[Blast_table['Cluster_hmmer'].str.contains("Cluster559-Cluster482",na=False),"Prot_name"]="lef-5"
sub_dNdS_tab=Blast_table[['Cluster_hmmer','Event','Tips4dNdS.array','Prot_name','Mean_dNdS','Pvalue_dNdS','SE_dNdS']]

#check
list_ORF_to_test=['LbFVorf72','LbFVorf83','LbFVorf108','LbFVorf87','LbFVorf94','LbFVorf5','LbFVorf60 (lcat)','LbFVorf107 (lef-4)','LbFVorf96 (lef-8)','LbFVorf78 (lef-9)','LbFVorf58 (DNApol)','LbFVorf68 (helicase)','LbFVorf92','LbFVorf85 (Ac81)','lef-5']


print("Checking dNdS for cynipoidea...")

for prot in  list_ORF_to_test:
	print(prot)
	#print(Blast_table.loc[Blast_table['Prot_name'].eq(prot) & Blast_table['Names2'].isin(['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV'])][['Cluster_hmmer','Nsp','Names2','Mean_dNdS','Pvalue_dNdS','SE_dNdS','Prot_name','Boot']]['Boot'].unique())
	Blast_table.loc[Blast_table['Prot_name'].eq('LbFVorf96 (lef-8)') & Blast_table['Names2'].isin(['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV'])][['Cluster_hmmer','Nsp','Names2','Mean_dNdS','Pvalue_dNdS','SE_dNdS','Prot_name','Boot']]
	print("\n")

sub_dNdS_tab.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/dNdS_analysis/All_EVE_dNdS_estimations.tab",sep=";",index=False)


# Load all BUSCO dN/dS estimations 
print("Adding BUSCO informations for EVEs ...")
# Add all dN/dS out tab results 
BUSCO_dNdS_tab = pd.DataFrame(columns=attribs)
#open the hmmer result tables
for filename in os.listdir(BUSCO_dNdS_path):
    if filename.endswith("_1.out"):
      if 'BUSCO' in filename:
        subtab=pd.read_csv(BUSCO_dNdS_path+filename,sep="\t")
        BUSCO_dNdS_tab=BUSCO_dNdS_tab.append(subtab)

BUSCO_dNdS_tab.columns=['Busco_name','Event','Mean_dNdS','Pvalue_dNdS','SE_dNdS']


BUSCO_dNdS_tab.to_csv(BUSCO_dNdS_path+"All_BUSCO_dNdS_estimations.tab",sep=";",index=False)

print("All BUSCO dN/dS estimations have been written to : ", BUSCO_dNdS_path+"All_BUSCO_dNdS_estimations.tab")

