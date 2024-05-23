import pandas as pd 
import sys 
import argparse
import os
from collections import defaultdict
import pandas as pd
from Bio import SearchIO
from Bio import SeqIO 
import re 
import subprocess
from statistics import mean
from multipy.fdr import tst
import numpy as np

# Print out a message when the program is initiated.
print('------------------------------------------------------------------------------\n')
print('                        Add Metaeuk and TE results                          .\n')
print('------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Gather HMMEr files')
parser.add_argument("-b", "--blast_cluster_file", help="The blast cluster file")
parser.add_argument("-out1", "--out1", help="The first output containing clusters informations")
parser.add_argument("-out2", "--out2", help="The first output containing clusters  and coverages informations ")
parser.add_argument("-mk1", "--mk1", help="The pred metaeuk file ")
parser.add_argument("-mk2", "--mk2", help="The contig metaeuk file ")
parser.add_argument("-r", "--repeat", help="The repeat file ")

args = parser.parse_args()

#Example usage : 

"""
#
python3 Add_cov_metaeuk_TE_score.py  -b /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs.tab \
-out1 /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov.tab \
-out2 /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score.tab \
-mk1 /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/Taxa_Metaeuk_preds_results_tax_per_pred.tsv \
-mk2 /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/Taxa_Metaeuk_preds_results_tax_per_contig.tsv \
-r /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/All_candidate_scaffold_run_repeat_result.m8

"""
# Variable that stores fasta sequences

#path=args.hmmer_path
blast_file=args.blast_cluster_file
output_file1=args.out1
output_file2=args.out2
Repeat_file=args.repeat
MEuk_file_pred=args.mk1
MEuk_file_contig=args.mk2

blast_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs.tab"
output_file1="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov.tab"
output_file2="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs_Cov_Metaeuk_TE_Score.tab"
Repeat_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/All_candidate_scaffold_run_repeat_result.m8"
MEuk_file_pred="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/Taxa_Metaeuk_preds_results_tax_per_pred.tsv"
MEuk_file_contig="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/Taxa_Metaeuk_preds_results_tax_per_contig.tsv"

## Add coverage information 
Blast_table = pd.read_csv(blast_file,sep=";")

Blast_table['Names']=Blast_table['Names'].str.replace("GCA_009025955","GCA_009025955.1")
Blast_table['Species_name']=Blast_table['Species_name'].str.replace("GCA_009025955","GCA_009025955.1")
Blast_table['Names']=Blast_table['Names'].str.replace("GCA_015476485","GCA_015476485.1")
Blast_table['Species_name']=Blast_table['Species_name'].str.replace("GCA_015476485","GCA_015476485.1")

list_cynipoidea =['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV']
list_filamentous = ['LbFV','EfFV','PoFV','PcFV','CcFV1','CcFV2','LhFV']
list_cov_scaff_viral=[]

if os.path.exists("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Coverage_informations.tab") :
  candidates_cov_scaff=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Coverage_informations.tab",sep=";")
else:
  for sp in list_cynipoidea:
          print("\n")
          print("Processing ", sp)
          sub_Blast_table=Blast_table.loc[Blast_table['Names2'].str.contains(sp)]
          species=sub_Blast_table['Names'].str.replace(".*:","").iloc[0]
          print(species)
          #In order to open the file correctly
          subprocess.call(" sed -i 's@# Busco@Busco@g' /beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/assembly/BUSCO_results_hymenoptera/run_hymenoptera_odb10/full_table.tsv", shell=True)
          if sp == "LhEFV":
                  BUSCO_table=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/assembly/BUSCO_results_hymenoptera/run_hymenoptera_odb10/full_table.tsv",comment="#",sep="\t")
                  heterotoma_tab=pd.DataFrame(columns=['Scaffold_name','Info'])
                  for scf in SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/Leptopilina_heterotoma_GCA_015476425.1/assembly/Leptopilina_heterotoma_GCA_015476425.1_final_assembly.fna", "fasta"):
                          heterotoma_tab=heterotoma_tab.reset_index(drop=True)
                          heterotoma_tab = heterotoma_tab.append({'Scaffold_name':scf.id, 'Info':scf.description},ignore_index=True)
                  BUSCO_table_bis=BUSCO_table.loc[BUSCO_table['Sequence'].isin(heterotoma_tab.loc[heterotoma_tab['Info'].str.contains("whole genome shotgun sequence")]['Scaffold_name'])]
                  #BUSCO_table=BUSCO_table.loc[~BUSCO_table['Species_name'].str.contains("hetero")]
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
          if sp == 'LhEFV':
            list_scaff_keep=[]
            for scf in SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/Leptopilina_heterotoma_GCA_015476425.1/assembly/Leptopilina_heterotoma_GCA_015476425.1_final_assembly.fna", "fasta"):
              list_scaff_keep.append(scf.id)
            BUSCO_table=BUSCO_table.loc[BUSCO_table['Scaffold_name'].isin(list_scaff_keep)]
          Mean_busco_GC=BUSCO_table.loc[~BUSCO_table['GC'].isna()]['GC'].mean()
          Mean_busco_cov=BUSCO_table.loc[~BUSCO_table['Mean_cov'].isna()]['Mean_cov'].mean()
          print("Mean_busco_cov : ", Mean_busco_cov)
          Median_busco_cov=BUSCO_table.loc[~BUSCO_table['Median_cov'].isna()]['Median_cov'].mean()
          print("Median_busco_cov : ", Median_busco_cov)
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
                  list_cov_scaff_viral.append({'Species_name':species,'Scaffold_name':row['Scaffold_name'],"Mean_GC_candidat":row['GC'], "Mean_GC_BUSCO":Mean_busco_GC,'Mean_cov_depth_candidat':row['Mean_cov'],'Median_cov_depth_candidat':row['Median_cov'],'Mean_cov_depth_BUSCO':Mean_busco_cov,'Median_cov_depth_BUSCO':Median_busco_cov,'pvalue_cov_mean':pvalue_mean, 'pvalue_cov_median':pvalue_median,})
          print(pd.DataFrame(list_cov_scaff_viral))
  candidates_cov_scaff=pd.DataFrame(list_cov_scaff_viral)  #
  candidates_cov_scaff.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Coverage_informations.tab",sep=";",index=False)


# FDR analysis
"""
list_nb=[0.1,0.09,0.08,0.07,0.06,0.05,0.04,0.03,0.02,0.01]

for i in list_nb:
         candidates_cov_scaff['FDR_pvalue_cov_mean']=tst(candidates_cov_scaff['pvalue_cov_mean'], q=i)
         print(i,len(candidates_cov_scaff.loc[candidates_cov_scaff['FDR_pvalue_cov_mean'].astype(str)=="False"]))
"""
#Here we take 0.04
FDR_with_cov_inf= candidates_cov_scaff.drop_duplicates(subset=['Scaffold_name'], keep='first')
FDR_with_cov_inf['FDR_pvalue_cov_mean']=tst(FDR_with_cov_inf['pvalue_cov_mean'], q=0.04)
FDR_with_cov_inf['FDR_pvalue_cov_median']=tst(FDR_with_cov_inf['pvalue_cov_median'], q=0.04)
FDR_with_cov_inf=FDR_with_cov_inf[['Scaffold_name','Species_name','FDR_pvalue_cov_mean','FDR_pvalue_cov_median']]
#FDR_with_cov_inf.to_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Coverage_informations.tab",sep=";",index=False)

Blast_table['Scaffold_name']=Blast_table['Names'].str.replace(":.*","")


Blast_table=Blast_table.merge(FDR_with_cov_inf,on=['Scaffold_name','Species_name'],how="left")

Blast_table = Blast_table.merge(candidates_cov_scaff, on=["Scaffold_name","Species_name"],how="left")   


Blast_table['FDR_cov_global']=np.nan

for index, row in Blast_table.iterrows():
  if row['FDR_pvalue_cov_mean']==False:
    Blast_table.at[index,'FDR_cov_global'] = "False"
  elif row['FDR_pvalue_cov_median']==False:
    Blast_table.at[index,'FDR_cov_global'] = "False"
  else:
    Blast_table.at[index,'FDR_cov_global'] = "True"
      
Blast_table = Blast_table.drop_duplicates()

Blast_table.to_csv(output_file1,sep=";",index=False)
print("File saved to : ",output_file1)

# Save the file for plotting heatmap gene content 

list_ORF_to_test=['LbFVorf72','LbFVorf83','LbFVorf108','LbFVorf87','LbFVorf94','LbFVorf5','LbFVorf60 (lcat)','LbFVorf107 (lef-4)','LbFVorf96 (lef-8)','LbFVorf78 (lef-9)','LbFVorf58 (DNApol)','LbFVorf68 (helicase)','LbFVorf92','LbFVorf85 (Ac81)']

for sp in ['GCA_019393585.1','GCA_015476425.1','Thrichoplasta','Trybliographa','Rhoptromeris','Leptopilina_clavipes']:
  for orf in list_ORF_to_test:
    if len(Blast_table.loc[Blast_table['Prot_name'].eq(orf) & Blast_table['Names'].str.contains(sp) & Blast_table['FDR_cov_global'].eq("False")]['Names'].unique()) ==0:
      print(orf ," - ", sp)
  print("\n")



#############################################################
###    Look for ET/EUK gene along scaffolds with EVEs   #####
#############################################################

Blast_table=pd.read_csv(output_file1,sep=";")
###########################
## Transposable Elements ##
###########################

# Filter repeat results 

TE_table=pd.read_csv(Repeat_file,sep="\t",header=None)
TE_table.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','qlen','tlen','tcov']
TE_table['query']=TE_table['query'].str.replace("GCA_009025955","GCA_009025955.1")
TE_table['query']=TE_table['query'].str.replace("GCA_015476485","GCA_015476485.1")
TE_table=TE_table.loc[TE_table['bits'].gt(60)]
TE_table=TE_table.loc[TE_table['tcov'].gt(0.25)]


sLength=TE_table.shape[0]
TE_table['strand']=np.where(TE_table["qstart"]>TE_table["qend"],'-','+')
TE_table[['qstart', 'qend']] = np.sort(TE_table[['qstart', 'qend']].values, axis=1)
m = TE_table['strand'].eq('-')
TE_table = TE_table.sort_values(['query', 'qend', 'qstart'])
c1 = TE_table['query'].shift() != TE_table['query']
c2 = TE_table['qend'].shift() - TE_table['qstart'] < 0
TE_table['overlap'] = (c1 | c2).cumsum()
TE_table=TE_table.sort_values(['tlen'], ascending=False).groupby('overlap').first()
TE_table_count=TE_table.groupby(['query']).size().reset_index(name='count')
TE_table_count.columns=['Scaffold_and_Species_name','repeat_count']

# Filter metaeuk results 
Metaeuk_per_pred_table=pd.read_csv(MEuk_file_pred,sep="\t",header=None)
Metaeuk_per_pred_table.columns=['Protein','taxid','rank','prediction','taxlineages']
Metaeuk_per_pred_table['Protein']=Metaeuk_per_pred_table['Protein'].str.replace("GCA_009025955","GCA_009025955.1")
Metaeuk_per_pred_table['Protein']=Metaeuk_per_pred_table['Protein'].str.replace("GCA_015476485","GCA_015476485.1")
Metaeuk_per_pred_table=Metaeuk_per_pred_table.loc[Metaeuk_per_pred_table['taxlineages'].str.contains("Euk",na=False)]
Metaeuk_per_pred_table['Scaffold_and_Species_name']=Metaeuk_per_pred_table['Protein'].str.replace("\\|\\+.*","")
Metaeuk_per_pred_table['Scaffold_and_Species_name']=Metaeuk_per_pred_table['Scaffold_and_Species_name'].str.replace("\\|-.*","")
Metaeuk_per_pred_table['Scaffold_and_Species_name']=Metaeuk_per_pred_table['Scaffold_and_Species_name'].str.split('|', n=1).str.get(-1)
Metaeuk_per_pred_table_count=Metaeuk_per_pred_table.groupby(['Scaffold_and_Species_name']).size().reset_index(name='euk_count')

Metaeuk_per_contig_table=pd.read_csv(MEuk_file_contig,sep="\t",header=None)
Metaeuk_per_contig_table[0]=Metaeuk_per_contig_table[0].str.replace("GCA_009025955","GCA_009025955.1")
Metaeuk_per_contig_table[0]=Metaeuk_per_contig_table[0].str.replace("GCA_015476485","GCA_015476485.1")
Metaeuk_per_contig_table=Metaeuk_per_contig_table.loc[Metaeuk_per_contig_table[8].str.contains("Eukaryot",na=False)]
Metaeuk_per_contig_table['Metaeuk_vote']="yes"
Metaeuk_per_contig_table=Metaeuk_per_contig_table[[0,'Metaeuk_vote']]
Metaeuk_per_contig_table.columns=['Scaffold_and_Species_name','Metaeuk_vote']
Metaeuk_per_pred_table_count=Metaeuk_per_pred_table_count.merge(Metaeuk_per_contig_table,on="Scaffold_and_Species_name",how="outer")



#Merge both 
TE_and_Metaeuk_per_pred_table_count=Metaeuk_per_pred_table_count.merge(TE_table_count,on="Scaffold_and_Species_name",how="outer")

# Merge with the main table result 
Blast_table['Scaffold_and_Species_name']=Blast_table['Species_name']+":"+Blast_table['Scaffold_name']
Blast_table=Blast_table.merge(TE_and_Metaeuk_per_pred_table_count,on="Scaffold_and_Species_name",how='outer')
Blast_table.repeat_count.fillna(0,inplace=True)
Blast_table.Metaeuk_vote.fillna("no",inplace=True)
Blast_table.euk_count.fillna(0,inplace=True)




# Calculate the scaffolds scores 

Blast_table['euk_count_plus_repeat']=Blast_table['euk_count']+Blast_table['repeat_count']
Blast_table.loc[Blast_table['euk_count'].eq(0) & Blast_table['FDR_pvalue_cov_median'].astype(str).str.contains("True") & Blast_table['FDR_pvalue_cov_mean'].astype(str).str.contains("True"),'Scaffold_score' ]='X'
Blast_table.loc[Blast_table['repeat_count'].eq(0) & Blast_table['FDR_pvalue_cov_median'].astype(str).str.contains("True") & Blast_table['FDR_pvalue_cov_mean'].astype(str).str.contains("True"),'Scaffold_score']='X'

Blast_table.loc[Blast_table['euk_count'].gt(0) & Blast_table['euk_count'].lt(5) & Blast_table['FDR_pvalue_cov_median'].astype(str).str.contains("True") & Blast_table['FDR_pvalue_cov_mean'].astype(str).str.contains("True"),'Scaffold_score']='F'
Blast_table.loc[Blast_table['repeat_count'].gt(0) & Blast_table['repeat_count'].lt(2) & Blast_table['FDR_pvalue_cov_median'].astype(str).str.contains("True") & Blast_table['FDR_pvalue_cov_mean'].astype(str).str.contains("True"),'Scaffold_score']='F'

Blast_table.loc[Blast_table['euk_count'].ge(5) & Blast_table['FDR_pvalue_cov_median'].astype(str).str.contains("True") & Blast_table['FDR_pvalue_cov_mean'].astype(str).str.contains("True"),'Scaffold_score']='D'
Blast_table.loc[Blast_table['repeat_count'].gt(0) & Blast_table['FDR_pvalue_cov_median'].astype(str).str.contains("True") & Blast_table['FDR_pvalue_cov_mean'].astype(str).str.contains("True"),'Scaffold_score']='D'

Blast_table.loc[Blast_table['euk_count'].eq(0) & Blast_table['FDR_pvalue_cov_median'].astype(str).str.contains("False") | Blast_table['euk_count'].eq(0) & Blast_table['FDR_pvalue_cov_mean'].astype(str).str.contains("False"),'Scaffold_score']='C'
Blast_table.loc[Blast_table['repeat_count'].eq(0) & Blast_table['FDR_pvalue_cov_median'].astype(str).str.contains("False") | Blast_table['repeat_count'].eq(0) & Blast_table['FDR_pvalue_cov_mean'].astype(str).str.contains("False"),'Scaffold_score']='C'

Blast_table.loc[Blast_table['euk_count'].ge(1) & Blast_table['FDR_pvalue_cov_median'].astype(str).str.contains("False") | Blast_table['euk_count'].ge(1) & Blast_table['FDR_pvalue_cov_mean'].astype(str).str.contains("False") ,'Scaffold_score']='A'
Blast_table.loc[Blast_table['repeat_count'].ge(1) &Blast_table['FDR_pvalue_cov_median'].astype(str).str.contains("False") |Blast_table['repeat_count'].ge(1) & Blast_table['FDR_pvalue_cov_mean'].astype(str).str.contains("False") ,'Scaffold_score']='A'

# Save results 



#Change name to fit phylogeny ORF label name and Blast_table names
Blast_table['Names3']=Blast_table['Names']
for index, row in Blast_table.iterrows():
	if row['Best_hit_ORF_perc']> 0.70:
           Blast_table.at[index,'Names3']= row['Scaffold_name']+":"+str(int(row['ORF_start']))+"-"+str(int(row['ORF_end']))+"("+str(row['ORF_strand'])+"):"+row['Species_name']+' [ORF]'

Blast_table.to_csv(output_file2,sep=";",index=False)

"""
Blast_table=pd.read_csv(output_file2,sep=";")

# Create cov gc plot table for R scripyt
Blast_table=pd.read_csv(output_file2,sep=";")

list_cynipoidea =['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV']
Blast_table['Species_name']=Blast_table['Names'].str.replace(".*:","")
Blast_table['Scaffold_name']=Blast_table['Names'].str.replace(":.*","")
BUSCO_table_and_EVE=pd.DataFrame(columns=['Cluster_hmmer','Species_name', 'Scaffold_name', 'start','end','strand','Median_cov', 'Mean_cov', 'Scaffold_length', 'GC','Type','euk_count','repeat_count','Filamentous_count','Nb_ORFS','count_Eukaryota','count_Bacteria','count_Virus','Best_NR','Cluster_Perc_viral_vs_other','Cluster_count_Viral','Cluster_count_Uncertain','Cluster_count_Not_Viral','FDR_pvalue_cov_mean', 'FDR_pvalue_cov_median','Prot_name','Names'])
for species in  list(Blast_table.loc[Blast_table['Names2'].isin(list_cynipoidea)]['Species_name'].unique()):
        BUSCO_table=pd.DataFrame(columns=['Cluster_hmmer','Species_name', 'Scaffold_name', 'start','end','strand','Median_cov', 'Mean_cov', 'Scaffold_length', 'GC','euk_count','repeat_count'])
        subBUSCO_table=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/assembly/BUSCO_results_hymenoptera/run_hymenoptera_odb10/full_table_cov_gc.tsv",sep=";")
        subBUSCO_table['Species_name']=species
        subBUSCO_table['Type']="BUSCO"
        subBUSCO_table['Filamentous_count']=0
        subBUSCO_table['Nb_ORFS']=1
        subBUSCO_table['Cluster_Perc_viral_vs_other']=1000000
        subBUSCO_table['FDR_pvalue_cov_mean']="False"
        subBUSCO_table['FDR_pvalue_cov_median']="False"
        subBUSCO_table['Prot_name']="BUSCO"
        subBUSCO_table['Names']="BUSCO"
        subBUSCO_table['Cluster_hmmer']="BUSCO"
        subBUSCO_table['euk_count']=0
        subBUSCO_table['repeat_count']=0
        subBUSCO_table=subBUSCO_table[['Cluster_hmmer','Species_name', 'Scaffold_name','Gene Start', 'Gene End', 'Strand', 'Median_cov', 'Mean_cov', 'Scaffold_length', 'GC','Type','euk_count','repeat_count','Filamentous_count','Nb_ORFS','Cluster_Perc_viral_vs_other','FDR_pvalue_cov_mean', 'FDR_pvalue_cov_median','Names']]
        subBUSCO_table.columns=['Cluster_hmmer','Species_name', 'Scaffold_name', 'start','end','strand','Median_cov', 'Mean_cov', 'Scaffold_length', 'GC','Type','euk_count','repeat_count','Filamentous_count','Nb_ORFS','Cluster_Perc_viral_vs_other','FDR_pvalue_cov_mean', 'FDR_pvalue_cov_median','Names']
        #print(subBUSCO_table)
        BUSCO_table=BUSCO_table.append(subBUSCO_table)
        subBlast_table=Blast_table.loc[Blast_table['Species_name']==species]
        subBlast_table['Names3']=subBlast_table['Names'].str.replace("\\(.*","")
        subBlast_table['end']=subBlast_table['Names3'].str.replace(".*-","")
        subBlast_table['start']=subBlast_table['Names3'].str.replace(".*:","")
        subBlast_table['start']=subBlast_table['start'].str.replace("-.*","")
        subBlast_table['strand']=subBlast_table['Names'].str.replace("\\).*","")
        subBlast_table['strand']=subBlast_table['strand'].str.replace(".*\\(","")
        subBlast_table['Type']="EVE"
        subBlast_table['Nb_ORFS'] = subBlast_table.groupby(['Scaffold_name'])['Names'].transform('nunique')
        subBlast_table=subBlast_table[['Cluster_hmmer','Species_name','Scaffold_name','start','end','strand','Mean_cov_depth_candidat', 'Median_cov_depth_candidat', 'Scaffold_length', 'Mean_GC_candidat','Type','euk_count','repeat_count','Filamentous_count','Nb_ORFS','count_Eukaryota','count_Bacteria','count_Virus','Best_NR','Cluster_Perc_viral_vs_other','Cluster_count_Viral','Cluster_count_Uncertain','Cluster_count_Not_Viral','FDR_pvalue_cov_mean', 'FDR_pvalue_cov_median','Prot_name','Names']]
        subBlast_table.columns=['Cluster_hmmer','Species_name', 'Scaffold_name', 'start','end','strand','Mean_cov', 'Median_cov', 'Scaffold_length', 'GC','Type','euk_count','repeat_count','Filamentous_count','Nb_ORFS','count_Eukaryota','count_Bacteria','count_Virus','Best_NR','Cluster_Perc_viral_vs_other','Cluster_count_Viral','Cluster_count_Uncertain','Cluster_count_Not_Viral','FDR_pvalue_cov_mean', 'FDR_pvalue_cov_median','Prot_name','Names']
        #print(subBlast_table)
        subBUSCO_table_and_EVE=subBlast_table.append(BUSCO_table)
        BUSCO_table_and_EVE=BUSCO_table_and_EVE.append(subBUSCO_table_and_EVE)
        #print(subBUSCO_table_and_EVE)

"""
