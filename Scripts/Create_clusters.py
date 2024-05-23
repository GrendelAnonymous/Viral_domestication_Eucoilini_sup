

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
parser.add_argument("-o", "--output_file", help="The hmmer cluster output  file")

#Example usage : 

"""
#
python3 Create_clusters.py -b /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster.tab \
    -o /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs.tab
"""

# Variable that stores fasta sequences
args = parser.parse_args()
blast=args.blast_cluster_file
output=args.output_file

blast_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster.tab"
output_file="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Blast_file_HSPs_NR_Mmseqs-Hmmer_cluster_ORFs.tab"
Blast_table=pd.read_csv(blast_file,sep=";")

Blast_table['Species_name']=Blast_table['Names'].str.replace(".*:","")
Blast_table['Species_scaffold_name']=Blast_table['Species_name']+":"+Blast_table['Scaffold_name']

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

############################################
# Extract all scaffold containing candidates #
############################################

Loaded_scaffolds=[]
if os.path.exists("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.dna") :
        Loaded_scaffolds_dic= to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.dna","fasta"))
        Loaded_scaffolds=list(Loaded_scaffolds_dic.keys())

# Manually add scaffolds 
"""
tot=len(Blast_table['Species_name'].unique())
n=0
with open ("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.dna","a") as output:
        for assembly in Blast_table.loc[~ (Blast_table['Species_scaffold_name'].isin(Loaded_scaffolds)) & ~ (Blast_table['Scaffold_name'].isna()) ]['Species_name'].unique():
                print(" recovering scaffold of the species : ", assembly)
                subBlast_table=Blast_table.loc[Blast_table['Species_name'].eq(assembly)]
                subBlast_table=subBlast_table.loc[~subBlast_table['Species_name'].isin(Loaded_scaffolds)]
                try:
                        assembly_genome=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+assembly+"/assembly/"+assembly+"_final_assembly.fna", "fasta"))
                except:
                        print(assembly , "  not found")
                for scaffold in subBlast_table['Scaffold_name'].unique():
                        #print(">",assembly,":",assembly_genome[scaffold].id,sep="",file=output)
                        assembly_genome[scaffold].id=assembly+":"+assembly_genome[scaffold].id
                        assembly_genome[scaffold].description=assembly+":"+assembly_genome[scaffold].description
                        #print(textwrap.fill(str(assembly_genome[scaffold].seq),width=80),file=output)
                        SeqIO.write(assembly_genome[scaffold],output,"fasta")
                Loaded_scaffolds.append(assembly+":"+scaffold)
                n+=1
                print(n,"/",tot)
"""
######


#############################################
# Extract all the EVE AA  sequences        ##
#############################################

# iterate over each group to get EVE loci sequences 

AA_records=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses.faa","fasta"))

List_viruses =['Chelonus_insularis','Trybliographa','Rhoptromeris','Thrichoplasta','Leptopilina_clavipes','Leptopilina_boulardi_GCA_003121605','Leptopilina_heterotoma_GCA_010016045','Leptopilina_heterotoma_GCA_009602685',
'Leptopilina_heterotoma_GCA_009026005','Leptopilina_heterotoma_GCA_009025955','Leptopilina_boulardi_GCA_015476485',
'Leptopilina_boulardi_GCA_011634795','Leptolamina','Phanerotoma','Dolichomitus','Melanaphis_sacchari','Eurytoma_brunniventris','Leptopilina_heterotoma_GCA_015476425.1','Leptopilina_boulardi_GCA_019393585.1','Platygaster_equestris','PoEFV']


with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_hits.aa","w") as output_aa:
  for index, row in Blast_table.iterrows():
    for virus in List_viruses:
      if virus in row['Names']:
        print('>',AA_records[row['Names']].id,sep="",file=output_aa)
        print(AA_records[row['Names']].seq,file=output_aa)


#####################################################
# Find ORFs along the scaffolds with candidate loci #
#####################################################

print("Looking for EVEs with ORFS to replace them by complete ORF with at least 70% of the the size of the best viral hit protein...")

Blast_table['Species_name']=Blast_table['Names'].str.replace(".*:","")
Blast_table['Species_scaffold_name']=Blast_table['Species_name']+":"+Blast_table['Scaffold_name']

#

# Open the fasta file containing all scaffolds with candidate EVEs
All_scaffolds=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.dna", "fasta"))


# Run python ORF finder 
import pathlib  

Candidate_scaffolds="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.dna"
Candidate_scaffolds_ORFs_bed="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.bed"
Candidate_scaffolds_ORFs_dna="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.fna"
Candidate_scaffolds_ORFs_aa="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.aa"
Candidate_scaffolds_ORFs_dna="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.fna"
outdir="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/"

if os.path.exists(Candidate_scaffolds_ORFs_bed):
  print(Candidate_scaffolds_ORFs_bed, "file already exists, opening previous run...")
else:
  subprocess.run("/beegfs/data/bguinet/Bguinet_conda/bin/orfipy " +Candidate_scaffolds+"  --ignore-case  --procs 20 --outdir "+outdir+" --bed "+Candidate_scaffolds_ORFs_bed+" --min 150 --start ATG --dna "+Candidate_scaffolds_ORFs_dna +" --pep "+Candidate_scaffolds_ORFs_aa, shell=True)

subprocess.run(" sed -i 's@).*@)@g' " +Candidate_scaffolds_ORFs_aa, shell=True)
subprocess.run("sed -i 's@ \[@;@g' " +Candidate_scaffolds_ORFs_aa, shell=True)
subprocess.run("sed -i 's@\]@@g' " +Candidate_scaffolds_ORFs_aa, shell=True)

subprocess.run(" sed -i 's@).*@)@g' " +Candidate_scaffolds_ORFs_dna, shell=True)
subprocess.run("sed -i 's@ \[@;@g' " +Candidate_scaffolds_ORFs_dna, shell=True)
subprocess.run("sed -i 's@\]@@g' " +Candidate_scaffolds_ORFs_dna, shell=True)

# Extract and translate the ORFs
ORF_bed=pd.read_csv(Candidate_scaffolds_ORFs_bed,sep="\t",header=None)
ORF_bed.columns=['Scaffold_name','ORF_start','ORF_end','ORF_name','zero','ORF_strand']

ORF_bed['ORF_name2']=ORF_bed['ORF_name'].str.replace(";.*","")
ORF_bed['ORF_name2']=ORF_bed['ORF_name2'].str.replace("ID=","")
ORF_bed['ORF_name']=ORF_bed['ORF_name'].str.replace("ID=","")
ORF_bed['dict_name']=ORF_bed['ORF_name2']+";"+ORF_bed['ORF_start'].astype(str)+"-"+ORF_bed['ORF_end'].astype(str)+"("+ORF_bed['ORF_strand']+")"


#Only keep ORFs around 5000bp from candidates to note keep to much ORFs to blast 
Selected_ORF_name=[]
n_tot=len(Blast_table['Species_scaffold_name'].unique())
n=0


print("Keep only ORFs near EVEs (< 5000bp)...")
for Species_scaffold_name in Blast_table['Species_scaffold_name'].unique():
       subORF_bed=ORF_bed.loc[ORF_bed['Scaffold_name'].eq(Species_scaffold_name)]
       subBlast_table=Blast_table.loc[Blast_table['Species_scaffold_name'].eq(Species_scaffold_name)]
       for index, row in subBlast_table.iterrows():
                min_start=int(row['start'])-5000
                max_start=int(row['end'])+5000
                for loci in subORF_bed.loc[subORF_bed['ORF_start'].lt(max_start) &  subORF_bed['ORF_end'].gt(min_start)]['dict_name'].unique():
                        Selected_ORF_name.append(loci)
       # n+=1
        #print(n, "/", n_tot)
#

Candidate_scaffolds_ORFs_aa_filtred="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds_filtrer.aa"

ORF_dict=to_dict_remove_dups(SeqIO.parse(Candidate_scaffolds_ORFs_aa, "fasta"))

#ORF_dict2=[ORF_dict.pop(x) for x in Selected_ORF_name if x in ORF_dict.keys()]


with open(Candidate_scaffolds_ORFs_aa_filtred,'w') as output:
        for loci in Selected_ORF_name:
                print(">",ORF_dict[loci].id,sep="",file=output)
                print(ORF_dict[loci].seq,file=output)
#

# Run mmseqs between ORFs and previously selected candidates

Filtred_loci="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_hits.aa"
subprocess.run("sed -i 's@ @_@g'  /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_hits.aa", shell=True)

Filtred_loci_db=outdir+"All_candidate_hits_db"

Filtred_loci_vs_Predicted_ORFs_result=outdir+"Filtred_loci_vs_Predicted_ORFs_result"
Filtred_loci_vs_Predicted_ORFs_temp=outdir+"Filtred_loci_vs_Predicted_ORFs_tpm"
Filtred_loci_vs_Predicted_ORFs_table=outdir+"Filtred_loci_vs_Predicted_ORFs_result.m8"

Candidate_scaffolds_ORFs_aa="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.aa"
Predicted_ORFs_db=outdir+"scaffold_orfipy_db"

if os.path.exists(Filtred_loci_vs_Predicted_ORFs_table):
  print(Filtred_loci_vs_Predicted_ORFs_table, "file already exists, opening previous run...")
else:
  subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ Candidate_scaffolds_ORFs_aa_filtred + " "+ Predicted_ORFs_db , shell=True)
  subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs createdb "+ Filtred_loci + " "+ Filtred_loci_db , shell=True)
  subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs search  "+  Filtred_loci_db  +" "+ Predicted_ORFs_db + " "+  Filtred_loci_vs_Predicted_ORFs_result + " "+ Filtred_loci_vs_Predicted_ORFs_temp + " -e 0.000001 --threads 30", shell=True)
  subprocess.run("/beegfs/data/bguinet/TOOLS/mmseqs/bin/mmseqs convertalis  --format-output 'query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,tcov,qcov'  "+  Filtred_loci_db +" "+ Predicted_ORFs_db + " "+  Filtred_loci_vs_Predicted_ORFs_result + " "+ Filtred_loci_vs_Predicted_ORFs_table + " --threads 10", shell=True)


ORF_vs_EVE_table=pd.read_csv(Filtred_loci_vs_Predicted_ORFs_table,sep="\t",header=None)
ORF_vs_EVE_table.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','qlen','tlen','tcov','qcov']
# Keep only self matching hits 
ORF_vs_EVE_table=ORF_vs_EVE_table.loc[ORF_vs_EVE_table['bits'].ge(50)]

#Blast_table = pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_results.tab",sep=";")

ORF_vs_EVE_table['query_EVE_species']=ORF_vs_EVE_table['query'].str.replace(":.*","")
#ORF_vs_EVE_table['query_EVE_species'].mask(ORF_vs_EVE_table['query'].str.contains("PoEFV"),ORF_vs_EVE_table['query'].str.replace(".*_",""),inplace=True)
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['target'].str.replace(".*;","")
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['ORF_start'].str.replace("\\(\\+\\)","")
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['ORF_start'].str.replace("\\(-\\)","")
ORF_vs_EVE_table['ORF_strand']="NA"
ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query'].str.contains("\\+\\)"),'ORF_strand'] ="+"
ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query'].str.contains("-\\)"),'ORF_strand'] ="-"
#ORF_vs_EVE_table['target'].str.replace(".*;","")
ORF_vs_EVE_table['ORF_end']=ORF_vs_EVE_table['ORF_start'].str.replace(".*-","").astype(int)
ORF_vs_EVE_table['ORF_start']=ORF_vs_EVE_table['ORF_start'].str.replace("-.*","").astype(int)
ORF_vs_EVE_table['target_ORF_species']=ORF_vs_EVE_table['target'].str.replace("_ORF.*","")
ORF_vs_EVE_table['target_ORF_species1']=ORF_vs_EVE_table['target_ORF_species'].str.replace(":.*","")
ORF_vs_EVE_table['target_ORF_scaffold']=ORF_vs_EVE_table['target_ORF_species'].str.replace(".*:","")
ORF_vs_EVE_table=ORF_vs_EVE_table.merge(Blast_table[['Scaffold_name','Species_name','Names']], left_on ="query", right_on="Names",how="left")
ORF_vs_EVE_table['target_ORF_species_scaffold']=ORF_vs_EVE_table['target_ORF_species1']+"-"+ORF_vs_EVE_table['target_ORF_scaffold']

ORF_vs_EVE_table['query_EVE_scaffold']=ORF_vs_EVE_table['query'].str.replace(":.*","")
ORF_vs_EVE_table['query_EVE_scaffold']=ORF_vs_EVE_table['query_EVE_scaffold'].str.replace(".*-","")
#ORF_vs_EVE_table['query_EVE_scaffold'].mask(ORF_vs_EVE_table['query'].str.contains("PoEFV") , ORF_vs_EVE_table['query'].str.replace("-.*",""),inplace=True)
#ORF_vs_EVE_table['query_EVE_scaffold'].mask(ORF_vs_EVE_table['query'].str.contains("PoEFV") , ORF_vs_EVE_table['query_EVE_scaffold'].str.split('_').str[:-1].str.join('_'),inplace=True)
ORF_vs_EVE_table['target_ORF_scaffold']=ORF_vs_EVE_table['target'].str.replace(".*:","")
ORF_vs_EVE_table['target_ORF_scaffold']=ORF_vs_EVE_table['target_ORF_scaffold'].str.replace("_ORF.*","")
ORF_vs_EVE_table['query_ORF_species_scaffold']=ORF_vs_EVE_table['Species_name']+'-'+ORF_vs_EVE_table['query_EVE_scaffold']

ORF_vs_EVE_table=ORF_vs_EVE_table.loc[ORF_vs_EVE_table['query_ORF_species_scaffold'] == ORF_vs_EVE_table['target_ORF_species_scaffold']]

# Remove duplicate 
ORF_vs_EVE_table = ORF_vs_EVE_table.drop_duplicates()

# find overlapping ORFs within the same scaffold
ORF_vs_EVE_table.rename(columns={'query': 'ORF_query',
                   'target': 'ORF_target'},
          inplace=True, errors='raise')

Blast_table_ORFs=Blast_table.merge(ORF_vs_EVE_table[['ORF_query','ORF_target','ORF_start','ORF_end','ORF_strand']],left_on="Names",right_on="ORF_query",how="left")


Blast_table_ORFs['Overlapp_ORF_EVEs']= np.nan
Blast_table_ORFs['ORF_perc']= np.nan



for index, row in Blast_table_ORFs.loc[Blast_table_ORFs['ORF_start'].ge(0)].iterrows():
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
                overlapp = "yes-inside"
        elif (int(row['start']) <= int(row['ORF_start'])) & (int(row['end']) >= int(row['ORF_end'])) :
                #print("outside")
                #print(row['start']," : ", row['end'])
                #print(row['ORF_start']," : ", row['ORF_end'])
                overlapp = "yes-outside"
        else:
                overlapp = "no"
        Blast_table_ORFs.loc[Blast_table_ORFs['ORF_query'].eq(row['ORF_query']) & Blast_table_ORFs['ORF_end'].eq(row['ORF_end']) & Blast_table_ORFs['end'].eq(row['end']),"Overlapp_ORF_EVEs"]=overlapp
        if overlapp =="no":
                continue
        else:
                Blast_table_ORFs.loc[Blast_table_ORFs['ORF_query'].eq(row['ORF_query']) & Blast_table_ORFs['ORF_end'].eq(row['ORF_end']) & Blast_table_ORFs['end'].eq(row['end']),"ORF_perc"]= (int(row['ORF_end'])- int(row['ORF_start']))/ (int(row['end'])- int(row['start']))
        #print("\n")
#

# remove non-overlapping ORFS with EVEs 
Blast_table_ORFs.loc[Blast_table_ORFs['Overlapp_ORF_EVEs'].eq("no"),"ORF_start"] = np.nan
Blast_table_ORFs.loc[Blast_table_ORFs['Overlapp_ORF_EVEs'].eq("no"),"ORF_end"] = np.nan

Blast_table_ORFs.loc[Blast_table_ORFs['Overlapp_ORF_EVEs'].eq("yes-outside") & Blast_table_ORFs['ORF_perc'].lt(0.5),"ORF_start"] = np.nan
Blast_table_ORFs.loc[Blast_table_ORFs['Overlapp_ORF_EVEs'].eq("yes-outside") & Blast_table_ORFs['ORF_perc'].lt(0.5),"ORF_end"] = np.nan

# find overlapping ORFs within the same scaffold

#Blast_table_ORFs['full_name']=Blast_table_ORFs['Species_name']+":"+Blast_table_ORFs['target']+":"+Blast_table_ORFs['start'].astype(int).astype(str)

is_overlapped = lambda x: x['ORF_start'] >= x['ORF_end'].shift(fill_value=-1)
Blast_table_ORFs['ORF_overlapp_group'] = Blast_table_ORFs.sort_values(['Names', 'ORF_start', 'ORF_end']) \
                .groupby(['Names'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()
#
#Blast_table_ORFs['count_ORF_overlapp_group'] = Blast_table_ORFs.groupby(Blast_table_ORFs['ORF_overlapp_group']).transform('sum')

Blast_table_ORFs = Blast_table_ORFs.sort_values(['Cluster_hmmer','Names', 'ORF_overlapp_group','ORF_perc'], ascending=[True, False,False,True]) \
    .groupby(['Cluster_hmmer','Names', 'ORF_overlapp_group']).head(1)

Blast_table_ORFs.loc[Blast_table_ORFs['ORF_start'].isna(),"ORF_overlapp_group"] = np.nan

# Calculate ORF_perc_best_hit 

# So first load for each Names, the best hits from the previous blast
Old_tab=pd.read_csv("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/NR_analyse/Blast_file_HSPs_NR.tab",sep=";")
# Keep first best hit info 
Old_tab=Old_tab.loc[~Old_tab['target'].str.contains("EFV")]
Old_tab=Old_tab.sort_values(['evalue', 'bits'], ascending=[True, False])
Old_tab = Old_tab.drop_duplicates(subset = "full_name")
Old_tab.rename(columns={'tlen': 'Best_hit_length', "full_name" : "Names"},
          inplace=True, errors='raise')

Blast_table_ORFs=Blast_table_ORFs.merge(Old_tab[["Names","Best_hit_length"]],on="Names", how="left")
Blast_table_ORFs['Best_hit_ORF_perc']= ((Blast_table_ORFs['ORF_end']-Blast_table_ORFs['ORF_start'])/3)/Blast_table_ORFs['Best_hit_length']

#######
# Remove cluster without virus inside 
# Count number filamentous in each clusters 
list_cynipoidea =['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV']
list_filamentous = ['LbFV','EfFV','PoFV','PcFV','CcFV1','CcFV2','LhFV']
list_virus = ['tom', 'OrNV', 'NlENV', 'GbNV', 'FaENV', 'PoEFV', 'LhFV', 'PoFV',
       'GpSGHV', 'DFV', 'EfFV', 'MdSGHV', 'CcFV2', 'PcFV', 'CcFV1',
       'LbFV', 'VcENV', 'DhNV', 'PmNV', 'AcMNPV', 'CpV', 'LdMNPV', 'ToNV',
       'CiENV', 'MdENV', 'CuniNPV', 'WSSV', 'AmFV',
       'HzNV-1', 'BtENV', 'NeseNPV', 'DnENV', 'TnENV', 'CoBV']
       
count_cluster= Blast_table_ORFs[['Cluster_hmmer','Names2']]
count_cluster = count_cluster.drop_duplicates(subset=['Cluster_hmmer', 'Names2'], keep='first')

count_cluster['Cynip_count'] = count_cluster['Names2'].isin(list_cynipoidea).groupby(count_cluster['Cluster_hmmer']).transform('sum')
count_cluster['Filamentous_count'] = count_cluster['Names2'].isin(list_filamentous).groupby(count_cluster['Cluster_hmmer']).transform('sum')
count_cluster['Virus_count'] = count_cluster['Names2'].isin(list_virus).groupby(count_cluster['Cluster_hmmer']).transform('sum')

count_cluster = count_cluster.drop_duplicates(subset=['Cluster_hmmer'], keep='first')
del Blast_table_ORFs['Cynip_count']
del Blast_table_ORFs['Filamentous_count']
Blast_table_ORFs=Blast_table_ORFs.merge(count_cluster[['Cluster_hmmer','Cynip_count','Filamentous_count','Virus_count']],on="Cluster_hmmer",how="left")


#Correct ORF coordinates for lef-5 with alternative start codons 

Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("JADEYJ010000038.1:4332741-4332950(+):Leptopilina_boulardi_GCA_019393585.1"),'Best_hit_ORF_perc']="na"
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("JADEYJ010000038.1:4332741-4332950(+):Leptopilina_boulardi_GCA_019393585.1"),'ORF_start']=4332741
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("JADEYJ010000038.1:4332741-4332950(+):Leptopilina_boulardi_GCA_019393585.1"),'ORF_end']=4332950
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("JADEYJ010000038.1:4332741-4332950(+):Leptopilina_boulardi_GCA_019393585.1"),'ORF_strand']="+"

Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("C2844038:1642-1860(-):Trybliographa"),'Best_hit_ORF_perc']="na"
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("C2844038:1642-1860(-):Trybliographa"),'ORF_start']=1642
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("C2844038:1642-1860(-):Trybliographa"),'ORF_end']=1860
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("C2844038:1642-1860(-):Trybliographa"),'ORF_strand']="-"  

Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("scaffold_532_2:17437-17637(-):Leptopilina_clavipes"),'Best_hit_ORF_perc']="na"
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("scaffold_532_2:17437-17637(-):Leptopilina_clavipes"),'ORF_start']=17432
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("scaffold_532_2:17437-17637(-):Leptopilina_clavipes"),'ORF_end']=1637
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("scaffold_532_2:17437-17637(-):Leptopilina_clavipes"),'ORF_strand']="-" 

Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("scaffold34593|size2940:1336-1542(-):Rhoptromeris"),'Best_hit_ORF_perc']="na"
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("scaffold34593|size2940:1336-1542(-):Rhoptromeris"),'ORF_start']=1334
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("scaffold34593|size2940:1336-1542(-):Rhoptromeris"),'ORF_end']=1641
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].eq("scaffold34593|size2940:1336-1542(-):Rhoptromeris"),'ORF_strand']="-" 


Blast_table_ORFs.to_csv(output_file,sep=";",index=False)

print ("ALL ORFs included and saved into the tab as : ", output)


print ("\n")

print("Writting all the cluster files...")

Blast_table_ORFs=pd.read_csv(output_file,sep=";")


Blast_table_ORFs=Blast_table_ORFs.loc[~Blast_table_ORFs['Species_name'].isin(['Leptopilina_boulardi_GCA_015476485','Leptopilina_boulardi_GCA_003121605',
       'Leptopilina_boulardi_GCA_011634795','Leptopilina_heterotoma_GCA_009602685',
       'Leptopilina_heterotoma_GCA_010016045',
       'Leptopilina_heterotoma_GCA_009026005',
       'Leptopilina_heterotoma_GCA_009025955'])]

#Blast_table_ORFs.loc[Blast_table_ORFs['Species_name'].str.contains("GCA_019393585"),"Species_name"]
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].str.contains("MdBV"),"Names2"]="MdENV"
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].str.contains("CcBV"),"Names2"]="CcENV"
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].str.contains("FaBV"),"Names2"]="FaENV"
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].str.contains("NlBV"),"Names2"]="NlENV"
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].str.contains("ApBV"),"Names2"]="ApENV"
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].str.contains("VcBV"),"Names2"]="VcENV"
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].str.contains("BtBV"),"Names2"]="BtENV"
Blast_table_ORFs.loc[Blast_table_ORFs['Names'].str.contains("DnBV"),"Names2"]="DnENV"


"""
this step is noew done in  Snakemake_phylogenies_for_plot_part7

#############################################
# Write the AA clusters                    ##
 #############################################

print("Writting AA cluster files...")
# iterate over each group to get EVE loci sequences 

AA_records_origin=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses.faa","fasta"))
AA_records_ORFs=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.aa","fasta"))
DNA_records_ORFs=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.fna","fasta"))

Blast_table_ORFS_grouped = Blast_table_ORFs.groupby('Cluster_hmmer')


list_cynipoidea = ['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV']
list_filamentous = ['LbFV','EfFV','PoFV','PcFV','CcFV1','CcFV2','LhFV']

Added_sequences1=[]
Added_sequences2=[]
for group_name, df_group in Blast_table_ORFS_grouped:
    # if there are viruses and more than two sequence within a cluster, create it !
    if df_group['Virus_count'].iloc[0]>0:
      if len(df_group['Names'].unique())>=2:
	#print("writting cluster : ", group_name)
        # EVE protein loci 
        # ALL sequences
        with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+".aa","w") as output_aa:
                for row_index, row in df_group.iterrows():
                  if row['Names'] not in Added_sequences1:
                    if row['Best_hit_ORF_perc']> 0.70:
                      New_names= row['Scaffold_name']+":"+str(int(row['ORF_start']))+"-"+str(int(row['ORF_end']))+"("+str(row['ORF_strand'])+"):"+row['Species_name']+' [ORF]'
                      print('>',New_names,sep="",file=output_aa)
                      print(AA_records_ORFs[row['ORF_target']].seq,file=output_aa)
                      Added_sequences1.append(row['Names'])
                    else:
                      print('>',AA_records_origin[row['Names']].id,sep="",file=output_aa)
                      print(AA_records_origin[row['Names']].seq,file=output_aa)
                      Added_sequences1.append(row['Names'])
                  else:
                    print("")
        # Only Cynipids and Filamentoviridae virus sequences 
        if pd.isnull(df_group['Prot_name'].iloc[0]):
          filename="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_dNdS.aa"
        else:
          prot_name=df_group['Prot_name'].iloc[0]
          prot_name=re.sub(" \\(","-",prot_name)
          prot_name=re.sub("\\)","",prot_name)
          filename="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_dNdS_"+prot_name+".aa"
        with open(str(filename),"w") as output_aa_dNdS:
                for row_index, row in df_group.iterrows():
                  if row['Filamentous_count']>0:
                    if row['Names'] not in Added_sequences2:
                      if len(df_group['Names'].unique())>=2:
                        if row['Cynip_count']>1:
                          if row['Names2'] in list_cynipoidea or  row['Names2'] in list_filamentous:
                            if row['Best_hit_ORF_perc']> 0.70:
                              New_names= row['Scaffold_name']+":"+str(int(row['ORF_start']))+"-"+str(int(row['ORF_end']))+"("+str(row['ORF_strand'])+"):"+row['Species_name']+' [ORF]'
                              print('>',New_names,sep="",file=output_aa_dNdS)
                              print(AA_records_ORFs[row['ORF_target']].seq,file=output_aa_dNdS)
                              Added_sequences2.append(row['Names'])
                            else:
                              print('>',AA_records_origin[row['Names']].id,sep="",file=output_aa_dNdS)
                              print(AA_records_origin[row['Names']].seq,file=output_aa_dNdS)
                              Added_sequences2.append(row['Names'])
                    else: 
                      print("")
#

print ("ALL AA files with all sequences written to : /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/ with extension .aa")
print ("ALL AA dN/dS files written to : /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/ with extension _dNdS_*.aa")

"""
#subprocess.run("find /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/ -size 0 -print -delete ", shell=True)      
#


##########################################################################################
# Write the AA clusters for phylogeny and datation with specific taxa                   ##
##########################################################################################

# iterate over each group to get EVE loci sequences 

AA_records_origin=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses.faa","fasta"))
AA_records_ORFs=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.aa","fasta"))

Blast_table_ORFS_grouped = Blast_table_ORFs.groupby('Cluster_hmmer')

list_cynipoidea = ['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV']
list_filamentous = ['LbFV','EfFV','PoFV','PcFV','CcFV1','CcFV2','LhFV','DoEFV']
list_virus = ['tom', 'OrNV', 'NlENV', 'GbNV', 'FaENV', 'PoEFV', 'LhFV', 'PoFV',
       'GpSGHV', 'DFV', 'EfFV', 'MdSGHV', 'CcFV2', 'PcFV', 'CcFV1',
       'LbFV', 'VcENV', 'DhNV', 'PmNV', 'AcMNPV', 'CpV', 'LdMNPV', 'ToNV',
       'CiENV', 'MdENV', 'CuniNPV', 'WSSV', 'AmFV',
       'HzNV-1', 'BtENV', 'NeseNPV', 'DnENV', 'TnENV', 'CoBV']


list_to_keep=['FaENV', 'OrNV', 'tom', 'NlENV', 'GbNV', 'PoEFV', 'BtENV', 'VcENV',
       'MsENV', 'WSSV', 'CoBV', 'LdMNPV', 'AcMNPV', 'CpV', 'DFV', 'LhFV',
       'EfFV', 'GpSGHV', 'AmFV', 'ToNV', 'DhNV', 'FaBV', 'HzNV-1',
       'CcFV1', 'NeseNPV', 'DoEFV', 'PcFV',
       'PmNV', 'CcFV2', 'MdSGHV', 'RhEFV', 'LhEFV', 'LbEFV', 'TrEFV', 'PoFV', 'LcEFV',
       'CuniNPV', 'LbFV', 'ThEFV', 'DnENV', 'PhENV', 'ApENV',
       'MdENV', 'Chelonus_insularis', 'EbENV', 'CiENV', 'TnENV', 'CcENV',
       'MdBV', 'CcBV', 'NlBV', 'VcBV', 'BtENV','CoENV', 'MdENV', 'FaENV', 'NlENV', 'ApENV', 'VcENV', 'BtENV', 'CcENV','DnENV']

All_species=[]
for i in list_cynipoidea:
	All_species.append(i)
for i in list_filamentous:
	All_species.append(i)
for i in list_virus:
	 All_species.append(i)


# To be a cluster for the phylogeny, there has to be at least 4 loci within the cluster :

for group_name, df_group in Blast_table_ORFs.groupby('Cluster_hmmer'):
        sub_df_group=df_group.loc[df_group['Names2'].isin(list_to_keep)]
        #Keep longest sequence per paralogs 
        sub_df_group['Length']=sub_df_group['end']-sub_df_group['start']
        sub_df_group=sub_df_group.sort_values(['Names2', 'Length'], ascending=[True, False])
        sub_df_group=sub_df_group.drop_duplicates(subset ="Names2",keep = "first") 
        if len(sub_df_group['Names2'].unique())>=4:
        # EVE protein loci 
        # ALL sequences
          with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_phylogeny.aa","w") as output_aa:
                for row_index, row in sub_df_group.iterrows():
                    if row['Best_hit_ORF_perc']> 0.70:
                      print('>',row['Names2'],sep="",file=output_aa)
                      print(AA_records_ORFs[row['ORF_target']].seq,file=output_aa)
                    else:
                      print('>',row['Names2'],sep="",file=output_aa)
                      print(AA_records_origin[row['Names']].seq,file=output_aa)


# To be a cluster for the phylogeny datation , there has to be at least 4 loci within the cluster excluding AmFV, DoEFV since we do not know their real taxonomic assignation in the phylogeny of in their free-living/integrated nature 

for group_name, df_group in Blast_table_ORFs.groupby('Cluster_hmmer'):
       sub_df_group=df_group.loc[df_group['Names2'].isin(list_to_keep)]
       #Keep longest sequence per paralogs 
       sub_df_group['Length']=sub_df_group['end']-sub_df_group['start']
       sub_df_group=sub_df_group.sort_values(['Names2', 'Length'], ascending=[True, False])
       sub_df_group=sub_df_group.drop_duplicates(subset ="Names2",keep = "first") 
       if (len(sub_df_group['Names2'].unique())>=10) or (len(sub_df_group.loc[sub_df_group['Names2'].isin(['CoBV','WSSV'])])>0):
        if len(sub_df_group['Names2'].unique())>=4:
        # EVE protein loci 
        # ALL sequences
          with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_datation.aa","w") as output_aa:
                for row_index, row in sub_df_group.iterrows():
                    if row['Names2'] in ['AmFV','DoEFV','DoFV']:
                       continue
                    else:
                       if row['Best_hit_ORF_perc']> 0.70:
                          print('>',row['Names2'],sep="",file=output_aa)
                          print(AA_records_ORFs[row['ORF_target']].seq,file=output_aa)
                       else:
                          print('>',row['Names2'],sep="",file=output_aa)
                          print(AA_records_origin[row['Names']].seq,file=output_aa)

print ("ALL AA files for phylogeny written to : /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/ with extension .aa")


print("Writting DNA cluster files...")
#
#############################################
# Write the DNA clusters                   ##
#############################################

# iterate over each group to get EVE loci sequences 

List_viruses =['Chelonus_insularis','Trybliographa','Rhoptromeris','Thrichoplasta','Leptopilina_clavipes','Leptopilina_boulardi_GCA_003121605','Leptopilina_heterotoma_GCA_010016045','Leptopilina_heterotoma_GCA_009602685',
'Leptopilina_heterotoma_GCA_009026005','Leptopilina_heterotoma_GCA_009025955','Leptopilina_boulardi_GCA_015476485',
'Leptopilina_boulardi_GCA_011634795','Leptolamina','Phanerotoma','Dolichomitus','Melanaphis_sacchari','Eurytoma_brunniventris','Leptopilina_heterotoma_GCA_015476425.1','Leptopilina_boulardi_GCA_019393585.1','Platygaster_equestris','Platygaster_orseoliae']

if os.path.exists("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses.fna"):
  print("file alreay exists")
else:
  for species in List_viruses:
    if os.path.exists("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/run_mmseqs2/Fasta_viral_loci.fna"): 
      subprocess.run(" cat /beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species+"/run_mmseqs2/Fasta_viral_loci.fna >> /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses.fna", shell="True")

DNA_records_origin=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ALL_candidate_loci_filtred_and_viruses.fna","fasta"))
DNA_records_ORFs=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/ORF_analysis/All_candidate_scaffolds.fna","fasta"))


list_ran=[]

for row_index, row in Blast_table_ORFs.iterrows():
	try:
		prot_name=row['Prot_name']
		prot_name=re.sub(" \\(","-",prot_name)
		prot_name=re.sub("\\)","",prot_name)
		if os.path.exists("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+str(row['Cluster_hmmer'])+"_dNdS_"+str(prot_name)+".dna"):
			list_ran.append(row['Cluster_hmmer'])
	except:
		if os.path.exists("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+str(row['Cluster_hmmer'])+"_dNdS.dna"):
			list_ran.append(row['Cluster_hmmer'])

Blast_table_ORFS_grouped = Blast_table_ORFs.loc[~Blast_table_ORFs['Cluster_hmmer'].isin(list_ran)].groupby('Cluster_hmmer')

list_cynipoidea = ['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV']
list_filamentous = ['LbFV','EfFV','PoFV','PcFV','CcFV1','CcFV2','LhFV']

from Bio import Entrez
Entrez.email = "benjamin.guinet95@gmail.com"
Entrez.api_key = "30bf99cff0e43d6827934fa6ab127f3b5f09"
 
Added_sequences1=[]
Added_sequences2=[]
for group_name, df_group in Blast_table_ORFS_grouped:
    # if there are viruses and more than two sequence within a cluster, create it !
    if df_group['Virus_count'].iloc[0]>0:
      if len(df_group['Names'].unique())>=2:
        print(group_name)
        # EVE protein loci 
        # Only Cynipids and Filamentoviridae virus sequences 
        if pd.isnull(df_group['Prot_name'].iloc[0]):
          filename="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_dNdS.dna"
        else:
          prot_name=df_group['Prot_name'].iloc[0]
          prot_name=re.sub(" \\(","-",prot_name)
          prot_name=re.sub("\\)","",prot_name)
          filename="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_dNdS_"+prot_name+".dna"
        with open(str(filename),"w") as output_dna_dNdS:
                for row_index, row in df_group.iterrows():
                 if row['Names'] in Added_sequences1:
                  print("")
                 else:
                  Added_sequences1.append(row['Names'])
                  if row['Filamentous_count']>0:
                    if row['Names'] not in Added_sequences2:
                      if len(df_group['Names'].unique())>=2:
                        if row['Cynip_count']>0:
                          if row['Names2'] in list_cynipoidea or  row['Names2'] in list_filamentous:
                            if row['Names4'] in list_filamentous:
                              if row['Names2'] =="LbFV":
                                  print(row['Names2'])
                                  handle=Entrez.efetch(db="protein", id=re.sub("_LbFV","",row['Names']), rettype="fasta_cds_na", retmode="text")
                                  record = SeqIO.read(handle, "fasta")
                                  print('>',row['Names'],sep="",file=output_dna_dNdS)
                                  print(str(record.seq),file=output_dna_dNdS)
                                  #print('>',row['Names'],sep="",file=output_aa_dNdS)
                                  #print(str(record.seq),file=output_aa_dNdS)
                              else:
                                DNA_records_virus=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Genomes/"+row['Names2']+"/Predicted_orfs/Final_ORF_prediction_"+row['Names2']+".fna","fasta"))
                                print('>',DNA_records_virus[row['Names']].id,sep="",file=output_dna_dNdS)
                                print(DNA_records_virus[row['Names']].seq,file=output_dna_dNdS)
                                Added_sequences2.append(row['Names'])
                            else:
                              if row['Best_hit_ORF_perc']> 0.70:
                                New_names= row['Scaffold_name']+":"+str(int(row['ORF_start']))+"-"+str(int(row['ORF_end']))+"("+str(row['ORF_strand'])+"):"+row['Species_name']+' [ORF]'
                                print('>',New_names,sep="",file=output_dna_dNdS)
                                print(DNA_records_ORFs[row['ORF_target']].seq,file=output_dna_dNdS)
                                Added_sequences2.append(row['Names'])
                              else:
                                try:
                                  print('>',DNA_records_origin[row['Names']].id,sep="",file=output_dna_dNdS)
                                  print(DNA_records_origin[row['Names']].seq,file=output_dna_dNdS)
                                except:
                                   new_DNA_records_origin=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+re.sub(".*:","",row['Names'])+"/assembly/"+re.sub(".*:","",row['Names'])+"_final_assembly.fna","fasta"))
                                   if '+' in row['Names']:
                                    print('>',row['Names'],sep="",file=output_dna_dNdS)
                                    print(str(new_DNA_records_origin[row['Scaffold_name']].seq[int(row['start']-1):int(row['end'])]),file=output_dna_dNdS)
                                   else:
                                    print('>',row['Names'],sep="",file=output_dna_dNdS)
                                    #print(row[['Names','start','end','Scaffold_name']])
                                    print(str(new_DNA_records_origin[row['Scaffold_name']].seq[int(row['start']-1):int(row['end'])].reverse_complement()),file=output_dna_dNdS)
                                Added_sequences2.append(row['Names'])
                    else: 
                      print(group_name, " not created...")
                      continue
#

print ("ALL DNA dN/dS files written to : /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/ with extension _dNdS_*.dna")
subprocess.run("find /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/ -size 0 -print -delete ", shell=True)      



####################################################################
## Create MSA AA fasta file to plot the MSA alignments  ############
####################################################################


from Bio import Entrez
Entrez.email = "benjamin.guinet95@gmail.com"
Entrez.api_key = "30bf99cff0e43d6827934fa6ab127f3b5f09"
 
Added_sequences1=[]
Added_sequences2=[]

list_cynipoidea = ['LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV']
list_filamentous = ['LbFV','EfFV','PoFV','PcFV','CcFV1','CcFV2','LhFV']

Blast_table_ORFs_msa=Blast_table_ORFs.loc[Blast_table_ORFs['Filamentous_count'].gt(0) & Blast_table_ORFs['Cynip_count'].gt(0)]
Blast_table_ORFs_msa=Blast_table_ORFs_msa.loc[Blast_table_ORFs_msa['Names2'].isin(['LbFV','EfFV','PoFV','PcFV','CcFV1','CcFV2','LhFV','LbEFV','LcEFV','LhEFV','RhEFV','ThEFV','TrEFV'])]
Blast_table_ORFs_msa['Names4']=Blast_table_ORFs_msa['Names2']
mask = Blast_table_ORFs_msa[['Cluster_hmmer','Names2']].duplicated(keep=False)
Blast_table_ORFs_msa.loc[mask, 'Names2'] += Blast_table_ORFs_msa.groupby(['Cluster_hmmer','Names2']).cumcount().add(1).astype(str)

Blast_table_ORFS_grouped = Blast_table_ORFs_msa.groupby('Cluster_hmmer')

for group_name, df_group in Blast_table_ORFS_grouped:
    # if there are viruses and more than two sequence within a cluster, create it !
    if df_group['Virus_count'].iloc[0]>0:
      if len(df_group['Names'].unique())>=2:
        # EVE protein loci 
        # Only Cynipids and Filamentoviridae virus sequences 
        if pd.isnull(df_group['Prot_name'].iloc[0]):
          filename="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_MSA.aa"
        else:
          prot_name=df_group['Prot_name'].iloc[0]
          prot_name=re.sub(" \\(","-",prot_name)
          prot_name=re.sub("\\)","",prot_name)
          filename="/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+group_name+"_MSA_"+prot_name+".aa"
        with open(str(filename),"w") as output_dna_dNdS:
                for row_index, row in df_group.iterrows():
                 if row['Names'] in Added_sequences1:
                  print("")
                 else:
                  Added_sequences1.append(row['Names'])
                  if row['Filamentous_count']>0:
                    if row['Names'] not in Added_sequences2:
                      if len(df_group['Names'].unique())>=2:
                        if row['Cynip_count']>0:
                          if row['Names4'] in list_cynipoidea or row['Names4'] in list_filamentous:
                            if row['Names4'] in list_filamentous:
                              if row['Names4'] =="LbFV":
                                  handle=Entrez.efetch(db="protein", id=re.sub("_LbFV","",row['Names']), rettype="fasta_cds_na", retmode="text")
                                  record = SeqIO.read(handle, "fasta")
                                  print('>',row['Names2'],sep="",file=output_dna_dNdS)
                                  print(str(record.seq.translate()),file=output_dna_dNdS)
                              else:
                                DNA_records_virus=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/LbFV_family_project/Genomes/"+row['Names4']+"/Predicted_orfs/Final_ORF_prediction_"+row['Names4']+".fna","fasta"))
                                print('>',row['Names2'],sep="",file=output_dna_dNdS)
                                print(DNA_records_virus[row['Names']].seq.translate(),file=output_dna_dNdS)
                                Added_sequences2.append(row['Names'])
                            else:
                              if row['Best_hit_ORF_perc']> 0.70:
                                New_names=row['Names2']
                                #New_names= row['Scaffold_name']+":"+str(int(row['ORF_start']))+"-"+str(int(row['ORF_end']))+"("+str(row['ORF_strand'])+"):"+row['Species_name']+' [ORF]'
                                print('>',New_names,sep="",file=output_dna_dNdS)
                                print(DNA_records_ORFs[row['ORF_target']].seq.translate(),file=output_dna_dNdS)
                                Added_sequences2.append(row['Names'])
                              else:
                                try:
                                  print('>',row['Names2'],sep="",file=output_dna_dNdS)
                                  print(DNA_records_origin[row['Names']].seq.translate(),file=output_dna_dNdS)
                                except:
                                   new_DNA_records_origin=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+re.sub(".*:","",row['Names'])+"/assembly/"+re.sub(".*:","",row['Names'])+"_final_assembly.fna","fasta"))
                                   if '+' in row['Names']:
                                    print('>',row['Names2'],sep="",file=output_dna_dNdS)
                                    print(str(new_DNA_records_origin[row['Scaffold_name']].seq[int(row['start']-1):int(row['end'])].translate()),file=output_dna_dNdS)
                                   else:
                                    print('>',row['Names2'],sep="",file=output_dna_dNdS)
                                    #print(row[['Names','start','end','Scaffold_name']])
                                    print(str(new_DNA_records_origin[row['Scaffold_name']].seq[int(row['start']-1):int(row['end'])].reverse_complement().translate()),file=output_dna_dNdS)
                                Added_sequences2.append(row['Names'])
                    else: 
                      print(group_name, " not created...")
                      continue
subprocess.run(" cp /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/*MSA* /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/All_MSA_files/ ", shell=True)   


from Bio import SeqIO 
import pandas as pd 
import re,os 

for fasta_file in os.listdir('/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/All_MSA_files/'):
  if 'MSA' in fasta_file:
    Cluster_name=re.sub(".*\\/","",fasta_file)
    Cluster_name=re.sub(".aa","",Cluster_name)
    record_dict=SeqIO.to_dict(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+fasta_file, "fasta"))
    list_record=[]
    for record in SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/"+fasta_file, "fasta"):
      list_record.append(record.id)
    list_order=['TrEFV','ThEFV','RhEFV','RhEFV1','RhEFV2','RhEFV3','RhEFV4','LhEFV','LcEFV','LcEFV1','LcEFV2','LbEFV','LbEFV1','LbEFV2','LbEFV3','LbEFV4','LbEFV5','LbEFV6','LbEFV7','LbEFV8','LbEFV9','LbEFV10','LbEFV11','LbEFV12','LbEFV13','LbEFV14','LbEFV15','LbEFV16',
    'LbEFV17','LbEFV18','LbEFV19','LbEFV20','LbEFV21','LbEFV22','LbEFV23','LbEFV24','LbEFV25','LbEFV26','LbEFV27','LbEFV28',
    'LbEFV29','LbEFV30','LbEFV31','LbEFV32','LbEFV33','LbEFV34','LbEFV35','LbEFV36','LbEFV37','LbEFV38','LbEFV39','LbEFV40',
    'LhFV','LhFV1','LhFV2','LbFV','LbFV1','LbFV2','DFV','EfFV','EfFV1','EfFV2','PcFV','PcFV1','PcFV2','PoFV','PoFV1','PoFV2','PoFV3','CcFV1','CcFV2','CcFV21','CcFV22',
    'CcFV11','CcFV12']
    list_record = sorted(list_record, key = lambda val : list_order.index(val))
    with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Cluster_final/All_MSA_files/"+Cluster_name+".aa","w") as output :
      for record in list_record:
        print('>',record_dict[record].id,sep="",file=output)
        print(record_dict[record].seq,file=output)

"""
# Save all the scaffold containing the actual EVEs for the next analyses 

############################################
# Extract all scaffold containing candidates #
############################################

Blast_table=pd.read_csv(output_file,sep=";")
Loaded_scaffolds=[]
if os.path.exists("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/All_candidate_scaffolds.dna") :
        Loaded_scaffolds_dic= to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/All_candidate_scaffolds.dna","fasta"))
        Loaded_scaffolds=list(Loaded_scaffolds_dic.keys())

# Manually add scaffolds 

Blast_table.loc[Blast_table['Names'].str.contains("NW_025111203.1"),"Scaffold_name"]="NW_025111203.1"
Blast_table.loc[Blast_table['Names'].str.contains("JADEYJ010000123.1"),"Scaffold_name"]="JADEYJ010000123.1"

tot=len(Blast_table['Species_name'].unique())
n=0
with open ("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/All_candidate_scaffolds.dna","a") as output:
        for assembly in Blast_table.loc[~ (Blast_table['Species_scaffold_name'].isin(Loaded_scaffolds)) & ~ (Blast_table['Scaffold_name'].isna()) ]['Species_name'].unique():
                print(" recovering scaffold of the species : ", assembly)
                subBlast_table=Blast_table.loc[Blast_table['Species_name'].eq(assembly)]
                subBlast_table=subBlast_table.loc[~subBlast_table['Species_name'].isin(Loaded_scaffolds)]
                try:
                        assembly_genome=to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+assembly+"/assembly/"+assembly+"_final_assembly.fna", "fasta"))
                except:
                        print(assembly , "  not found")
                for scaffold in subBlast_table['Scaffold_name'].unique():
                        #print(">",assembly,":",assembly_genome[scaffold].id,sep="",file=output)
                        assembly_genome[scaffold].id=assembly+":"+assembly_genome[scaffold].id
                        assembly_genome[scaffold].description=assembly+":"+assembly_genome[scaffold].description
                        #print(textwrap.fill(str(assembly_genome[scaffold].seq),width=80),file=output)
                        SeqIO.write(assembly_genome[scaffold],output,"fasta")
                Loaded_scaffolds.append(assembly+":"+scaffold)
                n+=1
                print(n,"/",tot)
######

Loaded_scaffolds_dic= to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/All_candidate_scaffolds.dna","fasta"))
with open("/beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/Metaeuk_TE_analysis/All_candidate_scaffolds.dna", 'w') as handle:
    SeqIO.write(Loaded_scaffolds_dic.values(), handle, 'fasta')

"""
