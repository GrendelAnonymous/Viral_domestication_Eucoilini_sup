# -*- coding:utf-8 -*-

import pandas as pd 
import numpy as np
import argparse
import re,os
import subprocess
from taxadb.taxid import TaxID
from Bio import SeqIO 


# Print out a message when the program is initiated.
print('----------------------------------------------------------------\n')
print('                        Filter output from Blast .\n')
print('----------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Allow add strand information and change coordinate position in the blast output file')
parser.add_argument("-f", "--file", help=" the mmseqs blast file in .m8 format")
parser.add_argument("-of", "--output_file", help=" the desired output filtred mmseqs blast file")
parser.add_argument("-oaa", "--output_aa", help=" the desired output filtred mmseqs blast file")
parser.add_argument("-odna", "--output_dna", help=" the desired output filtred mmseqs blast file")
parser.add_argument("-g", "--genome_assembly", help=" the genome assembly path")
parser.add_argument("-sp", "--species_name", help=" the name of the species")

args = parser.parse_args()

file =args.file
output_file= args.output_file
output_aa=args.output_aa
output_dna=args.output_dna
Genome_assembly=args.genome_assembly
species_name=args.species_name
"""
python3 Filter_Extract_and_translated_blast_loci.py -f /beegfs/data/bguinet/Cynipoidea_project/Genomes/Leptolamina/run_mmseqs2/Mmseqs_dsDNA_EVEs_results.m8 -of /beegfs/data/bguinet/Cynipoidea_project/Genomes/Leptolamina/run_mmseqs2/Mmseqs_dsDNA_EVEs_results_filtred.m8 \
  -oaa /beegfs/data/bguinet/Cynipoidea_project/Genomes/Leptolamina/run_mmseqs2/Fasta_viral_loci.faa \
  -odna /beegfs/data/bguinet/Cynipoidea_project/Genomes/Leptolamina/run_mmseqs2/Fasta_viral_loci.fna \
  -g /beegfs/data/bguinet/Cynipoidea_project/Genomes/Rhoptromeris/assembly/Rhoptromeris_final_assembly.fna \
  -sp Rhoptromeris

"""
"""
species_name="Chelonus_insularis"
species_name="Rhoptromeris"
file="/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species_name+"/run_mmseqs2/Mmseqs_dsDNA_EVEs_results.m8"
output_file="/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species_name+"/run_mmseqs2/Mmseqs_dsDNA_EVEs_results_filtred.m8"
Genome_assembly="/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species_name+"/assembly/Rhoptromeris_final_assembly.fna"
output_dna="/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species_name+"/run_mmseqs2/Fasta_viral_loci.fna"
output_aa="/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+species_name+"/run_mmseqs2/Fasta_viral_loci.faa"
"""

Blast_table=pd.read_csv(file,sep="\t",header=None)
Blast_table.columns=['query','target','pident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue','bits','tlen','qlen','qcov','tcov']
Blast_table.reset_index(drop=True, inplace=True)

# Filter 
Blast_table=Blast_table.loc[Blast_table['bits'].ge(40)]

Blast_table['new_start']=np.nan
Blast_table['new_end']=np.nan
Blast_table['strand']="+"
Blast_table.loc[Blast_table['qstart']> Blast_table['qend'],"strand"]="-"
m = Blast_table['strand'].eq('-')
Blast_table[['qstart','qend']] = np.where(m.to_numpy()[:, None], 
                                       Blast_table[['qend','qstart']], 
                                       Blast_table[['qstart','qend']])
#


# Deal with HSPs 

Blast_table['count_duplicate'] = Blast_table.groupby(['target','query'])['target'].transform('size')

Putative_HSPs=Blast_table.loc[Blast_table['count_duplicate'].ge(2)]

#Run Augustus on those scaffolds 

Augustus_scaffolds=Putative_HSPs['query'].unique()

Putative_HSPs=Putative_HSPs.sort_values(['target','query'], ascending=[True, True])
Putative_HSPs['diff_length_target']=np.nan
Putative_HSPs['diff_length_query']=np.nan
Putative_HSPs['HSP_group']=np.nan



Putative_HSPs = Putative_HSPs.reset_index(drop=True) 

Putative_HSPs['HSP_group']=Putative_HSPs.groupby([ 'query','target']).ngroup()
Putative_HSPs['HSP_min_target']=Putative_HSPs.groupby('HSP_group')["tend"].transform("min")
Putative_HSPs['HSP_max_target']=Putative_HSPs.groupby('HSP_group')["tstart"].transform("max")
Putative_HSPs['HSP_min_query']=Putative_HSPs.groupby('HSP_group')["qend"].transform("min")
Putative_HSPs['HSP_max_query']=Putative_HSPs.groupby('HSP_group')["qstart"].transform("max")
Putative_HSPs['diff_length_target']=abs(Putative_HSPs['HSP_min_target']-Putative_HSPs['HSP_max_target'])
Putative_HSPs['diff_length_query']=abs(Putative_HSPs['HSP_min_query']-Putative_HSPs['HSP_max_query'])



is_overlapped = lambda x: x['qstart'] >= x['qend'].shift(fill_value=-1)
Putative_HSPs['Overlapp_HSP_query_group'] = Putative_HSPs.sort_values(['query','target','HSP_group', 'qstart', 'qend']) \
                .groupby(['query','target','HSP_group'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()

###


Putative_HSPs["qstart"] = Putative_HSPs.groupby(['query', 'target','Overlapp_HSP_query_group'])["qstart"].transform("min")
Putative_HSPs["qend"] = Putative_HSPs.groupby(['query', 'target','Overlapp_HSP_query_group'])["qend"].transform("max")
Putative_HSPs["tstart"] = Putative_HSPs.groupby(['query', 'target','Overlapp_HSP_query_group'])["tstart"].transform("min")
Putative_HSPs["tend"] = Putative_HSPs.groupby(['query', 'target','Overlapp_HSP_query_group'])["tend"].transform("max")
Putative_HSPs["bits"] = Putative_HSPs.groupby(['query', 'target','Overlapp_HSP_query_group'])['bits'].transform('mean') #ok
Putative_HSPs["tcov"] = Putative_HSPs.groupby(['query', 'target','Overlapp_HSP_query_group'])['tcov'].transform('sum')


Putative_HSPs = Putative_HSPs.drop_duplicates(subset=['target','query','Overlapp_HSP_query_group'], keep='first')

Putative_HSPs['HSP_min_query']=Putative_HSPs.groupby('HSP_group')["qend"].transform("min")
Putative_HSPs['HSP_max_query']=Putative_HSPs.groupby('HSP_group')["qstart"].transform("max")
Putative_HSPs['HSP_min_target']=Putative_HSPs.groupby('HSP_group')["tend"].transform("min")
Putative_HSPs['HSP_max_target']=Putative_HSPs.groupby('HSP_group')["tstart"].transform("max")
Putative_HSPs['diff_length_target']=abs(Putative_HSPs['HSP_max_target']-Putative_HSPs['HSP_min_target'])
Putative_HSPs['diff_length_query']=abs(Putative_HSPs['HSP_max_query']-Putative_HSPs['HSP_min_query'])

#Group close qstart-qend coordinates < 5OObp away from each other to merge them
Putative_HSPs=Putative_HSPs.sort_values(['target','query','HSP_group', 'qstart'], ascending=[True, True,True,True])

def get_group(g):
    s = g['qend'].shift().rsub(g['qstart']).gt(5000)
    return s.cumsum().add(1-s.iloc[0]).astype(str).radd('G')

Putative_HSPs['HSP_close_query_groups'] = (Putative_HSPs
   .groupby(['target','query','HSP_group'], group_keys=False)
   .apply(get_group)
)


def get_group(g):
    s = g['tend'].shift().rsub(g['tstart']).gt(1000)
    return s.cumsum().add(1-s.iloc[0]).astype(str).radd('G')

Putative_HSPs['HSP_close_query_target_groups'] = (Putative_HSPs
   .groupby(['target','query','HSP_group','HSP_close_query_groups'], group_keys=False)
   .apply(get_group)
)

Putative_HSPs['HSP_close_target_groups'] = (Putative_HSPs
   .groupby(['target','query','HSP_group'], group_keys=False)
   .apply(get_group)
)


Putative_HSPs=Putative_HSPs.sort_values(['query', 'target','tstart','tend'], ascending=[True, True,True,True])

is_overlapped = lambda x: x['tstart'] > x['tend'].shift(fill_value=-1)
Putative_HSPs['Overlapp_HSP_target_group'] = Putative_HSPs.sort_values(['query','target','HSP_group', 'tstart', 'tend']) \
                .groupby(['query','target','HSP_group'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()

Putative_HSPs['Nb_overlapping_target'] = (Putative_HSPs['tend'].sub(Putative_HSPs.groupby(['HSP_group', 'Overlapp_HSP_target_group'])['tstart'].shift(-1))
                        .ffill(downcast='infer'))

Putative_HSPs=Putative_HSPs.sort_values(['query', 'target','qstart','qend'], ascending=[True, True,True,True])
Putative_HSPs['Nb_overlapping_query'] = (Putative_HSPs['qend'].sub(Putative_HSPs.groupby(['HSP_group', 'Overlapp_HSP_query_group'])['qstart'].shift(-1))
                        .ffill(downcast='infer'))


Putative_HSPs=Putative_HSPs.sort_values(['query', 'target','tstart','tend'], ascending=[True, True,True,True])
Putative_HSPs['Nb_sep_target'] = (Putative_HSPs['tend'].sub(Putative_HSPs.groupby(['HSP_group', 'HSP_close_target_groups'])['tstart'].shift(-1))
                        .ffill(downcast='infer'))


Putative_HSPs=Putative_HSPs.sort_values(['query', 'target','qstart','qend'], ascending=[True, True,True,True])
Putative_HSPs['Nb_sep_query'] = (Putative_HSPs['qend'].sub(Putative_HSPs.groupby(['HSP_group', 'HSP_close_target_groups'])['qstart'].shift(-1))
                        .ffill(downcast='infer'))

#Remove groups where there are overlapping targets 
Putative_HSPs['count_duplicate_target'] = Putative_HSPs.groupby(['target','query','HSP_close_query_groups','Overlapp_HSP_target_group'])['target'].transform('size')
Putative_HSPs['Groups']=Putative_HSPs['HSP_close_query_groups']
Putative_HSPs['Groups']=Putative_HSPs['Groups'].str.replace("G","").astype(int)
for index, row in Putative_HSPs.iterrows():
	if row['Nb_overlapping_target']>=30:
		Putative_HSPs.loc[Putative_HSPs['query'].eq(row['query']) & Putative_HSPs['target'].eq(row['target']) & Putative_HSPs['HSP_group'].eq(row['HSP_group']) & Putative_HSPs['Overlapp_HSP_query_group'].eq(row['Overlapp_HSP_query_group']),"Groups"]= Putative_HSPs['Groups'].max()+1


Putative_HSPs['perc_overlapping']=Putative_HSPs['Nb_overlapping_target']/(Putative_HSPs['tend']-Putative_HSPs['tstart'])
## Change Groups for non cons√cutives

var = False
prv_grp = 0
grp_name = "" # For storing prv. or current  group name

def func(grp, val):
    global grp_name
    global var
    global prv_grp
    if grp == grp_name: # if group hasn't changed
        if var:
            if val < 0.5:
                return "NG"+str(prv_grp)
            else:
                var = False
                prv_grp += 1
                return "NG"+str(prv_grp)
        else:
            if val < 0.5:
                var = True
            prv_grp += 1
            return "NG"+str(prv_grp)
    else: # If group name had changed
        grp_name = grp # Stroring the new Group name
        if val < 0.5:
            var = True
            prv_grp += 1
        else:
            prv_grp += 1
            var = False
        return "NG"+str(prv_grp)
    

#Putative_HSPs['HSP_group_HSP_close_query_groups'] =  Putative_HSPs['HSP_group']+"_"+Putative_HSPs['HSP_close_query_groups']
Putative_HSPs['new_Groups'] = Putative_HSPs.apply(lambda x: func(grp = x.HSP_group, val = x.perc_overlapping), axis=1)


def func(grp, val,val2):
    global grp_name
    global var
    global prv_grp
    if grp == grp_name: # if group hasn't changed
        if var:
            if ((val < 0) and (val2<0)) :
                return "NG"+str(prv_grp)
            else:
                var = False
                prv_grp += 1
                return "NG"+str(prv_grp)
        else:
            if ((val < 0) and (val2<0)) :
                var = True
            prv_grp += 1
            return "NG"+str(prv_grp)
    else: # If group name had changed
        grp_name = grp # Stroring the new Group name
        if ((val < 0) and (val2<0)) :
            var = True
            prv_grp += 1
        else:
            prv_grp += 1
            var = False
        return "NG"+str(prv_grp)

Putative_HSPs['Close_target_and_query_groups'] = Putative_HSPs.apply(lambda x: func(grp = x.HSP_group, val = x.Nb_sep_query, val2 = x.Nb_sep_target), axis=1)


#Putative_HSPs['count_duplicate'] = Putative_HSPs.groupby(['HSP_group','HSP_close_query_groups'])['HSP_close_query_groups'].transform('size')

#Putative_HSPs.loc[Putative_HSPs['count_duplicate'].gt(1)][['query','target','qstart','qend','tstart','tend','HSP_group','HSP_close_query_groups','count_duplicate']]

# Do not merge if the overlap of the target is to high > 50bp 

#is_overlapped = lambda x: x['tstart'] >= x['tend'].shift(fill_value=-1)
#Putative_HSPs['Overlapp_HSP_target_group'] = Putative_HSPs.sort_values(['query','target','HSP_group','HSP_close_query_groups', 'tstart', 'tend']) \
#                .groupby(['query','target','HSP_group','HSP_close_query_groups'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()

Putative_HSPs['HSP_close_query_new_groups'] = (Putative_HSPs
   .groupby(['target','query','HSP_group','new_Groups'], group_keys=False)
   .apply(get_group)
)

# Putative_HSPs.loc[Putative_HSPs['query'].str.contains("6363.1")& Putative_HSPs['target'].str.contains("")][['query','target','qstart','qend','tstart','tend','Overlapp_HSP_query_group','HSP_close_query_groups','HSP_close_query_new_groups','new_Groups','HSP_group']]
#Merge close targets 

Putative_HSPs["qstart"] = Putative_HSPs.groupby(['query', 'target','HSP_close_query_groups','HSP_close_query_new_groups','new_Groups'])["qstart"].transform("min")
Putative_HSPs["qend"] = Putative_HSPs.groupby(['query', 'target','HSP_close_query_groups','HSP_close_query_new_groups','new_Groups'])["qend"].transform("max")
Putative_HSPs["tstart"] = Putative_HSPs.groupby(['query', 'target','HSP_close_query_groups','HSP_close_query_new_groups','new_Groups'])["tstart"].transform("min")
Putative_HSPs["tend"] = Putative_HSPs.groupby(['query', 'target','HSP_close_query_groups','HSP_close_query_new_groups','new_Groups'])["tend"].transform("max")
Putative_HSPs["bits"] = Putative_HSPs.groupby(['query', 'target','HSP_close_query_groups','HSP_close_query_new_groups','new_Groups'])['bits'].transform('mean') #ok
Putative_HSPs["tcov"] = Putative_HSPs.groupby(['query', 'target','HSP_close_query_groups','HSP_close_query_new_groups','new_Groups'])['tcov'].transform('sum')
Putative_HSPs = Putative_HSPs.drop_duplicates(subset=['target','query','HSP_close_query_groups','HSP_close_query_new_groups','new_Groups'], keep='first')


# Merge also if close query and target 

Putative_HSPs["qstart"] = Putative_HSPs.groupby(['query', 'target','HSP_close_query_groups','HSP_close_query_new_groups','HSP_close_query_target_groups','Close_target_and_query_groups'])["qstart"].transform("min")
Putative_HSPs["qend"] = Putative_HSPs.groupby(['query', 'target','HSP_close_query_groups','HSP_close_query_new_groups','HSP_close_query_target_groups','Close_target_and_query_groups'])["qend"].transform("max")
Putative_HSPs["tstart"] = Putative_HSPs.groupby(['query', 'target','HSP_close_query_groups','HSP_close_query_new_groups','HSP_close_query_target_groups','Close_target_and_query_groups'])["tstart"].transform("min")
Putative_HSPs["tend"] = Putative_HSPs.groupby(['query', 'target','HSP_close_query_groups','HSP_close_query_new_groups','HSP_close_query_target_groups','Close_target_and_query_groups'])["tend"].transform("max")
Putative_HSPs["bits"] = Putative_HSPs.groupby(['query', 'target','HSP_close_query_groups','HSP_close_query_new_groups','HSP_close_query_target_groups','Close_target_and_query_groups'])['bits'].transform('mean') #ok
Putative_HSPs["tcov"] = Putative_HSPs.groupby(['query', 'target','HSP_close_query_groups','HSP_close_query_new_groups','HSP_close_query_target_groups','Close_target_and_query_groups'])['tcov'].transform('sum')
Putative_HSPs = Putative_HSPs.drop_duplicates(subset=['target','query','HSP_close_query_groups','HSP_close_query_new_groups','HSP_close_query_target_groups','Close_target_and_query_groups'], keep='first')


##
#Check if hit overlapp within the target between putative HSP that could then be duplicate genes 
is_overlapped = lambda x: x['tstart'] >= x['tend'].shift(fill_value=-1)
Putative_HSPs['Overlapp_HSP_target_group'] = Putative_HSPs.sort_values(['query','target','HSP_group', 'tstart', 'tend']) \
                .groupby(['query','target','HSP_group'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()
#
Putative_HSPs['count_Overlapp_HSP_target_group']=Putative_HSPs.groupby("Overlapp_HSP_target_group")["Overlapp_HSP_target_group"].transform("count")



# Run Augustus for putative intron inside 




#Remove all processed HSPs from the previous tab 
Blast_table=Blast_table.loc[Blast_table['count_duplicate'].lt(2)]

Blast_table2=pd.concat([Blast_table,Putative_HSPs])

Blast_table2 = Blast_table2.reset_index(drop=True)
#Keep the very best hit 
is_overlapped = lambda x: x['qstart'] >= x['qend'].shift(fill_value=-1)
Blast_table2['Overlapp_group'] = Blast_table2.sort_values(['query','qstart','qend']) \
                .groupby(['query'], as_index=False).apply(is_overlapped).droplevel(0).cumsum()

#
Blast_table2['Tot_length'] = abs(Blast_table2['qend'] - Blast_table2['qstart'])
Blast_table2 = Blast_table2.sort_values(['query', 'Tot_length'], ascending=[True, False]) \
        .groupby(Blast_table2['Overlapp_group']).head(1)

Blast_table2['full_name']=Blast_table2['query']+":"+Blast_table2['qstart'].astype(str)+'-'+Blast_table2['qend'].astype(str)+"("+Blast_table2['strand']+")"+":"+species_name


Blast_table2.to_csv(output_file,sep=";",index=False)
print("Blast file filtred printed to : ",output_file)
#

###########################
# Extract candidate loci ##
###########################

#Open records :
print("Opening the genome file...")
Genome_dict= SeqIO.to_dict(SeqIO.parse(Genome_assembly, "fasta"))

#Write new fasta loci aa file 
with open(output_aa,"w") as output_aa_file:
    for index, row in Blast_table2.iterrows():
        seq_name= row['query']+":"+str(row['qstart'])+'-'+str(row['qend'])+"("+row['strand']+")"+":"+species_name
        if row['strand'] =="+":
            print('>',seq_name,sep="",file=output_aa_file)
            print(str(Genome_dict[row['query']].seq[row['qstart']-1:row['qend']].translate()),file=output_aa_file) # since python uses 0-based 
        elif row['strand'] =="-":
            print('>',seq_name,sep="",file=output_aa_file)
            print((Genome_dict[row['query']][row['qstart']-1: row['qend']].seq.reverse_complement().translate()),file=output_aa_file)

#Write new fasta loci dna file 
with open(output_dna,"w") as output_dna_file:
    for index, row in Blast_table2.iterrows():
        seq_name= row['query']+":"+str(row['qstart'])+'-'+str(row['qend'])+"("+row['strand']+")"+":"+species_name
        if row['strand'] =="+":
            print('>',seq_name,sep="",file=output_dna_file)
            print(str(Genome_dict[row['query']].seq[row['qstart']-1:row['qend']]),file=output_dna_file) # since python uses 0-based 
        elif row['strand'] =="-":
            print('>',seq_name,sep="",file=output_dna_file)
            print((Genome_dict[row['query']][row['qstart']-1: row['qend']].seq.reverse_complement()),file=output_dna_file)


print("The translated file has been saved to :", output_aa)
print("Nucleotide file has been saved to :", output_dna)
