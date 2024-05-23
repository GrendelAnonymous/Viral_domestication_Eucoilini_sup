import pandas as pd
import argparse
import os
import io 
import os.path
from os import path
from Bio import SeqIO

# Print out a message when the program is initiated.
print('-----------------------------------------------------------------------------------\n')
print('                        Process BUSCO table to create BUSCO files for phylogeny                          .\n')
print('-----------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description=' Process BUSCO table to create BUSCO files for phylogeny')
parser.add_argument("-l", "--sp_list", help="The comma separated list of species to analyse",action='append',nargs='+') 
parser.add_argument("-out_busco_dir", "--out_busco_dir", help="The full busco desired path") 
parser.add_argument("-main_dir", "--main_dir", help="The full Genome paths") 
parser.add_argument("-n_missing", "--n_missing", help="Number of authorized missing busco species per IDs") 
parser.add_argument("-n_max", "--n_max", help="Number max of busco file to create") 
args = parser.parse_args()


#Usage example 

"""
python3 /beegfs/data/bguinet/Cynipoidea_project/Scripts/Create_BUSCO_files_for_dNdS.py -l Leptopilina_boulardi_GCA_019393585.1,Leptopilina_heterotoma_GCA_015476425.1,Leptopilina_clavipes,Rhoptromeris,Thrichoplasta,Trybliographa \
-out_busco_dir /beegfs/data/bguinet/Cynipoidea_project/Analyses/Clustering/BUSCO_dNdS -main_dir /beegfs/data/bguinet/Cynipoidea_project/Genomes/  -n_missing 0  -n_max 1000

"""
# paths 

list_sp=args.sp_list
list_sp=  list_sp[0][0].split(",")

"""
list_sp=["Leptolamina","Leptopilina_boulardi_GCA_015476485","Leptopilina_heterotoma_GCA_009025955","Leptopilina_heterotoma_GCA_010016045","Rhoptromeris","Thrichoplasta",
"Ganaspis_brasiliensis","Leptopilina_boulardi_GCA_003121605","Leptopilina_boulardi_lr","Leptopilina_heterotoma_GCA_009026005","Leptopilina_heterotoma_lr","Synergus_japonicus","Trissolcus_brochymenae",
"Ganaspis_sp","Leptopilina_boulardi_GCA_011634795","Leptopilina_clavipes","Leptopilina_heterotoma_GCA_009602685","Platygaster_equestris","Synergus_umbraculus","Trybliographa"]

list_sp=['Leptopilina_boulardi_GCA_019393585.1','Leptopilina_heterotoma_GCA_015476425.1','Leptopilina_clavipes','Rhoptromeris','Thrichoplasta','Trybliographa']
"""
Main_dir=args.main_dir
#Main_dir="/beegfs/data/bguinet/Cynipoidea_project/Genomes/"
#Output_BUSCO_dir="/beegfs/data/bguinet/Cynipoidea_project/BUSCO_dNdS/"


Output_BUSCO_dir = args.out_busco_dir

print("List of species to analyse:")
print(list_sp)
print("\n")
n_missing=args.n_missing
n_missing=int(n_missing)
#n_missing=0

n_max=args.n_max
n_max=int(n_max)


BUSCO_tab=pd.DataFrame(columns=['Species','Busco_id', 'Status', 'Sequence', 'Gene Start', 'Gene End', 'Strand', 'Score', 'Length', 'OrthoDB url', 'Description'])

# Open all BUSCO tables
for sp in list_sp:
			#print(cat,sp)
			#print(Main_dir,cat,"/results/Stat_Results/",sp,"/BUSCO_results/run_diptera_odb10/full_table.tsv",sep="")
			try:
				#print(Main_dir+sp+"assembly/BUSCO_results_hymenoptera/run_hymenoptera_odb10/full_table.tsv")
				subBUSCO_tab=pd.read_csv(Main_dir+"/"+sp+"/assembly/BUSCO_results_hymenoptera_with_nuc/run_hymenoptera_odb10/full_table.tsv", sep='\t', header=2, index_col=0)
				subBUSCO_tab=subBUSCO_tab.reset_index().rename({'index':'BUSCO_name'}, axis = 'columns')
				subBUSCO_tab.columns=['Busco_id', 'Status', 'Sequence', 'Gene Start', 'Gene End', 'Strand', 'Score', 'Length', 'OrthoDB url', 'Description']
				subBUSCO_tab['Species']=sp
				if sp == 'Leptopilina_heterotoma_GCA_015476425.1':
					subBUSCO_tab=subBUSCO_tab.loc[subBUSCO_tab['Sequence'].astype(str).str.contains("NW_")]
				BUSCO_tab=BUSCO_tab.append(subBUSCO_tab)
			except:
				print(sp)
				continue 

#print(BUSCO_tab)
# Keep only complete BUSCOs
BUSCO_tab=BUSCO_tab.loc[BUSCO_tab['Status'].isin(["Complete"])] # | (BUSCO_tab['Status'].isin(["Complete","Duplicated"]) & BUSCO_tab['Species'].str.contains("hetero"))]

# Count number of BUSCO IDs found in each species 
BUSCO_tab[BUSCO_tab.duplicated(['Busco_id', 'Status', 'Sequence', 'Gene Start', 'Gene End', 'Strand', 'Score', 'Length', 'OrthoDB url', 'Description'])].groupby("Busco_id").size().reset_index(name='Duplicates')
BUSCO_tab_id = BUSCO_tab.groupby('Busco_id')['Species'].nunique().to_frame()

# allow n_missing taxa per BUSCO ids 
BUSCO_tab_id=BUSCO_tab_id.loc[BUSCO_tab_id['Species'].ge(len(list_sp)-n_missing)]
BUSCO_tab_id['Busco_id'] = BUSCO_tab_id.index

print("Number of BUSCO in common on ",len(list_sp)-n_missing,"/",len(list_sp) ," species from the dataset : ", len(BUSCO_tab_id),sep="")

# Create all the retained BUSCO files 
print ("Creating ",len(BUSCO_tab_id),' BUSCO files...',sep="")

BUSCO_tab=BUSCO_tab.loc[BUSCO_tab['Busco_id'].isin(BUSCO_tab_id['Busco_id'])]



#'Leptopilina_boulardi_GCA_019393585.1','Leptopilina_heterotoma_GCA_015476425.1','Leptopilina_clavipes','Rhoptromeris','Thrichoplasta','Trybliographa'
"""
from Bio import SeqIO 
def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}

#Load the fasta assembly files 
for full_sp_name in list_sp:
  if full_sp_name == "Leptopilina_boulardi_GCA_019393585.1":
   record_dict_LbEFV = to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+full_sp_name+"/assembly/"+full_sp_name+"_final_assembly.fna","fasta"))
  if full_sp_name == "Leptopilina_clavipes":
    record_dict_LcEFV = to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+full_sp_name+"/assembly/"+full_sp_name+"_final_assembly.fna","fasta"))
  if full_sp_name == "Leptopilina_heterotoma_GCA_015476425.1":
    record_dict_LhEFV = to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+full_sp_name+"/assembly/"+full_sp_name+"_final_assembly.fna","fasta"))
  if full_sp_name == "Rhoptromeris":
   record_dict_RhEFV = to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+full_sp_name+"/assembly/"+full_sp_name+"_final_assembly.fna","fasta"))
  if full_sp_name == "Thrichoplasta":
    record_dict_ThEFV = to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+full_sp_name+"/assembly/"+full_sp_name+"_final_assembly.fna","fasta"))
  if full_sp_name == "Trybliographa":
    record_dict_TrEFV = to_dict_remove_dups(SeqIO.parse("/beegfs/data/bguinet/Cynipoidea_project/Genomes/"+full_sp_name+"/assembly/"+full_sp_name+"_final_assembly.fna","fasta"))


grouped = BUSCO_tab.groupby('Busco_id')

for group_name, group in grouped:
	print("\n")
	print(group_name)
	for index, row in group.iterrows():
		#Select the good assembly 
		full_sp_name= row['Species']
		if full_sp_name=="Leptopilina_boulardi_GCA_019393585.1":
			record_dict=record_dict_LbEFV
		if full_sp_name=="Leptopilina_heterotoma_GCA_015476425.1":
			record_dict=record_dict_LhEFV
		if full_sp_name=="Leptopilina_clavipes":
			record_dict=record_dict_LcEFV
		if full_sp_name=="Rhoptromeris":
			record_dict=record_dict_RhEFV
		if full_sp_name=="Thrichoplasta":
			record_dict=record_dict_ThEFV
		if full_sp_name=="Trybliographa":
 			record_dict=record_dict_TrEFV
		if full_sp_name=="Trybliographa":
			Scaffold_name=re.sub("\\|.*","",row['Sequence'])
		else:
			Scaffold_name=row['Sequence']
		if row['Strand']=="+":
			seq=str(record_dict[Scaffold_name].seq[int(row['Gene Start']):int(row['Gene End'])].translate())
			#The first exon 
			sequence=re.sub("\\*.*","",seq)
			#print('>',full_sp_name)
			#print(sequence)
		if row['Strand']=="-":
			seq=str(record_dict[Scaffold_name].seq[int(row['Gene Start']):int(row['Gene End'])].reverse_complement().translate())
                        #The first exon 
                        sequence=re.sub("\\*.*","",seq)
                        print('>',full_sp_name)
                        print(sequence)
"""
# Check directory 
if path.exists(Output_BUSCO_dir) ==False:
  os.mkdir(Output_BUSCO_dir)
if path.exists(Output_BUSCO_dir+"/BUSCO_files") ==False:  
  os.mkdir(Output_BUSCO_dir+"/BUSCO_files") 

n=1
for busco in BUSCO_tab_id['Busco_id'].unique():
  if n <= n_max:  
   with open(Output_BUSCO_dir+"/BUSCO_files/BUSCO_"+busco+".fna","w") as output:
    #print(Output_BUSCO_dir+"/BUSCO_files/BUSCO_"+busco+".faa")
    for sp in list_sp :
      try:
        for seq in SeqIO.parse(Main_dir+"/"+sp+"/assembly/BUSCO_results_hymenoptera_with_nuc/run_hymenoptera_odb10/busco_sequences/single_copy_busco_sequences/"+busco+".fna","fasta"):
          print(">",sp,sep="",file=output)
          print(seq.seq,file=output)
      except:
          continue
   n+=1
  else:
    continue

with open(Output_BUSCO_dir+"/BUSCO_files/Control_check1.txt", 'w') as fp:
	print("Rule ok",file=fp) 

print("\n")
print("All BUSCO files written to : ",Output_BUSCO_dir+"/BUSCO_files")
print("Total number of BUSCO file written : ", n-1 ) 




