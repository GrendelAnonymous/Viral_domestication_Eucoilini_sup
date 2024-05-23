import os
import argparse
import re
import pandas as pd 
from Bio import SeqIO

# Print out a message when the program is initiated.
print('------------------------------------------------------------------------------------------\n')
print('                        Remove duplicate sequence (duplicate name) in cluster ali files .\n')
print('------------------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Remove duplicate sequence within cluster')
parser.add_argument("-c", "--cluster_name", help="The cluster name")
args = parser.parse_args()

#Example : python3 /beegfs/data/bguinet/LbFV_family_project/Clustering/Remove_dup_clusters.py -c Cluster591
Cluster_name=args.cluster_name
Cluster_path="/beegfs/data/bguinet/Cynipoidea_project/Clustering/Cluster_seq_NR/"

def to_dict_remove_dups(sequences):
    return {record.id: record for record in sequences}


#df = pd.DataFrame(columns=['short_name', 'full_name', 'length'])

list_to_keep=[]
#record_seq=SeqIO.parse(Cluster_path+Cluster_name+".faa", "fasta")
for record in SeqIO.parse(Cluster_path+Cluster_name+".faa", "fasta"): 
  list_to_keep.append(record.id)


list_to_keep=list(dict.fromkeys(list_to_keep))

record_dict = to_dict_remove_dups(SeqIO.parse(Cluster_path+Cluster_name+".faa", "fasta"))
with open(Cluster_path+Cluster_name+"_nodup.faa","w") as output:
 for record in list_to_keep:
  print('>',record,sep="",file=output)
  print(record_dict[record].seq,file=output)


os.rename(Cluster_path+Cluster_name+"_nodup.faa",Cluster_path+Cluster_name+".faa") 
print("output written to : ", Cluster_path)





