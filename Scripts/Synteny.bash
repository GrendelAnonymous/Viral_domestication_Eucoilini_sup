#!/bin/bash
#SBATCH --mem 10G
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=1
#SBATCH -e /beegfs/data/bguinet/Cynipoidea_project/Scripts/Synteny_job.error
#SBATCH -o /beegfs/data/bguinet/Cynipoidea_project/Scripts/Synteny_job.out
#SBATCH -J Synteny job

python3 -u /beegfs/data/bguinet/Cynipoidea_project/Scripts/Synteny.py
