#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --time=24:00:00
#SBATCH --mem=5G
#SBATCH --partition=pritch,owners,hns

./sm_script.sh $1
