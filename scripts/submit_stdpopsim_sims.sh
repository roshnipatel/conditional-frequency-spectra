#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --time=4:00:00
#SBATCH --mem=300M
#SBATCH --partition=pritch,owners,hns
#SBATCH --array=0-999%100
#SBATCH -o slurm/slurm-%A-%a.out
#SBATCH -e slurm/slurm-%A-%a.err

model=$1 # jouganous_wo_migration, gutenkunst_wo_migration, ragsdale_w_migration
h=$2 # 0.5, 5e6
s=$3 # 0, -1e-3, -5e-4, -1e-4, -1e-10, -5e-11, -1e-11
dir=data/simulations/$model/h${h}_s$s

source /home/users/rpatel7/.bashrc
conda activate slim
for freq in 0.00003 0.0001 0.0005 0.001 0.002 0.005 0.015 0.025
do
mkdir -p $dir/ancestor_$freq
slim -t -m \
  -d "s=$s" \
  -d "h=$h" \
  -d "freq=$freq" \
  -d "out='$dir/ancestor_$freq/output_$SLURM_ARRAY_TASK_ID.txt'" \
  scripts/slim/$model.slim
done

for i in {1..9}
do
  mkdir -p $dir/ancestor_0.0$i
  slim -t -m \
    -d "s=$s" \
    -d "h=$h" \
    -d "freq=0.0$i" \
    -d "out='$dir/ancestor_0.0$i/output_$SLURM_ARRAY_TASK_ID.txt'" \
    scripts/slim/$model.slim
done

for i in {10..99}
do
  mkdir -p $dir/ancestor_0.$i
  slim -t -m \
    -d "s=$s" \
    -d "h=$h" \
    -d "freq=0.$i" \
    -d "out='$dir/ancestor_0.$i/output_$SLURM_ARRAY_TASK_ID.txt'" \
    scripts/slim/$model.slim
done
conda deactivate
