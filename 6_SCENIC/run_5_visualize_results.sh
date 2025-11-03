#!/bin/bash
#$ -l h_vmem=64G
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=2:00:00
#$ -l os=RedHat7
#$ -R y
#$ -cwd

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse
use Python-3.9
use Anaconda3
source activate /broad/vangalenlab/safina/tools/conda/envs/scenic

mkdir figures
mkdir final_loom_files


for i in Zhang2022 Hourigan2018 Safina2024 Ainciburu2023 Weng2024 Li2024 Adelman2019
do
  echo $i
  python plot_scenic_results.py $i
done
