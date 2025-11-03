#!/bin/bash
#$ -l h_vmem=16G
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=1:00:00
#$ -l os=RedHat7
#$ -R y
#$ -cwd

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse
use Python-3.9
use Anaconda3
source activate /broad/vangalenlab/safina/tools/conda/envs/scenic

python prepare_input_loom_file.py
