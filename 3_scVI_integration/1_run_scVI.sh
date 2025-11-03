#!/bin/bash
#$ -l h_vmem=24G
#$ -pe smp 1
#$ -binding linear:1
#$ -l os=RedHat7
#$ -l h_rt=10:00:00
#$ -R y
#$ -cwd

source /broad/software/scripts/useuse
use Anaconda3
use GCC-5.2
source activate /broad/vangalenlab/safina/tools/conda/envs/scvi-env-new
export LD_LIBRARY_PATH=LD_LIBRARY_PATH:/broad/vangalenlab/safina/tools/conda/scvi-env/lib/python3.9/site-packages/nvidia/nvjitlink/lib/

mkdir hsc_integration
cd hsc_integration
python ../integrate_HSCs_with_scVI.py

mkdir hspc_integration
cd hspc_integration
python ../integrate_HSPCs_with_scVI.py
