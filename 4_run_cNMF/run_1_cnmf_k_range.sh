#!/bin/bash
#$ -l h_vmem=1G
#$ -pe smp 16
#$ -binding linear:16
#$ -l h_rt=2:00:00
#$ -l os=RedHat7
#$ -R y
#$ -cwd
 #
# qsub launch
# for i in `cat sample_list.txt`; do qsub -v sample=$i,lowestk=3,highestk=10 run_1_cnmf_k_range.sh; done

source /broad/software/scripts/useuse
use Anaconda3
source activate /broad/vangalenlab/safina/tools/conda/envs/cnmf

echo $sample

mkdir cNMF_out
cd cNMF_out

rm -r $sample
mkdir $sample
cd $sample
python /broad/vangalenlab/safina/projects/aging_hsc/prepare_for_publication/4_run_cNMF/run_cnmf_k_range.py $sample $lowestk $highestk

use R-4.1
use .gcc-13.2.0

# Plot usage heatmaps
cd sample
mkdir postplots
python /broad/vangalenlab/safina/projects/aging_hsc/prepare_for_publication/4_run_cNMF/plot_heatmap.py
