#!/bin/bash
#$ -l h_vmem=64G
#$ -pe smp 2
#$ -binding linear:2
#$ -l h_rt=15:00:00
#$ -l os=RedHat7
#$ -R y
#$ -cwd

# for i in Ainciburu2023 Li2024 Weng2024 Adelman2019 Hourigan2018 Safina2024 Zhang2022; do qsub -v dataset=$i run_2_grnboost.sh; done
# for i in Ainciburu2023 Li2024 Weng2024; do qsub -v dataset=$i run_grnboost.sh; done

# Ainciburu run stats: 82.349G,wallclock~4h,cpu~2cores. Hojun: 37.515G,wallclock~5h,cpu~2cores. Chen: 52.147G,~9h,~2 cores

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse
use Python-3.9
use Anaconda3
source activate /broad/vangalenlab/safina/tools/conda/envs/scenic

mkdir grnboost2_out

arboreto_with_multiprocessing.py --num_worker 2 --method grnboost2 -o grnboost2_out/${dataset}.csv input_loom_objects/${dataset}.seu.loom scenic_db/allTFs_hg38.txt
