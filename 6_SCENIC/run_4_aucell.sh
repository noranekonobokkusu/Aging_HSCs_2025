#!/bin/bash
#$ -l h_vmem=24G
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_rt=1:00:00
#$ -l os=RedHat7
#$ -cwd
#$ -R y

# for i in Adelman2019 Hourigan2018 Safina2024 Zhang2022 Ainciburu2023 Li2024 Weng2024; do qsub -v dataset=$i run_4_aucell.sh; done

# !!! it doesn't work properly if it doesn't have enough memory, producing ambiguous warnings
# maxvmem<24G, cpu/wallclock<4 cores, wallclock<3 minutes


# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse
use Python-3.9
use Anaconda3
source activate /broad/vangalenlab/safina/tools/conda/envs/scenic


mkdir aucell_out

pyscenic aucell \
    input_loom_objects/${dataset}.seu.loom \
    ctx_out/${dataset}.thr_0.01.csv \
    --output aucell_out/${dataset}.aucell.thr_0.01.loom \
    --num_workers 4
