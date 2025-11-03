#!/bin/bash
#$ -l h_vmem=24G
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_rt=5:00:00
#$ -l os=RedHat7
#$ -R y
#$ -cwd


# for i in Adelman2019 Hourigan2018 Safina2024 Zhang2022 Ainciburu2023 Li2024 Weng2024; do qsub -v dataset=$i run_ctx.sh; done

# maxvmem=<90G, cpu/wallclock=~4 cores, wallclock<2h


# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse
use Python-3.9
use Anaconda3
source activate /broad/vangalenlab/safina/tools/conda/envs/scenic


mkdir ctx_out

pyscenic ctx grnboost2_out/${dataset}.csv \
    scenic_db/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    scenic_db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
    --annotations_fname scenic_db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname input_loom_objects/${dataset}.seu.loom \
    --output ctx_out/${dataset}.thr_0.01.csv \
    --mask_dropouts \
    --auc_threshold 0.01 \
    --num_workers 4
