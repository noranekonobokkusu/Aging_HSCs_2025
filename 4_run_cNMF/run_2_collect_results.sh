#!/bin/bash
#$ -l h_vmem=1G
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=1:00:00
#$ -R y
#$ -cwd

# for i in `cat sample_list.txt`; do qsub -v sample=$i,lowestk=3,highestk=10 run_collect_results.sh; done

source /broad/software/scripts/useuse
use Anaconda3
source activate cnmf

base_dir=/broad/vangalenlab/safina/projects/aging_hsc/prepare_for_publication/4_run_cNMF

# Collect usage files
mkdir usages
cd usages
mkdir $sample
cd $sample

cp ${base_dir}/cNMF_out/${sample}/sample/sample.usages.k_*.consensus.txt .

for i in sample.usages.*
do
  new_name=`echo $i | sed 's/sample.usages.k_//' | sed 's/consensus.//' | sed 's/dt_//'`
  mv $i $new_name
done

cd ${base_dir}

# Collect spectra files
mkdir Z_scores
cd Z_scores
mkdir $sample
cd $sample

cp ${base_dir}/cNMF_out/${sample}/sample/sample.gene_spectra_score.k_*.dt_*.txt .

for i in sample.gene_spectra_score.k_*.dt_0_*.txt
do
  new_name=`echo $i | sed 's/.*.gene_spectra_score.k_//' | sed 's/dt_//'`
  mv $i $new_name
done

cd ${base_dir}

# Collect and renormalize TPMs
mkdir TPMs
cd TPMs
mkdir $sample
cd $sample

python ${base_dir}/renormalize_tpm.py ${base_dir}/cNMF_out/${sample}/sample
