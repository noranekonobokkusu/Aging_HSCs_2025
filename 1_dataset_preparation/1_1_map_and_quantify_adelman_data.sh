#!/bin/bash
#$ -l h_vmem=8G
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_rt=01:00:00
#$ -cwd

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse
use R-4.1
use .homer-4.10
use .star-2.7.2a
use Samtools

cd STAR_output
rm -r $sample
mkdir $sample
cd $sample

fastq1=/broad/vangalenlab/safina/projects/aging_hsc/adelman_sc_data/SRR/${sample}_1.fastq
fastq2=/broad/vangalenlab/safina/projects/aging_hsc/adelman_sc_data/SRR/${sample}_2.fastq

STAR \
--readFilesIn $fastq1 $fastq2 \
--outFileNamePrefix ${sample}_ \
--genomeDir /broad/vangalenlab/vangalen/Genomes/GRCh38.221223/GRCh38/star \
--runThreadN 4 \
--quantMode GeneCounts \
--genomeLoad NoSharedMemory \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 20000000000

Rscript ../../quantify_bam.R ${sample}_Aligned.sortedByCoord.out.bam ${sample}_gene_counts.txt

# Individual ${sample}_gene_counts.txt should then be merged into a single file "total.gene.counts.txt"
