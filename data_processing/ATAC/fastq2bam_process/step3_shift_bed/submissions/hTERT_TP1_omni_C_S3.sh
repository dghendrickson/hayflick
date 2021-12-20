#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_3
#SBATCH --output=hTERT_TP1_omni_C_S3.%j.out
samtools sort -n -T hTERT_TP1_omni_C_S3_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/hTERT_TP1_omni_C_S3_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP1_omni_C_S3_SE.bed
