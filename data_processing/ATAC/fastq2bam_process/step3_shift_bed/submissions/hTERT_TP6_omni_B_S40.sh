#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_19
#SBATCH --output=hTERT_TP6_omni_B_S40.%j.out
samtools sort -n -T hTERT_TP6_omni_B_S40_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/hTERT_TP6_omni_B_S40_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/hTERT_TP6_omni_B_S40_SE.bed
