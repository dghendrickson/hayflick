#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_49
#SBATCH --output=PDL50_TP7_omni_B_S49.%j.out
samtools sort -n -T PDL50_TP7_omni_B_S49_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/PDL50_TP7_omni_B_S49_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL50_TP7_omni_B_S49_SE.bed
