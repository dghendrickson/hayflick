#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_33
#SBATCH --output=PDL25_TP2_std_A_S21.%j.out
samtools sort -n -T PDL25_TP2_std_A_S21_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/PDL25_TP2_std_A_S21_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL25_TP2_std_A_S21_SE.bed
