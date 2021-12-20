#! /bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --job-name=step3_shift_bed_35
#SBATCH --output=PDL25_TP2_std_C_S23.%j.out
samtools sort -n -T PDL25_TP2_std_C_S23_nsorted /home/yuanh/analysis/WI38_ATAC/process/bam/PDL25_TP2_std_C_S23_R1_001.filter.rmdup.bam |bedtools bamtobed -i - | /home/yuanh/programs/myscripts/NGS_scripts/adjustBedTn5.sh - > /home/yuanh/analysis/WI38_ATAC/process/shift_bed/PDL25_TP2_std_C_S23_SE.bed
